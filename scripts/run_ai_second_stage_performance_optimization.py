from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_ai_deepsurv_signature_gene_benchmark as base  # noqa: E402
import run_ai_signature_gene_optimization as opt  # noqa: E402


OUT = ROOT / "outputs" / "ai_second_stage_performance_optimization"
OUT.mkdir(parents=True, exist_ok=True)


def sanitize(v):
    if isinstance(v, dict):
        return {str(k): sanitize(x) for k, x in v.items()}
    if isinstance(v, list):
        return [sanitize(x) for x in v]
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, float):
        return None if not math.isfinite(v) else v
    return v


def transform_external(expr, csub, genes, mode, tcga_mean=None, tcga_sd=None):
    X = expr.loc[genes, csub["geo_accession"]].T.to_numpy(float)
    if mode == "tcga_locked":
        return base.standardize_apply(X, tcga_mean, tcga_sd)
    if mode in {"cohort_zscore", "cohort_rank_normal"}:
        return opt.cohort_transform(X, mode)
    raise ValueError(mode)


def fit_core_score(tcga_expr, tcga_clin, genes_common, stability, genes, weight_mode, transform_mode):
    idx = [genes_common.index(g) for g in genes]
    Xraw = tcga_expr.T.to_numpy(float)[:, idx]
    if transform_mode == "tcga_locked":
        X, mean, sd = base.standardize_fit(Xraw)
    else:
        X = opt.cohort_transform(Xraw, transform_mode)
        mean = sd = None
    stab = stability.set_index("gene").loc[genes]
    if weight_mode == "tcga_beta":
        w = stab["tcga_beta"].to_numpy(float)
    elif weight_mode == "direction_only":
        w = np.sign(stab["tcga_beta"].to_numpy(float))
    elif weight_mode == "external_median_logHR":
        w = np.log(stab["external_median_HR_per_SD"].to_numpy(float))
    elif weight_mode == "hybrid_tcga_external":
        w = 0.5 * stab["tcga_beta"].to_numpy(float) + 0.5 * np.log(stab["external_median_HR_per_SD"].to_numpy(float))
    else:
        raise ValueError(weight_mode)
    w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
    norm = np.sqrt(np.sum(w**2)) or 1.0
    return {"genes": genes, "weights": w / norm, "mean": mean, "sd": sd, "transform_mode": transform_mode, "weight_mode": weight_mode}


def score_core(model, expr, csub):
    X = transform_external(expr, csub, model["genes"], model["transform_mode"], model["mean"], model["sd"])
    return X @ model["weights"]


def external_results_for_score(strategy_name, score_func, ready, geo_clin, geo_exprs):
    rows = []
    for _, row in ready.iterrows():
        gse = row["GSE"]
        expr = geo_exprs.get(gse)
        if expr is None:
            continue
        csub = opt.geo_cohort_subset(geo_clin, gse, expr)
        if len(csub) < 30:
            continue
        score = score_func(gse, expr, csub)
        t = pd.to_numeric(csub["OS_time_months"], errors="coerce").to_numpy(float)
        e = pd.to_numeric(csub["OS_event"], errors="coerce").to_numpy(float)
        fit = base.cox_univariate(t, e, score, True)
        rows.append(
            {
                "strategy": strategy_name,
                "GSE": gse,
                "N": len(csub),
                "Events": int(np.nansum(e)),
                "C_index": base.c_index(t, e, score),
                "HR_per_SD": fit.get("HR", np.nan),
                "Cox_p": fit.get("p", np.nan),
            }
        )
    return pd.DataFrame(rows)


def summarize(df, meta):
    ok = df.copy()
    return {
        **meta,
        "external_tests": int(len(ok)),
        "external_samples": int(ok["N"].sum()),
        "external_events": int(ok["Events"].sum()),
        "median_external_C_index": float(ok["C_index"].median()),
        "mean_external_C_index": float(ok["C_index"].mean()),
        "min_external_C_index": float(ok["C_index"].min()),
        "max_external_C_index": float(ok["C_index"].max()),
        "direction_HR_gt1_n": int((pd.to_numeric(ok["HR_per_SD"], errors="coerce") > 1).sum()),
        "significant_p05_n": int((pd.to_numeric(ok["Cox_p"], errors="coerce") < 0.05).sum()),
        "cohorts_C_index_ge_055_n": int((pd.to_numeric(ok["C_index"], errors="coerce") >= 0.55).sum()),
    }


def z(x):
    x = np.asarray(x, dtype=float)
    sd = np.nanstd(x)
    return (x - np.nanmean(x)) / (sd if sd and np.isfinite(sd) else 1.0)


def rank01(x):
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1)
    return (ranks - 0.5) / len(x)


def main():
    _, genes_common, tcga_expr, tcga_clin, ready, geo_clin, geo_exprs = opt.prepare_data()
    stability = pd.read_csv(ROOT / "outputs" / "ai_signature_gene_optimization" / "AI_core_feature_external_stability.csv")
    X = tcga_expr.T.to_numpy(float)
    time = tcga_clin["OS_time_months"].to_numpy(float)
    event = tcga_clin["OS_event"].to_numpy(float)

    all_details = []
    rows = []
    core_models = {}
    for transform in ["tcga_locked", "cohort_zscore", "cohort_rank_normal"]:
        for k in [5, 10, 15, 20, 25, 30, 40, 50]:
            genes = stability.head(k)["gene"].tolist()
            genes = [g for g in genes if g in genes_common]
            for weight_mode in ["direction_only", "tcga_beta", "external_median_logHR", "hybrid_tcga_external"]:
                model = fit_core_score(tcga_expr, tcga_clin, genes_common, stability, genes, weight_mode, transform)
                name = f"core_{transform}_{weight_mode}_top{k}"
                core_models[name] = model
                details = external_results_for_score(name, lambda _g, expr, csub, m=model: score_core(m, expr, csub), ready, geo_clin, geo_exprs)
                all_details.append(details)
                rows.append(
                    summarize(
                        details,
                        {
                            "strategy": name,
                            "strategy_family": "core_feature_score",
                            "transform": transform,
                            "weight_mode": weight_mode,
                            "feature_k": len(genes),
                            "ensemble_alpha_model": np.nan,
                        },
                    )
                )

    # Compare against the previous best linear and deep configs, then build
    # score-level ensembles using cohort-standardized or rank-normalized scores.
    best_linear_cfg = {"preprocess": "tcga_locked", "feature_k": 150, "hidden": 0, "l2": 0.1}
    best_deep_cfg = {"preprocess": "cohort_zscore", "feature_k": 100, "hidden": 4, "l2": 0.1}
    linear_model = opt.fit_final_model(X, time, event, genes_common, best_linear_cfg)
    deep_model = opt.fit_final_model(X, time, event, genes_common, best_deep_cfg)

    def model_score(model, cfg, expr, csub):
        score, missing = opt.apply_model_to_matrix(expr, csub, model, cfg)
        if missing:
            raise ValueError(";".join(missing))
        return score

    model_specs = [
        ("optimized_linear_tcga_locked_150_l2_0p1", linear_model, best_linear_cfg),
        ("optimized_deep_cohort_zscore_100_h4_l2_0p1", deep_model, best_deep_cfg),
    ]
    base_scores_by_model = {}
    for name, model, cfg in model_specs:
        def f(_g, expr, csub, m=model, c=cfg):
            return model_score(m, c, expr, csub)
        details = external_results_for_score(name, f, ready, geo_clin, geo_exprs)
        all_details.append(details)
        rows.append(
            summarize(
                details,
                {
                    "strategy": name,
                    "strategy_family": "optimized_model",
                    "transform": cfg["preprocess"],
                    "weight_mode": "cox_partial_likelihood",
                    "feature_k": int(cfg["feature_k"]),
                    "ensemble_alpha_model": np.nan,
                },
            )
        )
        base_scores_by_model[name] = (model, cfg)

    best_core_pre = pd.DataFrame(rows).sort_values(
        ["median_external_C_index", "direction_HR_gt1_n", "significant_p05_n", "min_external_C_index"],
        ascending=[False, False, False, False],
    )
    top_core_names = best_core_pre[best_core_pre["strategy_family"].eq("core_feature_score")].head(6)["strategy"].tolist()

    for model_name, model, cfg in model_specs:
        for core_name in top_core_names:
            core_model = core_models[core_name]
            for blend in ["zscore_blend", "rank_blend"]:
                for alpha in [0.25, 0.5, 0.75]:
                    def f(_g, expr, csub, m=model, c=cfg, cm=core_model, a=alpha, b=blend):
                        s1 = model_score(m, c, expr, csub)
                        s2 = score_core(cm, expr, csub)
                        if b == "rank_blend":
                            return a * rank01(s1) + (1 - a) * rank01(s2)
                        return a * z(s1) + (1 - a) * z(s2)
                    name = f"ensemble_{blend}_alpha{alpha:g}_{model_name}__{core_name}"
                    details = external_results_for_score(name, f, ready, geo_clin, geo_exprs)
                    all_details.append(details)
                    rows.append(
                        summarize(
                            details,
                            {
                                "strategy": name,
                                "strategy_family": "score_ensemble",
                                "transform": "mixed",
                                "weight_mode": blend,
                                "feature_k": int(cfg["feature_k"]),
                                "ensemble_alpha_model": alpha,
                            },
                        )
                    )

    summary = pd.DataFrame(rows).sort_values(
        ["median_external_C_index", "direction_HR_gt1_n", "significant_p05_n", "min_external_C_index"],
        ascending=[False, False, False, False],
    )
    details = pd.concat(all_details, ignore_index=True)
    summary.to_csv(OUT / "AI_second_stage_strategy_summary.csv", index=False, encoding="utf-8-sig")
    details.to_csv(OUT / "AI_second_stage_strategy_by_cohort.csv", index=False, encoding="utf-8-sig")
    best_strategy = summary.iloc[0].to_dict()
    best_detail = details[details["strategy"].eq(best_strategy["strategy"])]
    best_detail.to_csv(OUT / "AI_second_stage_best_strategy_by_cohort.csv", index=False, encoding="utf-8-sig")

    manifest = {
        "analysis": "second-stage AI performance optimization",
        "scope": "post hoc upper-bound analysis within published-signature gene universe",
        "strategies_evaluated": int(len(summary)),
        "best_strategy": best_strategy,
        "top_10_strategies": summary.head(10).to_dict(orient="records"),
        "interpretation": "Performance optimization used external labels and is hypothesis-generating. It identifies a candidate route for improving transportability, not a locked clinical model.",
    }
    manifest = sanitize(manifest)
    (OUT / "AI_second_stage_optimization_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
