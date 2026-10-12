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


OUT = ROOT / "outputs" / "ai_signature_gene_optimization"
OUT.mkdir(parents=True, exist_ok=True)


def json_sanitize(value):
    if isinstance(value, dict):
        return {str(k): json_sanitize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [json_sanitize(v) for v in value]
    if isinstance(value, tuple):
        return [json_sanitize(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, float):
        return None if not math.isfinite(value) else value
    return value


def rank_normal_matrix(X: np.ndarray) -> np.ndarray:
    try:
        from scipy.stats import norm

        use_norm = True
    except Exception:
        norm = None
        use_norm = False
    X = np.asarray(X, dtype=float)
    Z = np.zeros_like(X, dtype=float)
    for j in range(X.shape[1]):
        x = X[:, j]
        mask = np.isfinite(x)
        z = np.zeros(len(x), dtype=float)
        if mask.sum() <= 2:
            Z[:, j] = z
            continue
        order = np.argsort(x[mask], kind="mergesort")
        ranks = np.empty(mask.sum(), dtype=float)
        ranks[order] = np.arange(1, mask.sum() + 1)
        p = (ranks - 0.5) / mask.sum()
        if use_norm:
            vals = norm.ppf(p)
        else:
            vals = np.log(p / (1 - p)) / 1.8138
        z[mask] = vals
        Z[:, j] = z
    return Z


def cohort_transform(X: np.ndarray, mode: str, fit_mean=None, fit_sd=None):
    X = np.asarray(X, dtype=float)
    if mode == "tcga_locked":
        if fit_mean is None:
            return base.standardize_fit(X)
        return base.standardize_apply(X, fit_mean, fit_sd)
    if mode == "cohort_zscore":
        mean = np.nanmean(X, axis=0)
        sd = np.nanstd(X, axis=0)
        sd[~np.isfinite(sd) | (sd == 0)] = 1.0
        return (np.where(np.isfinite(X), X, mean) - mean) / sd
    if mode == "cohort_rank_normal":
        return rank_normal_matrix(X)
    raise ValueError(mode)


def deterministic_split(samples):
    vals = np.array([sum(bytearray(str(s).encode("utf-8"))) % 5 for s in samples])
    return vals != 0, vals == 0


def prepare_data():
    gene_pool, _ = base.load_signature_gene_pool()
    tcga_expr = base.read_tcga_star_symbols(gene_pool)
    tcga_clin = base.load_tcga_survival(list(tcga_expr.columns))
    tcga_expr = tcga_expr[tcga_clin["sample"]]

    ready = pd.read_csv(base.GEO_READY)
    ready = ready[ready["Machine_ready_LUAD_like"].eq("yes") & ready["Matrix_status"].eq("ok")].copy()
    geo_clin = pd.read_csv(base.GEO_CLIN)
    geo_exprs = {}
    common = set(tcga_expr.index)
    for _, row in ready.iterrows():
        gse = row["GSE"]
        p = Path(str(row["Scoring_expression_file"]))
        if not p.exists():
            continue
        expr = base.load_geo_expression(p)
        geo_exprs[gse] = expr
        common &= set(expr.index)
    genes_common = sorted(common)
    tcga_expr = tcga_expr.loc[genes_common]
    return gene_pool, genes_common, tcga_expr, tcga_clin, ready, geo_clin, geo_exprs


def rank_features(X, time, event, genes):
    rows = []
    for j, g in enumerate(genes):
        fit = base.cox_univariate(time, event, X[:, j], scale_x=True)
        rows.append(
            {
                "gene": g,
                "tcga_beta": fit.get("beta", np.nan),
                "tcga_HR_per_SD": fit.get("HR", np.nan),
                "tcga_p": fit.get("p", np.nan),
                "abs_beta": abs(fit.get("beta", 0)) if fit else 0,
            }
        )
    return pd.DataFrame(rows).sort_values(["tcga_p", "abs_beta"], ascending=[True, False])


def train_eval_config(X_mode, time, event, samples, genes_common, cfg):
    train_mask, val_mask = deterministic_split(samples)
    # For cohort-wise transforms, apply transform at cohort level before split.
    if cfg["preprocess"] == "tcga_locked":
        ranking_raw = X_mode
        rank_df = rank_features(ranking_raw[train_mask], time[train_mask], event[train_mask], genes_common)
        genes = rank_df.head(cfg["feature_k"])["gene"].tolist()
        idx = [genes_common.index(g) for g in genes]
        Xtr, mean, sd = base.standardize_fit(X_mode[train_mask][:, idx])
        Xva = base.standardize_apply(X_mode[val_mask][:, idx], mean, sd)
        final_fit_meta = {"mean": mean, "sd": sd}
    else:
        Xt = cohort_transform(X_mode, cfg["preprocess"])
        rank_df = rank_features(Xt[train_mask], time[train_mask], event[train_mask], genes_common)
        genes = rank_df.head(cfg["feature_k"])["gene"].tolist()
        idx = [genes_common.index(g) for g in genes]
        Xtr = Xt[train_mask][:, idx]
        Xva = Xt[val_mask][:, idx]
        final_fit_meta = {"mean": None, "sd": None}
    params, loss = base.train_model(
        Xtr,
        time[train_mask],
        event[train_mask],
        hidden=cfg["hidden"],
        l2=cfg["l2"],
        lr=0.01 if cfg["hidden"] > 0 else 0.03,
        epochs=450,
        seed=700 + cfg["feature_k"] + cfg["hidden"],
    )
    eta_tr, _ = base.forward(Xtr, params, cfg["hidden"])
    eta_va, _ = base.forward(Xva, params, cfg["hidden"])
    cox_va = base.cox_univariate(time[val_mask], event[val_mask], eta_va, True)
    return {
        **cfg,
        "genes": genes,
        "train_loss": loss,
        "TCGA_train_C_index": base.c_index(time[train_mask], event[train_mask], eta_tr),
        "TCGA_validation_C_index": base.c_index(time[val_mask], event[val_mask], eta_va),
        "TCGA_validation_HR_per_SD": cox_va.get("HR", np.nan),
        "TCGA_validation_p": cox_va.get("p", np.nan),
    }


def fit_final_model(X, time, event, genes_common, cfg):
    if cfg["preprocess"] == "tcga_locked":
        rank_df = rank_features(X, time, event, genes_common)
        genes = rank_df.head(cfg["feature_k"])["gene"].tolist()
        idx = [genes_common.index(g) for g in genes]
        Xfit, mean, sd = base.standardize_fit(X[:, idx])
    else:
        Xt = cohort_transform(X, cfg["preprocess"])
        rank_df = rank_features(Xt, time, event, genes_common)
        genes = rank_df.head(cfg["feature_k"])["gene"].tolist()
        idx = [genes_common.index(g) for g in genes]
        Xfit = Xt[:, idx]
        mean = sd = None
    params, loss = base.train_model(
        Xfit,
        time,
        event,
        hidden=cfg["hidden"],
        l2=cfg["l2"],
        lr=0.01 if cfg["hidden"] > 0 else 0.03,
        epochs=650,
        seed=1700 + cfg["feature_k"] + cfg["hidden"],
    )
    return {"params": params, "genes": genes, "mean": mean, "sd": sd, "loss": loss}


def apply_model_to_matrix(expr, csub, model, cfg):
    genes = model["genes"]
    missing = [g for g in genes if g not in expr.index]
    if missing:
        return None, missing
    Xraw = expr.loc[genes, csub["geo_accession"]].T.to_numpy(float)
    if cfg["preprocess"] == "tcga_locked":
        Xg = base.standardize_apply(Xraw, model["mean"], model["sd"])
    else:
        Xg = cohort_transform(Xraw, cfg["preprocess"])
    score, _ = base.forward(Xg, model["params"], cfg["hidden"])
    return score, []


def geo_cohort_subset(geo_clin, gse, expr):
    csub = geo_clin[geo_clin["GSE"].eq(gse)].copy()
    csub["LUAD_like"] = csub.apply(base.is_luad_like, axis=1)
    sample_cols = set(expr.columns)
    csub = csub[
        csub["LUAD_like"]
        & csub["geo_accession"].isin(sample_cols)
        & csub["OS_time_months"].notna()
        & csub["OS_event"].notna()
        & (pd.to_numeric(csub["OS_time_months"], errors="coerce") > 0)
    ].copy()
    return csub


def external_eval(model, cfg, ready, geo_clin, geo_exprs):
    rows = []
    for _, row in ready.iterrows():
        gse = row["GSE"]
        if gse not in geo_exprs:
            continue
        expr = geo_exprs[gse]
        csub = geo_cohort_subset(geo_clin, gse, expr)
        if len(csub) < 30:
            continue
        score, missing = apply_model_to_matrix(expr, csub, model, cfg)
        if missing:
            rows.append({"GSE": gse, "status": "missing_model_genes", "Missing_genes": "; ".join(missing)})
            continue
        t = pd.to_numeric(csub["OS_time_months"], errors="coerce").to_numpy(float)
        e = pd.to_numeric(csub["OS_event"], errors="coerce").to_numpy(float)
        fit = base.cox_univariate(t, e, score, True)
        rows.append(
            {
                "GSE": gse,
                "status": "ok",
                "N": len(csub),
                "Events": int(np.nansum(e)),
                "C_index": base.c_index(t, e, score),
                "HR_per_SD": fit.get("HR", np.nan),
                "Cox_p": fit.get("p", np.nan),
                "platform_class": row.get("GPL", ""),
                "Scoring_transform": row.get("Scoring_transform", ""),
            }
        )
    return pd.DataFrame(rows)


def summarize_external(df):
    ok = df[df["status"].fillna("ok").eq("ok")].copy()
    return {
        "external_tests": int(len(ok)),
        "external_samples": int(ok["N"].sum()),
        "external_events": int(ok["Events"].sum()),
        "median_external_C_index": float(ok["C_index"].median()),
        "min_external_C_index": float(ok["C_index"].min()),
        "max_external_C_index": float(ok["C_index"].max()),
        "direction_HR_gt1_n": int((pd.to_numeric(ok["HR_per_SD"], errors="coerce") > 1).sum()),
        "significant_p05_n": int((pd.to_numeric(ok["Cox_p"], errors="coerce") < 0.05).sum()),
    }


def gene_external_stability(genes_common, tcga_X_mode, time, event, ready, geo_clin, geo_exprs, preprocess):
    if preprocess == "tcga_locked":
        Xtcga, mean, sd = base.standardize_fit(tcga_X_mode)
    else:
        Xtcga = cohort_transform(tcga_X_mode, preprocess)
        mean = sd = None
    tcga_rank = rank_features(Xtcga, time, event, genes_common).set_index("gene")
    rows = []
    for j, gene in enumerate(genes_common):
        vals = {
            "gene": gene,
            "tcga_beta": tcga_rank.loc[gene, "tcga_beta"],
            "tcga_HR_per_SD": tcga_rank.loc[gene, "tcga_HR_per_SD"],
            "tcga_p": tcga_rank.loc[gene, "tcga_p"],
        }
        ext_hr = []
        ext_c = []
        ext_p = []
        for _, row in ready.iterrows():
            gse = row["GSE"]
            expr = geo_exprs.get(gse)
            if expr is None or gene not in expr.index:
                continue
            csub = geo_cohort_subset(geo_clin, gse, expr)
            if len(csub) < 30:
                continue
            xraw = expr.loc[gene, csub["geo_accession"]].to_numpy(float).reshape(-1, 1)
            if preprocess == "tcga_locked":
                idx = genes_common.index(gene)
                x = base.standardize_apply(xraw, np.array([mean[idx]]), np.array([sd[idx]])).reshape(-1)
            else:
                x = cohort_transform(xraw, preprocess).reshape(-1)
            t = pd.to_numeric(csub["OS_time_months"], errors="coerce").to_numpy(float)
            e = pd.to_numeric(csub["OS_event"], errors="coerce").to_numpy(float)
            fit = base.cox_univariate(t, e, x, True)
            if fit:
                ext_hr.append(fit["HR"])
                ext_p.append(fit["p"])
                ext_c.append(base.c_index(t, e, x))
        tcga_dir = 1 if vals["tcga_HR_per_SD"] > 1 else -1
        same = sum((h > 1 and tcga_dir > 0) or (h < 1 and tcga_dir < 0) for h in ext_hr)
        vals.update(
            {
                "external_tested_cohorts": len(ext_hr),
                "external_same_direction_n": int(same),
                "external_same_direction_rate": same / len(ext_hr) if ext_hr else np.nan,
                "external_median_HR_per_SD": float(np.nanmedian(ext_hr)) if ext_hr else np.nan,
                "external_median_C_index_as_risk": float(np.nanmedian(ext_c)) if ext_c else np.nan,
                "external_p05_n": int(sum(p < 0.05 for p in ext_p)),
                "stability_score": (same if ext_hr else 0) + int(vals["tcga_p"] < 0.05) + int(sum(p < 0.05 for p in ext_p) >= 2),
            }
        )
        rows.append(vals)
    return pd.DataFrame(rows).sort_values(
        ["stability_score", "external_same_direction_rate", "external_p05_n", "tcga_p"],
        ascending=[False, False, False, True],
    )


def main():
    gene_pool, genes_common, tcga_expr, tcga_clin, ready, geo_clin, geo_exprs = prepare_data()
    X = tcga_expr.T.to_numpy(float)
    time = tcga_clin["OS_time_months"].to_numpy(float)
    event = tcga_clin["OS_event"].to_numpy(float)
    samples = tcga_clin["sample"].tolist()

    grid = []
    for preprocess in ["tcga_locked", "cohort_zscore", "cohort_rank_normal"]:
        for feature_k in [10, 25, 50, 100, min(150, len(genes_common))]:
            for hidden in [0, 4, 8, 16]:
                for l2 in [0.1, 0.01, 0.001]:
                    grid.append({"preprocess": preprocess, "feature_k": int(feature_k), "hidden": hidden, "l2": l2})

    tune_path = OUT / "AI_parameter_optimization_TCGA_internal_grid.csv"
    if tune_path.exists():
        tune = pd.read_csv(tune_path)
        tune = tune.sort_values(["TCGA_validation_C_index", "TCGA_validation_HR_per_SD"], ascending=[False, False])
    else:
        tune_rows = []
        for cfg in grid:
            row = train_eval_config(X, time, event, samples, genes_common, cfg)
            tune_rows.append({k: v for k, v in row.items() if k != "genes"})
        tune = pd.DataFrame(tune_rows).sort_values(["TCGA_validation_C_index", "TCGA_validation_HR_per_SD"], ascending=[False, False])
        tune.to_csv(tune_path, index=False, encoding="utf-8-sig")

    # Refit and evaluate selected configs: best by TCGA validation for each
    # preprocess/model family plus a compact external landscape over all configs.
    external_rows = []
    selected = []
    for preprocess in tune["preprocess"].unique():
        sub = tune[tune["preprocess"].eq(preprocess)]
        for family, cond in [("linear", sub["hidden"].eq(0)), ("deep", sub["hidden"].gt(0))]:
            pick = sub[cond].iloc[0].to_dict()
            pick["selection_rule"] = f"best_TCGA_validation_{preprocess}_{family}"
            selected.append(pick)

    # External optimization landscape is explicitly post hoc and reported as such.
    # Evaluate a compact candidate set rather than every low-performing grid
    # member: top TCGA candidates overall, top linear candidates, and the top
    # three within each preprocessing/model-family stratum.
    tune_candidates = tune.copy()
    tune_candidates["family"] = np.where(tune_candidates["hidden"].astype(int).eq(0), "linear", "deep")
    candidate_tune = pd.concat(
        [
            tune_candidates.head(24),
            tune_candidates[tune_candidates["hidden"].astype(int).eq(0)].head(12),
            tune_candidates.groupby(["preprocess", "family"], group_keys=False).head(3),
        ],
        ignore_index=True,
    ).drop_duplicates(["preprocess", "feature_k", "hidden", "l2"])
    for _, pick_row in candidate_tune.iterrows():
        cfg = {
            "preprocess": pick_row["preprocess"],
            "feature_k": int(pick_row["feature_k"]),
            "hidden": int(pick_row["hidden"]),
            "l2": float(pick_row["l2"]),
        }
        model = fit_final_model(X, time, event, genes_common, cfg)
        ext = external_eval(model, cfg, ready, geo_clin, geo_exprs)
        summ = summarize_external(ext)
        external_rows.append({**cfg, **summ, "TCGA_validation_C_index": pick_row["TCGA_validation_C_index"]})
    landscape = pd.DataFrame(external_rows).sort_values(
        ["median_external_C_index", "direction_HR_gt1_n", "TCGA_validation_C_index"],
        ascending=[False, False, False],
    )
    landscape.to_csv(OUT / "AI_external_transportability_optimization_landscape.csv", index=False, encoding="utf-8-sig")

    selected_rows = []
    selected_details = []
    # Add best external post-hoc configs by model family to show a realistic ceiling.
    for family, cond in [("linear", landscape["hidden"].eq(0)), ("deep", landscape["hidden"].gt(0))]:
        row = landscape[cond].iloc[0].to_dict()
        row["selection_rule"] = f"posthoc_best_external_{family}"
        selected.append(row)
    seen = set()
    for pick in selected:
        cfg = {
            "preprocess": pick["preprocess"],
            "feature_k": int(pick["feature_k"]),
            "hidden": int(pick["hidden"]),
            "l2": float(pick["l2"]),
        }
        key = (cfg["preprocess"], cfg["feature_k"], cfg["hidden"], cfg["l2"], pick["selection_rule"])
        if key in seen:
            continue
        seen.add(key)
        model = fit_final_model(X, time, event, genes_common, cfg)
        ext = external_eval(model, cfg, ready, geo_clin, geo_exprs)
        summ = summarize_external(ext)
        selected_rows.append({**cfg, "selection_rule": pick["selection_rule"], **summ, "selected_genes": "; ".join(model["genes"])})
        ext.insert(0, "selection_rule", pick["selection_rule"])
        ext.insert(1, "preprocess", cfg["preprocess"])
        ext.insert(2, "feature_k", cfg["feature_k"])
        ext.insert(3, "hidden", cfg["hidden"])
        ext.insert(4, "l2", cfg["l2"])
        selected_details.append(ext)
    selected_df = pd.DataFrame(selected_rows).sort_values(["median_external_C_index", "selection_rule"], ascending=[False, True])
    selected_df.to_csv(OUT / "AI_selected_strategy_external_results_summary.csv", index=False, encoding="utf-8-sig")
    pd.concat(selected_details, ignore_index=True).to_csv(OUT / "AI_selected_strategy_external_results_by_cohort.csv", index=False, encoding="utf-8-sig")

    best_preprocess = str(landscape.iloc[0]["preprocess"])
    stability = gene_external_stability(genes_common, X, time, event, ready, geo_clin, geo_exprs, best_preprocess)
    stability.to_csv(OUT / "AI_core_feature_external_stability.csv", index=False, encoding="utf-8-sig")
    core_genes = stability.head(20)["gene"].tolist()

    manifest = {
        "analysis": "AI parameter optimization and core-feature stability analysis",
        "published_signature_gene_pool_n": len(gene_pool),
        "common_tcga_geo_gene_pool_n": len(genes_common),
        "tcga_samples": int(len(time)),
        "tcga_events": int(np.nansum(event)),
        "tcga_parameter_grid_n": int(len(tune)),
        "external_candidate_strategies_n": int(len(candidate_tune)),
        "best_external_strategy": landscape.iloc[0].to_dict(),
        "best_tcga_validation_strategy": tune.iloc[0].to_dict(),
        "top_20_core_genes": core_genes,
        "interpretation": "External optimization is post hoc and is used to identify transportability bottlenecks and candidate parameters, not to claim a locked deployable model.",
    }
    manifest = json_sanitize(manifest)
    (OUT / "AI_optimization_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
