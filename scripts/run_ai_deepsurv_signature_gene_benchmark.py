from __future__ import annotations

import gzip
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
DATA = ROOT / "data"
TCGA = DATA / "raw" / "tcga_luad_xena"
OUT = ROOT / "outputs" / "ai_deepsurv_signature_benchmark"
OUT.mkdir(parents=True, exist_ok=True)

GENE_COEF = WORK / "luad_signature_gene_coefficients.csv"
HGNC = WORK / "hgnc_complete_set.txt"
STAR_TPM = TCGA / "TCGA-LUAD.star_tpm.tsv.gz"
SURV = TCGA / "LUAD_survival.txt"
GEO_READY = WORK / "geo_expression_clinical_readiness.csv"
GEO_CLIN = WORK / "geo_clinical_samples.csv"


def p_norm_two_sided(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else np.nan


def c_index(time, event, score):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    score = np.asarray(score, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(score) & (time > 0)
    time, event, score = time[mask], event[mask], score[mask]
    n = conc = ties = 0
    for i in range(len(time)):
        if event[i] != 1:
            continue
        for j in range(len(time)):
            if time[i] >= time[j]:
                continue
            n += 1
            if score[i] > score[j]:
                conc += 1
            elif score[i] == score[j]:
                ties += 1
    return (conc + 0.5 * ties) / n if n else np.nan


def cox_univariate(time, event, x, scale_x=True):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    x = np.asarray(x, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(x) & (time > 0)
    time, event, x = time[mask], event[mask], x[mask]
    if len(time) < 30 or event.sum() < 10 or np.nanstd(x) == 0:
        return {}
    if scale_x:
        x = (x - np.mean(x)) / np.std(x)
    beta = 0.0
    event_times = np.sort(np.unique(time[event == 1]))
    for _ in range(80):
        grad = hess = 0.0
        xb = beta * x
        expxb = np.exp(np.clip(xb, -50, 50))
        for t in event_times:
            deaths = (time == t) & (event == 1)
            risk = time >= t
            d = deaths.sum()
            if d == 0:
                continue
            s0 = expxb[risk].sum()
            s1 = (x[risk] * expxb[risk]).sum()
            s2 = ((x[risk] ** 2) * expxb[risk]).sum()
            if s0 <= 0:
                continue
            grad += x[deaths].sum() - d * s1 / s0
            hess -= d * (s2 / s0 - (s1 / s0) ** 2)
        if hess == 0 or not math.isfinite(hess):
            return {}
        step = grad / hess
        beta -= step
        if abs(step) < 1e-8:
            break
    xb = beta * x
    expxb = np.exp(np.clip(xb, -50, 50))
    hess = 0.0
    for t in event_times:
        deaths = (time == t) & (event == 1)
        risk = time >= t
        d = deaths.sum()
        s0 = expxb[risk].sum()
        s1 = (x[risk] * expxb[risk]).sum()
        s2 = ((x[risk] ** 2) * expxb[risk]).sum()
        if s0 > 0:
            hess -= d * (s2 / s0 - (s1 / s0) ** 2)
    if hess >= 0:
        return {}
    se = math.sqrt(-1 / hess)
    z = beta / se if se else np.nan
    return {"beta": beta, "HR": math.exp(beta), "p": p_norm_two_sided(z)}


def sample_type(sample: str) -> str:
    parts = str(sample).split("-")
    return parts[3][:2] if len(parts) >= 4 else ""


def load_hgnc_map():
    h = pd.read_csv(HGNC, sep="\t", dtype=str)
    out = {}
    for _, row in h.iterrows():
        symbol = str(row.get("symbol", "")).strip()
        ensg = str(row.get("ensembl_gene_id", "")).strip()
        if symbol and ensg and ensg != "nan":
            out[symbol] = ensg
    return out


def load_signature_gene_pool():
    coef = pd.read_csv(GENE_COEF)
    genes = sorted(
        {
            str(g).strip()
            for g in coef["Approved_symbol"].dropna()
            if str(g).strip() and str(g).strip().lower() != "nan"
        }
    )
    return genes, coef


def read_tcga_star_symbols(symbols):
    hgnc = load_hgnc_map()
    wanted_ensg = {hgnc[g]: g for g in symbols if g in hgnc}
    rows = {}
    with gzip.open(STAR_TPM, "rt", encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        tumor_idx = [i for i, s in enumerate(samples) if sample_type(s) == "01"]
        tumor_samples = [samples[i] for i in tumor_idx]
        for line in f:
            parts = line.rstrip("\n").split("\t")
            ensg = parts[0].split(".")[0]
            if ensg in wanted_ensg:
                rows[wanted_ensg[ensg]] = [parts[i + 1] for i in tumor_idx]
    expr = pd.DataFrame(rows, index=tumor_samples).T
    return expr.apply(pd.to_numeric, errors="coerce")


def load_tcga_survival(samples):
    surv = pd.read_csv(SURV, sep="\t")
    surv["patient"] = surv["sample"].astype(str).str[:12]
    surv["OS_event"] = pd.to_numeric(surv["OS"], errors="coerce")
    surv["OS_time_months"] = pd.to_numeric(surv["OS.time"], errors="coerce") / 30.4375
    sample_df = pd.DataFrame({"sample": samples, "patient": [s[:12] for s in samples]})
    out = sample_df.merge(surv[["patient", "OS_event", "OS_time_months"]], on="patient", how="inner")
    out = out.dropna(subset=["OS_event", "OS_time_months"]).drop_duplicates("sample")
    out = out[out["OS_time_months"] > 0]
    return out


def is_luad_like(row):
    text = " ".join(str(row.get(c, "")) for c in ["histology", "sample_type", "stage"]).lower()
    if any(x in text for x in ["normal", "squamous", "large cell", "ntl"]):
        return False
    if any(x in text for x in ["adeno", "adenocarcinoma", "lung tumour", "lung tumor", "primary lung tumor", "nsclc"]):
        return True
    return False


def load_geo_expression(path):
    df = pd.read_csv(path, sep="\t")
    return df.set_index("Gene_symbol").apply(pd.to_numeric, errors="coerce")


def cox_grad_eta(time, event, eta):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    eta = np.asarray(eta, dtype=float).reshape(-1)
    eta = np.clip(eta, -50, 50)
    exp_eta = np.exp(eta)
    event_idx = np.where(event == 1)[0]
    grad = np.zeros_like(eta)
    loss = 0.0
    for i in event_idx:
        risk = time >= time[i]
        denom = exp_eta[risk].sum()
        if denom <= 0:
            continue
        loss -= eta[i] - math.log(denom)
        grad[risk] += exp_eta[risk] / denom
        grad[i] -= 1.0
    denom_events = max(1, len(event_idx))
    return loss / denom_events, grad / denom_events


@dataclass
class Model:
    hidden: int
    params: dict
    mean: np.ndarray
    sd: np.ndarray
    genes: list[str]
    config: dict


def init_params(rng, p, hidden):
    if hidden <= 0:
        return {"w": rng.normal(0, 0.02, size=(p,)), "b": np.array(0.0)}
    return {
        "W1": rng.normal(0, math.sqrt(2 / max(1, p)), size=(p, hidden)),
        "b1": np.zeros(hidden),
        "W2": rng.normal(0, math.sqrt(2 / max(1, hidden)), size=(hidden,)),
        "b2": np.array(0.0),
    }


def forward(X, params, hidden):
    if hidden <= 0:
        cache = {"X": X}
        return X @ params["w"] + params["b"], cache
    z = X @ params["W1"] + params["b1"]
    h = np.tanh(z)
    eta = h @ params["W2"] + params["b2"]
    return eta, {"X": X, "z": z, "h": h}


def backward(deta, params, cache, hidden, l2):
    X = cache["X"]
    if hidden <= 0:
        return {"w": X.T @ deta + l2 * params["w"], "b": np.array(deta.sum())}
    h = cache["h"]
    dz = (deta[:, None] * params["W2"][None, :]) * (1 - h**2)
    return {
        "W2": h.T @ deta + l2 * params["W2"],
        "b2": np.array(deta.sum()),
        "W1": X.T @ dz + l2 * params["W1"],
        "b1": dz.sum(axis=0),
    }


def train_model(X, time, event, hidden, l2, lr=0.01, epochs=900, seed=1):
    rng = np.random.default_rng(seed)
    params = init_params(rng, X.shape[1], hidden)
    m = {k: np.zeros_like(v) for k, v in params.items()}
    v = {k: np.zeros_like(v) for k, v in params.items()}
    best = None
    best_loss = np.inf
    patience = 80
    stale = 0
    for ep in range(1, epochs + 1):
        eta, cache = forward(X, params, hidden)
        loss, deta = cox_grad_eta(time, event, eta)
        reg = 0.0
        for k, val in params.items():
            if k.startswith("W") or k == "w":
                reg += 0.5 * l2 * float((val * val).sum())
        grads = backward(deta, params, cache, hidden, l2)
        loss = loss + reg
        if loss < best_loss - 1e-7:
            best_loss = loss
            best = {k: v.copy() for k, v in params.items()}
            stale = 0
        else:
            stale += 1
        if stale >= patience:
            break
        b1, b2 = 0.9, 0.999
        for k in params:
            m[k] = b1 * m[k] + (1 - b1) * grads[k]
            v[k] = b2 * v[k] + (1 - b2) * (grads[k] ** 2)
            mh = m[k] / (1 - b1**ep)
            vh = v[k] / (1 - b2**ep)
            params[k] = params[k] - lr * mh / (np.sqrt(vh) + 1e-8)
    return best or params, float(best_loss)


def standardize_fit(X):
    mean = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    Xs = (np.where(np.isfinite(X), X, mean) - mean) / sd
    return Xs, mean, sd


def standardize_apply(X, mean, sd):
    return (np.where(np.isfinite(X), X, mean) - mean) / sd


def rank_features_by_tcga(X, time, event, genes):
    rows = []
    for j, g in enumerate(genes):
        fit = cox_univariate(time, event, X[:, j], scale_x=True)
        rows.append(
            {
                "gene": g,
                "univ_beta": fit.get("beta", np.nan),
                "univ_HR_per_SD": fit.get("HR", np.nan),
                "univ_p": fit.get("p", np.nan),
                "abs_logHR": abs(fit.get("beta", 0)) if fit else 0,
            }
        )
    df = pd.DataFrame(rows)
    return df.sort_values(["univ_p", "abs_logHR"], ascending=[True, False])


def deterministic_split(samples):
    vals = []
    for s in samples:
        vals.append(sum(bytearray(str(s).encode("utf-8"))) % 5)
    vals = np.array(vals)
    return vals != 0, vals == 0


def evaluate_scores(label, time, event, score):
    cidx = c_index(time, event, score)
    cox = cox_univariate(time, event, score, scale_x=True)
    return {
        "Model": label,
        "N": int(np.isfinite(score).sum()),
        "Events": int(np.nansum(event[np.isfinite(score)])),
        "C_index": cidx,
        "HR_per_SD": cox.get("HR", np.nan),
        "Cox_p": cox.get("p", np.nan),
    }


def main():
    gene_pool, coef = load_signature_gene_pool()
    tcga_expr = read_tcga_star_symbols(gene_pool)
    tcga_clin = load_tcga_survival(list(tcga_expr.columns))
    tcga_expr = tcga_expr[tcga_clin["sample"]]

    ready = pd.read_csv(GEO_READY)
    ready = ready[ready["Machine_ready_LUAD_like"].eq("yes") & ready["Matrix_status"].eq("ok")].copy()
    geo_clin = pd.read_csv(GEO_CLIN)
    geo_exprs = {}
    common = set(tcga_expr.index)
    for _, row in ready.iterrows():
        gse = row["GSE"]
        p = Path(str(row["Scoring_expression_file"]))
        if not p.exists():
            continue
        expr = load_geo_expression(p)
        geo_exprs[gse] = expr
        common &= set(expr.index)
    genes_common = sorted(common)
    tcga_expr = tcga_expr.loc[genes_common]

    X = tcga_expr.T.to_numpy(dtype=float)
    time = tcga_clin["OS_time_months"].to_numpy(float)
    event = tcga_clin["OS_event"].to_numpy(float)
    samples = tcga_clin["sample"].tolist()
    train_mask, val_mask = deterministic_split(samples)

    ranking = rank_features_by_tcga(X[train_mask], time[train_mask], event[train_mask], genes_common)
    ranking.to_csv(OUT / "AI_feature_ranking_TCGA_training_split.csv", index=False, encoding="utf-8-sig")

    configs = []
    for k in [25, 50, 100, min(150, len(genes_common))]:
        for hidden in [0, 8, 16, 32]:
            for l2 in [1e-2, 1e-3, 1e-4]:
                configs.append({"feature_k": int(k), "hidden": hidden, "l2": l2})

    tune_rows = []
    trained = []
    for cfg in configs:
        genes = ranking.head(cfg["feature_k"])["gene"].tolist()
        idx = [genes_common.index(g) for g in genes]
        Xtr_raw = X[train_mask][:, idx]
        Xva_raw = X[val_mask][:, idx]
        Xtr, mean, sd = standardize_fit(Xtr_raw)
        Xva = standardize_apply(Xva_raw, mean, sd)
        params, loss = train_model(
            Xtr,
            time[train_mask],
            event[train_mask],
            hidden=cfg["hidden"],
            l2=cfg["l2"],
            lr=0.01 if cfg["hidden"] > 0 else 0.03,
            epochs=900,
            seed=17 + cfg["feature_k"] + cfg["hidden"],
        )
        eta_tr, _ = forward(Xtr, params, cfg["hidden"])
        eta_va, _ = forward(Xva, params, cfg["hidden"])
        row = {
            **cfg,
            "train_loss": loss,
            "TCGA_train_C_index": c_index(time[train_mask], event[train_mask], eta_tr),
            "TCGA_validation_C_index": c_index(time[val_mask], event[val_mask], eta_va),
            "TCGA_validation_HR_per_SD": cox_univariate(time[val_mask], event[val_mask], eta_va, True).get("HR", np.nan),
            "TCGA_validation_p": cox_univariate(time[val_mask], event[val_mask], eta_va, True).get("p", np.nan),
        }
        tune_rows.append(row)
        trained.append((row, params, mean, sd, genes))

    tune = pd.DataFrame(tune_rows).sort_values(
        ["TCGA_validation_C_index", "TCGA_validation_HR_per_SD"], ascending=[False, False]
    )
    tune.to_csv(OUT / "AI_deepsurv_TCGA_hyperparameter_tuning.csv", index=False, encoding="utf-8-sig")
    best_deep_row = tune[tune["hidden"].gt(0)].iloc[0].to_dict()
    best_linear_row = tune[tune["hidden"].eq(0)].iloc[0].to_dict()

    def final_fit(best_row, label):
        # Re-rank on all TCGA only after hyperparameter selection. GEO remains untouched.
        full_rank = rank_features_by_tcga(X, time, event, genes_common)
        genes = full_rank.head(int(best_row["feature_k"]))["gene"].tolist()
        idx = [genes_common.index(g) for g in genes]
        Xraw = X[:, idx]
        Xs, mean, sd = standardize_fit(Xraw)
        params, loss = train_model(
            Xs,
            time,
            event,
            hidden=int(best_row["hidden"]),
            l2=float(best_row["l2"]),
            lr=0.01 if int(best_row["hidden"]) > 0 else 0.03,
            epochs=1200,
            seed=101 + int(best_row["feature_k"]) + int(best_row["hidden"]),
        )
        return Model(int(best_row["hidden"]), params, mean, sd, genes, {**best_row, "label": label, "final_loss": loss})

    deep_model = final_fit(best_deep_row, "AI_DeepSurv_signature_gene_model")
    linear_model = final_fit(best_linear_row, "Penalized_linear_Cox_signature_gene_baseline")

    model_defs = []
    eval_rows = []
    sample_rows = []
    for model in [deep_model, linear_model]:
        idx = [genes_common.index(g) for g in model.genes]
        Xtcga = standardize_apply(X[:, idx], model.mean, model.sd)
        eta, _ = forward(Xtcga, model.params, model.hidden)
        base = evaluate_scores(model.config["label"], time, event, eta)
        eval_rows.append({"Dataset": "TCGA_training_refit", "GSE": "TCGA-LUAD", **base})
        model_defs.append(
            {
                "Model": model.config["label"],
                "hidden_units": model.hidden,
                "feature_k": len(model.genes),
                "l2": model.config["l2"],
                "TCGA_validation_C_index_selected": model.config["TCGA_validation_C_index"],
                "selected_genes": "; ".join(model.genes),
                "interpretation": "Exploratory AI benchmark trained on TCGA-LUAD using published-signature gene universe; not a deployable clinical model.",
            }
        )
        for _, row in ready.iterrows():
            gse = row["GSE"]
            if gse not in geo_exprs:
                continue
            expr = geo_exprs[gse]
            csub = geo_clin[geo_clin["GSE"].eq(gse)].copy()
            csub["LUAD_like"] = csub.apply(is_luad_like, axis=1)
            sample_cols = set(expr.columns)
            csub = csub[
                csub["LUAD_like"]
                & csub["geo_accession"].isin(sample_cols)
                & csub["OS_time_months"].notna()
                & csub["OS_event"].notna()
                & (pd.to_numeric(csub["OS_time_months"], errors="coerce") > 0)
            ].copy()
            if len(csub) < 30:
                continue
            missing = [g for g in model.genes if g not in expr.index]
            if missing:
                eval_rows.append(
                    {
                        "Dataset": "GEO_external",
                        "GSE": gse,
                        "Model": model.config["label"],
                        "status": "missing_model_genes",
                        "Missing_genes": "; ".join(missing),
                    }
                )
                continue
            Xg = expr.loc[model.genes, csub["geo_accession"]].T.to_numpy(float)
            Xg = standardize_apply(Xg, model.mean, model.sd)
            score, _ = forward(Xg, model.params, model.hidden)
            ev = pd.to_numeric(csub["OS_event"], errors="coerce").to_numpy(float)
            tm = pd.to_numeric(csub["OS_time_months"], errors="coerce").to_numpy(float)
            base = evaluate_scores(model.config["label"], tm, ev, score)
            eval_rows.append(
                {
                    "Dataset": "GEO_external",
                    "GSE": gse,
                    "status": "ok",
                    "platform_class": row.get("GPL", ""),
                    "Scoring_transform": row.get("Scoring_transform", ""),
                    **base,
                }
            )
            for acc, s, t, e in zip(csub["geo_accession"], score, tm, ev):
                sample_rows.append(
                    {
                        "Dataset": "GEO_external",
                        "GSE": gse,
                        "sample": acc,
                        "Model": model.config["label"],
                        "AI_risk_score": float(s),
                        "OS_time_months": float(t),
                        "OS_event": int(e),
                    }
                )
    eval_df = pd.DataFrame(eval_rows)
    eval_df.to_csv(OUT / "AI_deepsurv_external_validation_results.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(sample_rows).to_csv(OUT / "AI_deepsurv_external_sample_scores.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(model_defs).to_csv(OUT / "AI_deepsurv_model_definitions.csv", index=False, encoding="utf-8-sig")

    ok = eval_df[(eval_df["Dataset"].eq("GEO_external")) & (eval_df["status"].fillna("ok").eq("ok"))].copy()
    summary = (
        ok.groupby("Model")
        .agg(
            external_tests=("GSE", "count"),
            external_samples=("N", "sum"),
            external_events=("Events", "sum"),
            median_external_C_index=("C_index", "median"),
            min_external_C_index=("C_index", "min"),
            max_external_C_index=("C_index", "max"),
            direction_HR_gt1_n=("HR_per_SD", lambda s: int((pd.to_numeric(s, errors="coerce") > 1).sum())),
            significant_p05_n=("Cox_p", lambda s: int((pd.to_numeric(s, errors="coerce") < 0.05).sum())),
        )
        .reset_index()
    )
    summary.to_csv(OUT / "AI_deepsurv_external_validation_summary.csv", index=False, encoding="utf-8-sig")

    manifest = {
        "analysis": "Exploratory AI/deep-learning augmentation benchmark",
        "training_source": "TCGA-LUAD GDC STAR TPM tumor samples with OS",
        "external_sources": sorted(ok["GSE"].unique().tolist()),
        "published_signature_gene_pool_n": len(gene_pool),
        "common_train_external_gene_pool_n": len(genes_common),
        "tcga_training_samples": int(len(time)),
        "tcga_training_events": int(event.sum()),
        "best_deep_config": best_deep_row,
        "best_linear_config": best_linear_row,
        "output_files": [p.name for p in OUT.glob("AI_*")],
    }
    (OUT / "AI_deepsurv_benchmark_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({"manifest": manifest, "summary": summary.to_dict(orient="records")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
