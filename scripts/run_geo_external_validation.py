import csv
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
PREP = WORK / "luad_reproduction_preparation.json"
CLIN = WORK / "geo_clinical_samples.csv"
READY = WORK / "geo_expression_clinical_readiness.csv"
OUT_RESULTS = WORK / "geo_external_validation_results.csv"
OUT_SCORES = WORK / "geo_external_validation_sample_scores.csv"
OUT_SUMMARY = WORK / "geo_external_validation_signature_summary.csv"


def norm(s):
    return re.sub(r"\s+", " ", str(s or "")).strip()


def split_cohorts(v):
    out = []
    for token in re.split(r"[;/|,]", str(v or "")):
        t = norm(token).upper()
        if not t or t == "NAN":
            continue
        t = t.replace("TCGA LUAD", "TCGA-LUAD")
        out.append(t)
    return set(out)


def is_luad_like(row):
    text = " ".join(str(row.get(c, "")) for c in ["histology", "sample_type", "stage"])
    text_l = text.lower()
    if any(x in text_l for x in ["normal", "squamous", "large cell", "ntl"]):
        return False
    if any(x in text_l for x in ["adeno", "adenocarcinoma", "lung tumour", "lung tumor", "primary lung tumor", "nsclc"]):
        return True
    return False


def parse_cutoff(value):
    s = str(value or "")
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except Exception:
        return None


def p_norm_two_sided(z):
    return math.erfc(abs(z) / math.sqrt(2.0))


def p_chisq1(x):
    if x < 0 or not math.isfinite(x):
        return float("nan")
    return math.erfc(math.sqrt(x / 2.0))


def cox_univariate(time, event, x, scale_x=False):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    x = np.asarray(x, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(x)
    time, event, x = time[mask], event[mask], x[mask]
    if len(time) < 10 or event.sum() < 3 or np.nanstd(x) == 0:
        return {}
    if scale_x:
        x = (x - np.mean(x)) / np.std(x)
    beta = 0.0
    event_times = np.sort(np.unique(time[event == 1]))
    for _ in range(50):
        grad = 0.0
        hess = 0.0
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
            break
        step = grad / hess
        beta -= step
        if abs(step) < 1e-7:
            break
    # Final Hessian for SE
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
    se = math.sqrt(-1.0 / hess)
    z = beta / se if se else float("nan")
    hr = math.exp(beta)
    return {
        "beta": beta,
        "HR": hr,
        "CI95_low": math.exp(beta - 1.96 * se),
        "CI95_high": math.exp(beta + 1.96 * se),
        "p": p_norm_two_sided(z),
    }


def logrank(time, event, group):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    group = np.asarray(group, dtype=int)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(group)
    time, event, group = time[mask], event[mask], group[mask]
    if len(set(group)) < 2 or event.sum() < 3:
        return {}
    o1 = e1 = v1 = 0.0
    for t in np.sort(np.unique(time[event == 1])):
        risk = time >= t
        deaths = (time == t) & (event == 1)
        n = risk.sum()
        n1 = (risk & (group == 1)).sum()
        d = deaths.sum()
        d1 = (deaths & (group == 1)).sum()
        if n <= 1:
            continue
        exp1 = d * n1 / n
        var1 = n1 * (n - n1) * d * (n - d) / (n * n * (n - 1))
        o1 += d1
        e1 += exp1
        v1 += var1
    if v1 <= 0:
        return {}
    chisq = (o1 - e1) ** 2 / v1
    return {"logrank_chisq": chisq, "logrank_p": p_chisq1(chisq), "observed_high": o1, "expected_high": e1}


def c_index(time, event, score):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    score = np.asarray(score, dtype=float)
    n = conc = ties = 0
    for i in range(len(time)):
        if event[i] != 1:
            continue
        for j in range(len(time)):
            if time[i] >= time[j]:
                continue
            if not (np.isfinite(score[i]) and np.isfinite(score[j])):
                continue
            n += 1
            if score[i] > score[j]:
                conc += 1
            elif score[i] == score[j]:
                ties += 1
    return (conc + 0.5 * ties) / n if n else float("nan")


def load_signatures():
    payload = json.loads(PREP.read_text(encoding="utf-8"))
    by_sig = {}
    for row in payload["signature_gene_coefficients"]:
        sig = row["Signature_ID"]
        by_sig.setdefault(sig, {"genes": [], "meta": row})
        gene = row.get("Approved_symbol") or row.get("Gene_original")
        coef = row.get("Coefficient_numeric")
        if gene and coef is not None:
            by_sig[sig]["genes"].append((gene, float(coef)))
    # Pull readiness metadata with cutoff/training.
    readiness = {r["Signature_ID"]: r for r in payload["signature_readiness"]}
    for sig in by_sig:
        by_sig[sig]["readiness"] = readiness.get(sig, {})
    return by_sig


def load_expression(path):
    df = pd.read_csv(path, sep="\t")
    return df.set_index("Gene_symbol")


def main():
    sigs = load_signatures()
    clin = pd.read_csv(CLIN)
    ready = pd.read_csv(READY)
    ready = ready[ready["Machine_ready_LUAD_like"].eq("yes") & ready["Matrix_status"].eq("ok")]
    expr_cache = {}
    results = []
    scores_out = []

    for _, cohort in ready.iterrows():
        gse = cohort["GSE"]
        expr_path = cohort["Scoring_expression_file"]
        if not expr_path or not Path(expr_path).exists():
            continue
        expr = expr_cache.setdefault(gse, load_expression(expr_path))
        sample_cols = set(expr.columns)
        csub = clin[clin["GSE"].eq(gse)].copy()
        csub["LUAD_like"] = csub.apply(is_luad_like, axis=1)
        csub = csub[csub["LUAD_like"] & csub["geo_accession"].isin(sample_cols) & csub["OS_time_months"].notna() & csub["OS_event"].notna()]
        if len(csub) < 30:
            continue
        for sig, obj in sigs.items():
            meta = obj["meta"]
            readiness = obj["readiness"]
            train = split_cohorts(readiness.get("Training_cohort", meta.get("Training_cohort", "")))
            if gse.upper() in train:
                results.append(
                    {
                        "Signature_ID": sig,
                        "Article_ID": meta.get("Article_ID", ""),
                        "PMID": meta.get("PMID", ""),
                        "Title": meta.get("Title", ""),
                        "GSE": gse,
                        "status": "skipped_training_overlap",
                        "Training_cohort": readiness.get("Training_cohort", ""),
                    }
                )
                continue
            genes = obj["genes"]
            needed = [g for g, c in genes]
            present = [g for g in needed if g in expr.index]
            missing = [g for g in needed if g not in expr.index]
            coverage = len(present) / len(needed) if needed else 0
            base = {
                "Signature_ID": sig,
                "Article_ID": meta.get("Article_ID", ""),
                "PMID": meta.get("PMID", ""),
                "Title": meta.get("Title", ""),
                "GSE": gse,
                "Training_cohort": readiness.get("Training_cohort", ""),
                "Validation_cohort": readiness.get("Validation_cohort", ""),
                "Gene_count": len(needed),
                "Present_gene_count": len(present),
                "Missing_gene_count": len(missing),
                "Gene_coverage": round(coverage, 4),
                "Missing_genes": "; ".join(missing),
                "N": len(csub),
                "Events": int(csub["OS_event"].sum()),
                "Scoring_transform": cohort.get("Scoring_transform", ""),
            }
            if coverage < 1.0:
                base["status"] = "insufficient_gene_coverage"
                results.append(base)
                continue
            xmat = expr.loc[needed, csub["geo_accession"].tolist()].astype(float)
            coefs = np.array([c for g, c in genes], dtype=float)
            score = np.dot(coefs, xmat.to_numpy(dtype=float))
            if np.isnan(score).any():
                base["status"] = "score_contains_missing"
                results.append(base)
                continue
            cutoff_value = parse_cutoff(meta.get("Cutoff_value", "") or readiness.get("Cutoff_value", ""))
            cutoff_type = meta.get("Cutoff_type", "") or readiness.get("Cutoff_type", "")
            if cutoff_value is None:
                cutoff = float(np.median(score))
                cutoff_used = "cohort_median_fallback"
            else:
                cutoff = cutoff_value
                cutoff_used = "reported_numeric_cutoff"
            group = (score > cutoff).astype(int)
            if len(set(group)) < 2:
                base["status"] = "single_group_after_cutoff"
                base["Cutoff_used"] = cutoff
                base["Cutoff_rule_used"] = cutoff_used
                results.append(base)
                continue
            time = csub["OS_time_months"].astype(float).to_numpy()
            event = csub["OS_event"].astype(float).to_numpy()
            cox_cont = cox_univariate(time, event, score, scale_x=True)
            cox_group = cox_univariate(time, event, group, scale_x=False)
            lr = logrank(time, event, group)
            cin = c_index(time, event, score)
            high_n = int(group.sum())
            high_events = int(event[group == 1].sum())
            low_n = int((group == 0).sum())
            low_events = int(event[group == 0].sum())
            base.update(
                {
                    "status": "ok",
                    "Cutoff_type_original": cutoff_type,
                    "Cutoff_used": cutoff,
                    "Cutoff_rule_used": cutoff_used,
                    "High_n": high_n,
                    "Low_n": low_n,
                    "High_events": high_events,
                    "Low_events": low_events,
                    "Score_mean": float(np.mean(score)),
                    "Score_sd": float(np.std(score)),
                    "Score_median": float(np.median(score)),
                    "C_index": cin,
                    "Cox_continuous_HR_per_SD": cox_cont.get("HR", ""),
                    "Cox_continuous_CI95_low": cox_cont.get("CI95_low", ""),
                    "Cox_continuous_CI95_high": cox_cont.get("CI95_high", ""),
                    "Cox_continuous_p": cox_cont.get("p", ""),
                    "Cox_high_vs_low_HR": cox_group.get("HR", ""),
                    "Cox_high_vs_low_CI95_low": cox_group.get("CI95_low", ""),
                    "Cox_high_vs_low_CI95_high": cox_group.get("CI95_high", ""),
                    "Cox_high_vs_low_p": cox_group.get("p", ""),
                    "Logrank_p": lr.get("logrank_p", ""),
                }
            )
            results.append(base)
            for sample, t, e, s, grp in zip(csub["geo_accession"], time, event, score, group):
                scores_out.append(
                    {
                        "Signature_ID": sig,
                        "GSE": gse,
                        "geo_accession": sample,
                        "OS_time_months": t,
                        "OS_event": int(e),
                        "Risk_score": s,
                        "Risk_group": "high" if grp == 1 else "low",
                        "Cutoff_used": cutoff,
                    }
                )

    pd.DataFrame(results).to_csv(OUT_RESULTS, index=False, encoding="utf-8-sig")
    pd.DataFrame(scores_out).to_csv(OUT_SCORES, index=False, encoding="utf-8-sig")
    rdf = pd.DataFrame(results)
    ok = rdf[rdf["status"].eq("ok")].copy() if not rdf.empty else pd.DataFrame()
    if not ok.empty:
        for col in ["Cox_high_vs_low_p", "Cox_high_vs_low_HR", "C_index"]:
            ok[col] = pd.to_numeric(ok[col], errors="coerce")
        summary = (
            ok.assign(success_p05=ok["Cox_high_vs_low_p"] < 0.05)
            .groupby(["Signature_ID", "Article_ID", "PMID", "Title"], dropna=False)
            .agg(
                Tested_cohorts=("GSE", "nunique"),
                Significant_cohorts_p05=("success_p05", "sum"),
                Median_C_index=("C_index", "median"),
                Best_C_index=("C_index", "max"),
                Median_group_HR=("Cox_high_vs_low_HR", "median"),
            )
            .reset_index()
        )
    else:
        summary = pd.DataFrame()
    summary.to_csv(OUT_SUMMARY, index=False, encoding="utf-8-sig")
    print(f"wrote {OUT_RESULTS} ({len(results)} rows)")
    print(f"wrote {OUT_SCORES} ({len(scores_out)} rows)")
    print(f"wrote {OUT_SUMMARY} ({len(summary)} rows)")
    if not ok.empty:
        print(ok[["Signature_ID", "GSE", "N", "Events", "Cox_high_vs_low_HR", "Cox_high_vs_low_p", "C_index"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
