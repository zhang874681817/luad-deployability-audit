import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
BASE_SCRIPT = ROOT / "scripts" / "run_geo_external_validation.py"

OUT_RESULTS = WORK / "geo_external_validation_results.manual_overrides.csv"
OUT_SCORES = WORK / "geo_external_validation_sample_scores.manual_overrides.csv"
OUT_SUMMARY = WORK / "geo_external_validation_signature_summary.manual_overrides.csv"


spec = importlib.util.spec_from_file_location("geo_validation_base", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)


MANUAL_OVERRIDES = {
    "A010_PMID40591077": [
        ("LDHA", ["LDHA"], 0.357),
        ("PLEC", ["PLEC"], 0.201),
        ("OGDH1_manual_as_OGDH", ["OGDH", "OGDH1"], 0.142),
        ("SOD1", ["SOD1"], 0.135),
        ("RXYLT1_manual_alias_TMEM5", ["RXYLT1", "TMEM5"], 0.100),
        ("COL4A1", ["COL4A1"], 0.044),
        ("DARS2", ["DARS2"], 0.032),
        ("C1QBP", ["C1QBP"], 0.031),
        ("ACAT2", ["ACAT2"], 0.024),
        ("HLA-DRB1", ["HLA-DRB1"], -0.020),
        ("ACAT1", ["ACAT1"], -0.083),
        ("TK2", ["TK2"], -0.096),
        ("SLC25A42", ["SLC25A42"], -0.141),
        ("TANGO2", ["TANGO2"], -0.177),
    ],
    "A104_PMID34760888": [
        ("NDUFS2", ["NDUFS2"], 0.1738),
        ("ATP8A2", ["ATP8A2"], -0.4258),
        ("SCO1", ["SCO1"], 0.4404),
        ("COX14", ["COX14"], -0.7904),
        ("COA6", ["COA6"], 0.6436),
        ("RRM2B", ["RRM2B"], -0.0686),
        ("TFAM", ["TFAM"], 0.1818),
        ("DARS2", ["DARS2"], 0.2117),
        ("GARS_manual_alias_GARS1", ["GARS", "GARS1"], 0.0779),
        ("YARS2", ["YARS2"], 0.4276),
        ("EFG1_manual_as_GFM1", ["GFM1", "EFG1"], 0.3687),
        ("GFM1", ["GFM1"], 0.0001),
        ("MRPL3", ["MRPL3"], -0.8567),
        ("MRPL44", ["MRPL44"], 0.1498),
        ("ISCU", ["ISCU"], -0.8863),
        ("CABC1_manual_as_COQ8A", ["COQ8A", "CABC1"], -0.4291),
        ("HSPD1", ["HSPD1"], 0.3229),
        ("ETHE1", ["ETHE1"], 0.5369),
    ],
}


def load_signatures_with_overrides():
    sigs = base.load_signatures()
    for sig, genes in MANUAL_OVERRIDES.items():
        if sig not in sigs:
            continue
        sigs[sig]["genes"] = [
            {"display": display, "aliases": aliases, "coef": coef}
            for display, aliases, coef in genes
        ]
    for sig, obj in sigs.items():
        if sig in MANUAL_OVERRIDES:
            continue
        obj["genes"] = [
            {"display": gene, "aliases": [gene], "coef": coef}
            for gene, coef in obj["genes"]
        ]
    return sigs


def resolve_genes(genes, expr_index):
    resolved = []
    missing = []
    for gene in genes:
        hit = next((alias for alias in gene["aliases"] if alias in expr_index), None)
        if hit is None:
            missing.append(gene["display"])
        else:
            resolved.append((gene["display"], hit, gene["coef"]))
    return resolved, missing


def main():
    sigs = load_signatures_with_overrides()
    clin = pd.read_csv(base.CLIN)
    ready = pd.read_csv(base.READY)
    ready = ready[ready["Machine_ready_LUAD_like"].eq("yes") & ready["Matrix_status"].eq("ok")]
    expr_cache = {}
    results = []
    scores_out = []

    for _, cohort in ready.iterrows():
        gse = cohort["GSE"]
        expr_path = cohort["Scoring_expression_file"]
        if not expr_path or not Path(expr_path).exists():
            continue
        expr = expr_cache.setdefault(gse, base.load_expression(expr_path))
        sample_cols = set(expr.columns)
        csub = clin[clin["GSE"].eq(gse)].copy()
        csub["LUAD_like"] = csub.apply(base.is_luad_like, axis=1)
        csub = csub[
            csub["LUAD_like"]
            & csub["geo_accession"].isin(sample_cols)
            & csub["OS_time_months"].notna()
            & csub["OS_event"].notna()
        ]
        if len(csub) < 30:
            continue
        for sig, obj in sigs.items():
            meta = obj["meta"]
            readiness = obj["readiness"]
            train = base.split_cohorts(readiness.get("Training_cohort", meta.get("Training_cohort", "")))
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
                        "Override_note": "manual_formula_symbol_override" if sig in MANUAL_OVERRIDES else "",
                    }
                )
                continue
            genes = obj["genes"]
            resolved, missing = resolve_genes(genes, set(expr.index))
            coverage = len(resolved) / len(genes) if genes else 0
            base_row = {
                "Signature_ID": sig,
                "Article_ID": meta.get("Article_ID", ""),
                "PMID": meta.get("PMID", ""),
                "Title": meta.get("Title", ""),
                "GSE": gse,
                "Training_cohort": readiness.get("Training_cohort", ""),
                "Validation_cohort": readiness.get("Validation_cohort", ""),
                "Gene_count": len(genes),
                "Present_gene_count": len(resolved),
                "Missing_gene_count": len(missing),
                "Gene_coverage": round(coverage, 4),
                "Missing_genes": "; ".join(missing),
                "Resolved_gene_symbols": "; ".join(f"{display}->{symbol}" for display, symbol, _ in resolved),
                "N": len(csub),
                "Events": int(csub["OS_event"].sum()),
                "Scoring_transform": cohort.get("Scoring_transform", ""),
                "Override_note": "manual_formula_symbol_override" if sig in MANUAL_OVERRIDES else "",
            }
            if coverage < 1.0:
                base_row["status"] = "insufficient_gene_coverage"
                results.append(base_row)
                continue
            needed = [symbol for _, symbol, _ in resolved]
            xmat = expr.loc[needed, csub["geo_accession"].tolist()].astype(float)
            coefs = np.array([coef for _, _, coef in resolved], dtype=float)
            score = np.dot(coefs, xmat.to_numpy(dtype=float))
            if np.isnan(score).any():
                base_row["status"] = "score_contains_missing"
                results.append(base_row)
                continue
            cutoff_value = base.parse_cutoff(meta.get("Cutoff_value", "") or readiness.get("Cutoff_value", ""))
            cutoff_type = meta.get("Cutoff_type", "") or readiness.get("Cutoff_type", "")
            if cutoff_value is None:
                cutoff = float(np.median(score))
                cutoff_used = "cohort_median_fallback"
            else:
                cutoff = cutoff_value
                cutoff_used = "reported_numeric_cutoff"
            group = (score > cutoff).astype(int)
            if len(set(group)) < 2:
                base_row["status"] = "single_group_after_cutoff"
                base_row["Cutoff_used"] = cutoff
                base_row["Cutoff_rule_used"] = cutoff_used
                results.append(base_row)
                continue
            time = csub["OS_time_months"].astype(float).to_numpy()
            event = csub["OS_event"].astype(float).to_numpy()
            cox_cont = base.cox_univariate(time, event, score, scale_x=True)
            cox_group = base.cox_univariate(time, event, group, scale_x=False)
            lr = base.logrank(time, event, group)
            cin = base.c_index(time, event, score)
            base_row.update(
                {
                    "status": "ok",
                    "Cutoff_type_original": cutoff_type,
                    "Cutoff_used": cutoff,
                    "Cutoff_rule_used": cutoff_used,
                    "High_n": int(group.sum()),
                    "Low_n": int((group == 0).sum()),
                    "High_events": int(event[group == 1].sum()),
                    "Low_events": int(event[group == 0].sum()),
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
            results.append(base_row)
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
                        "Override_note": base_row["Override_note"],
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
        print(
            ok[ok["Signature_ID"].isin(MANUAL_OVERRIDES)]
            [["Signature_ID", "GSE", "N", "Events", "Cox_high_vs_low_HR", "Cox_high_vs_low_p", "C_index"]]
            .to_string(index=False)
        )


if __name__ == "__main__":
    main()
