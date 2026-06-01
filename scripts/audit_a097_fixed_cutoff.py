from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
OUT = WORK / "a097_fixed_cutoff_1p66_stress_test.csv"

SIGNATURE_ID = "A097_PMID36971680"
CUTOFF = 1.66


def main():
    scores = pd.read_csv(WORK / "geo_external_validation_sample_scores.csv")
    results = pd.read_csv(WORK / "geo_external_validation_results.csv")
    sub = scores[scores["Signature_ID"].eq(SIGNATURE_ID)].copy()
    res = results[results["Signature_ID"].eq(SIGNATURE_ID)].copy()
    if sub.empty:
        raise SystemExit(f"No sample scores found for {SIGNATURE_ID}")

    rows = []
    for gse, g in sub.groupby("GSE"):
        risk = pd.to_numeric(g["Risk_score"], errors="coerce").dropna()
        source = res[res["GSE"].eq(gse)]
        original_status = source["status"].iloc[0] if len(source) else ""
        original_rule = source["Cutoff_rule_used"].iloc[0] if len(source) and "Cutoff_rule_used" in source else ""
        high = int((risk >= CUTOFF).sum())
        low = int((risk < CUTOFF).sum())
        if high == 0:
            stress_status = "all_low_with_fixed_cutoff"
        elif low == 0:
            stress_status = "all_high_with_fixed_cutoff"
        else:
            stress_status = "two_groups_with_fixed_cutoff"
        rows.append(
            {
                "Signature_ID": SIGNATURE_ID,
                "Fixed_cutoff_from_uploaded_pdf": CUTOFF,
                "GSE": gse,
                "Original_first_pass_status": original_status,
                "Original_first_pass_cutoff_rule": original_rule,
                "N_scored": len(risk),
                "Score_min": risk.min(),
                "Score_p25": risk.quantile(0.25),
                "Score_median": risk.median(),
                "Score_p75": risk.quantile(0.75),
                "Score_max": risk.max(),
                "Low_n_score_lt_1p66": low,
                "High_n_score_ge_1p66": high,
                "Fixed_cutoff_stress_status": stress_status,
                "Interpretation": (
                    "Fixed TCGA/training cutoff 1.66 is not portable to this processed GEO scale."
                    if stress_status != "two_groups_with_fixed_cutoff"
                    else "Fixed cutoff yields two groups; survival testing would still require matching original preprocessing."
                ),
            }
        )

    out = pd.DataFrame(rows).sort_values("GSE")
    out.to_csv(OUT, index=False, encoding="utf-8")
    print(f"Wrote {OUT} ({len(out)} rows)")
    print(out[["GSE", "N_scored", "Score_min", "Score_median", "Score_max", "Low_n_score_lt_1p66", "High_n_score_ge_1p66", "Fixed_cutoff_stress_status"]].to_string(index=False))


if __name__ == "__main__":
    main()
