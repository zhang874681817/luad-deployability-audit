from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
OUT = WORK / "luad_geo_manuscript_ready_signature_selection.csv"


STRICT_PRIMARY = {"A131_PMID39735553"}

PRIMARY_MEDIAN = {
    "A067_PMID40211263",
    "A129_PMID37583701",
    "A097_PMID36971680",
    "A168_PMID35711827",
    "A098_PMID32592319",
    "A065_PMID33858449",
    "A001_PMID34499810",
    "A173_PMID34434660",
    "A059_PMID38965505",
}

SECONDARY_MEDIAN = {
    "A153_PMID37143472",
    "A074_PMID41822324",
    "A130_PMID37492563",
    "A086_PMID38854940",
    "A110_PMID35096011",
    "A164_PMID34517518",
    "A112_PMID36313444",
    "A154_PMID38254708",
    "A138_PMID32850406",
    "A178_PMID33221751",
}

CONTINUOUS_SUPPORT = {
    "A005_PMID37458771",
    "A030_PMID39953119",
    "A069_PMID38111060",
    "A148_PMID36110953",
    "A062_PMID35840978",
}

MANUAL_OVERRIDE_SENSITIVITY = {
    "A010_PMID40591077",
    "A104_PMID34760888",
}


def label(row):
    sig = row["Signature_ID"]
    if sig in STRICT_PRIMARY:
        return "Tier_1_strict_reported_cutoff"
    if sig in PRIMARY_MEDIAN:
        return "Tier_2_primary_median_fallback"
    if sig in SECONDARY_MEDIAN:
        return "Tier_3_secondary_median_fallback"
    if sig in CONTINUOUS_SUPPORT:
        return "Tier_4_continuous_score_support"
    if sig in MANUAL_OVERRIDE_SENSITIVITY:
        return "Tier_5_manual_override_sensitivity"
    if row.get("Tested_cohorts_ok", 0) and row.get("Tested_cohorts_ok", 0) > 0:
        return "Do_not_advance_first_pass"
    return "Not_testable_first_pass"


def role(row):
    tier = row["Final_analysis_tier"]
    if tier == "Tier_1_strict_reported_cutoff":
        return "Main strict grouped-risk analysis using reported numeric cutoff."
    if tier == "Tier_2_primary_median_fallback":
        return "Main score-portability analysis using cohort median fallback; not original-cutoff validation."
    if tier == "Tier_3_secondary_median_fallback":
        return "Secondary/sensitivity score-portability analysis using cohort median fallback."
    if tier == "Tier_4_continuous_score_support":
        return "Continuous-risk-score support only; grouped result not robust enough for main analysis."
    if tier == "Tier_5_manual_override_sensitivity":
        return "Sensitivity-only after manual formula/symbol correction."
    if tier == "Do_not_advance_first_pass":
        return "Do not advance based on first-pass GEO evidence."
    return "Not externally testable in first-pass GEO screen."


def caveats(row):
    sig = row["Signature_ID"]
    notes = []
    def has_text(value):
        return pd.notna(value) and str(value).strip() != ""

    if row.get("All_group_tests_use_median_fallback") == "yes":
        notes.append("Grouped results use cohort-specific median fallback.")
    if has_text(row.get("Original_validation_GSEs_tested")):
        notes.append("Some GEO tests reuse article validation cohort(s).")
    if row.get("GSE50081_in_ok_tests") == "yes":
        notes.append("GSE50081 survival time verified as years and converted to months.")
    if sig == "A097_PMID36971680":
        notes.append("PDF confirms TCGA cutoff 1.66, but fixed cutoff collapses GEO cohorts to one group.")
    if sig == "A005_PMID37458771":
        notes.append("Reported cutoff available, but grouped GEO validation is weak; continuous score performs better.")
    if sig == "A010_PMID40591077":
        notes.append("Requires 14-gene formula override and OGDH/TMEM5 aliases.")
    if sig == "A104_PMID34760888":
        notes.append("Requires GARS/GARS1, EFG1/GFM1, and CABC1/COQ8A aliases.")
    return " ".join(notes)


def main():
    review = pd.read_csv(WORK / "geo_external_validation_manuscript_signature_review.manual_decisions.csv")
    override = pd.read_csv(WORK / "geo_external_validation_manual_override_comparison.csv")
    override = override[override["Manual_override_applied"].eq("yes")].copy()

    for sig in MANUAL_OVERRIDE_SENSITIVITY:
        if sig in set(review["Signature_ID"]):
            continue
        # Defensive fallback; current review has A010/A104 as not-testable rows.
        meta = override[override["Signature_ID"].eq(sig)].iloc[0].to_dict()
        review = pd.concat([review, pd.DataFrame([meta])], ignore_index=True)

    override_map = override.set_index("Signature_ID").to_dict("index")
    for col in [
        "Tested_cohorts_override",
        "Significant_cohorts_p05_override",
        "Median_C_index_override",
        "Median_group_HR_override",
    ]:
        review[col] = pd.Series([pd.NA] * len(review), dtype="object")

    for idx, row in review.iterrows():
        sig = row["Signature_ID"]
        if sig in override_map:
            for col in [
                "Tested_cohorts_override",
                "Significant_cohorts_p05_override",
                "Median_C_index_override",
                "Median_group_HR_override",
            ]:
                review.loc[idx, col] = override_map[sig].get(col, "")

    review["Final_analysis_tier"] = review.apply(label, axis=1)
    review["Recommended_manuscript_role"] = review.apply(role, axis=1)
    review["Key_caveats_for_manuscript"] = review.apply(caveats, axis=1)

    tier_order = {
        "Tier_1_strict_reported_cutoff": 1,
        "Tier_2_primary_median_fallback": 2,
        "Tier_3_secondary_median_fallback": 3,
        "Tier_4_continuous_score_support": 4,
        "Tier_5_manual_override_sensitivity": 5,
        "Do_not_advance_first_pass": 6,
        "Not_testable_first_pass": 7,
    }
    review["_tier_order"] = review["Final_analysis_tier"].map(tier_order).fillna(99)
    for col in ["Group_sig_p05_n", "Novel_external_group_sig_p05_n", "Continuous_sig_p05_n", "Median_C_index"]:
        review[col] = pd.to_numeric(review[col], errors="coerce")

    ordered_cols = [
        "Final_analysis_tier",
        "Recommended_manuscript_role",
        "Signature_ID",
        "Article_ID",
        "PMID",
        "Title",
        "Endpoint",
        "Omics_type",
        "Cutoff_type",
        "Cutoff_value",
        "Tested_cohorts_ok",
        "Group_sig_p05_n",
        "Novel_external_group_sig_p05_n",
        "Continuous_sig_p05_n",
        "Median_C_index",
        "Median_group_HR",
        "Reported_numeric_cutoff_tests",
        "Cohort_median_fallback_tests",
        "All_group_tests_use_median_fallback",
        "Original_validation_GSEs_tested",
        "GSE50081_in_ok_tests",
        "Tested_cohorts_override",
        "Significant_cohorts_p05_override",
        "Median_C_index_override",
        "Median_group_HR_override",
        "Manual_decision",
        "Manual_decision_reason",
        "Key_caveats_for_manuscript",
        "Training_cohort",
        "Validation_cohort",
    ]

    out = (
        review.sort_values(
            ["_tier_order", "Group_sig_p05_n", "Novel_external_group_sig_p05_n", "Continuous_sig_p05_n", "Median_C_index"],
            ascending=[True, False, False, False, False],
        )
        .drop(columns=["_tier_order"])
    )
    out = out[[c for c in ordered_cols if c in out.columns] + [c for c in out.columns if c not in ordered_cols]]
    out.to_csv(OUT, index=False, encoding="utf-8")
    print(f"Wrote {OUT} ({len(out)} rows)")
    print(out.groupby("Final_analysis_tier").size().to_string())


if __name__ == "__main__":
    main()
