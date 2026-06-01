import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
OUT_SIG = WORK / "geo_external_validation_manuscript_signature_review.csv"
OUT_PAIR = WORK / "geo_external_validation_manuscript_pairwise_review.csv"
OUT_FLAGS = WORK / "geo_external_validation_manuscript_manual_audit_flags.csv"


def to_num(df, cols):
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def gse_set(value):
    return set(re.findall(r"GSE\d+", str(value or "").upper()))


def join_sorted(values):
    vals = sorted({str(v) for v in values if pd.notna(v) and str(v)})
    return "; ".join(vals)


def median(values):
    vals = pd.to_numeric(values, errors="coerce").dropna()
    return float(np.median(vals)) if len(vals) else np.nan


def rate(mask):
    return float(np.mean(mask)) if len(mask) else np.nan


def safe_int(value):
    if pd.isna(value):
        return 0
    return int(value)


def recommendation(row):
    tested = safe_int(row["Tested_cohorts_ok"])
    group_sig = safe_int(row["Group_sig_p05_n"])
    novel_group_sig = safe_int(row["Novel_external_group_sig_p05_n"])
    cont_sig = safe_int(row["Continuous_sig_p05_n"])
    numeric_tests = safe_int(row["Reported_numeric_cutoff_tests"])
    fallback_tests = safe_int(row["Cohort_median_fallback_tests"])
    median_c = row["Median_C_index"]
    hr_rate = row["Group_HR_gt1_rate"]
    cont_rate = row["Continuous_HR_gt1_rate"]

    if tested == 0:
        return "Not_testable_first_pass"
    if numeric_tests >= 4 and group_sig >= 2 and median_c >= 0.58 and cont_sig >= 4:
        return "Primary_numeric_cutoff_candidate"
    if numeric_tests >= 3 and group_sig == 0 and cont_sig >= 3 and median_c >= 0.58 and cont_rate >= 0.8:
        return "Continuous_score_only_numeric_cutoff_not_group_ready"
    if fallback_tests >= 5 and group_sig >= 4 and novel_group_sig >= 3 and median_c >= 0.58 and hr_rate >= 0.85 and cont_rate >= 0.85:
        return "Primary_median_fallback_candidate"
    if fallback_tests >= 4 and group_sig >= 3 and median_c >= 0.56 and hr_rate >= 0.80 and cont_rate >= 0.80:
        return "Secondary_median_fallback_candidate"
    if cont_sig >= 3 and median_c >= 0.56 and cont_rate >= 0.80:
        return "Continuous_score_support_only"
    return "Do_not_advance_first_pass"


def use_notes(row):
    notes = []
    if row["Recommendation_Tier"] == "Primary_numeric_cutoff_candidate":
        notes.append("Best fit for formula-plus-reported-cutoff manuscript analysis.")
    if row["Recommendation_Tier"] == "Primary_median_fallback_candidate":
        notes.append("Strong cross-GEO signal, but grouped survival depends on cohort median fallback.")
    if row["Recommendation_Tier"] == "Secondary_median_fallback_candidate":
        notes.append("Potential sensitivity/secondary set; review cohort roles and cutoff before main claims.")
    if row["Recommendation_Tier"] == "Continuous_score_only_numeric_cutoff_not_group_ready":
        notes.append("Continuous score is promising; reported cutoff did not separate external risk groups.")
    if row["Original_validation_GSEs_tested"]:
        notes.append("Some tested GEO cohorts were also reported as article validation cohorts.")
    if row["GSE50081_in_ok_tests"] == "yes":
        notes.append("Contains GSE50081; survival time unit needs final manuscript verification.")
    if row["Single_group_after_cutoff_rows"] > 0:
        notes.append("At least one cohort collapsed to a single risk group under reported cutoff.")
    if row["Score_contains_missing_rows"] > 0:
        notes.append("At least one cohort had missing sample-level scores despite nominal gene coverage.")
    if row["Insufficient_gene_coverage_rows"] > 0:
        notes.append("Some GEO platforms lacked full signature gene coverage.")
    return " ".join(notes)


def main():
    results = pd.read_csv(WORK / "geo_external_validation_results.csv")
    readiness = pd.read_csv(WORK / "luad_signature_readiness.csv")
    plan = pd.read_csv(WORK / "luad_external_validation_plan.csv")

    numeric_cols = [
        "Gene_count",
        "Present_gene_count",
        "Missing_gene_count",
        "Gene_coverage",
        "N",
        "Events",
        "Cutoff_used",
        "High_n",
        "Low_n",
        "High_events",
        "Low_events",
        "Score_mean",
        "Score_sd",
        "Score_median",
        "C_index",
        "Cox_continuous_HR_per_SD",
        "Cox_continuous_p",
        "Cox_high_vs_low_HR",
        "Cox_high_vs_low_p",
        "Logrank_p",
    ]
    results = to_num(results, numeric_cols)

    ready_meta = readiness.set_index("Signature_ID", drop=False).to_dict("index")
    plan_meta = plan.set_index("Signature_ID", drop=False).to_dict("index")

    pair_rows = []
    sig_rows = []
    flag_rows = []

    for sig in sorted(readiness["Signature_ID"].unique()):
        meta = ready_meta.get(sig, {})
        pmeta = plan_meta.get(sig, {})
        rows = results[results["Signature_ID"].eq(sig)].copy()
        ok = rows[rows["status"].eq("ok")].copy()
        article_validation_gses = gse_set(meta.get("Usable_article_validation_cohorts", ""))

        if len(ok):
            ok["Cohort_role_for_signature"] = ok["GSE"].apply(
                lambda g: "original_article_validation" if str(g).upper() in article_validation_gses else "novel_external"
            )
            ok["Group_p_lt_0_05"] = ok["Cox_high_vs_low_p"] < 0.05
            ok["Continuous_p_lt_0_05"] = ok["Cox_continuous_p"] < 0.05
            ok["GSE50081_time_unit_flag"] = ok["GSE"].eq("GSE50081")
            pair_rows.extend(ok.to_dict("records"))

        status_counts = rows["status"].value_counts().to_dict()
        tested = len(ok)
        original_validation_tested = sorted(set(ok["GSE"]) & article_validation_gses) if tested else []
        novel_ok = ok[ok["Cohort_role_for_signature"].eq("novel_external")] if tested else ok
        original_ok = ok[ok["Cohort_role_for_signature"].eq("original_article_validation")] if tested else ok

        row = {
            "Signature_ID": sig,
            "Article_ID": meta.get("Article_ID", ""),
            "PMID": meta.get("PMID", ""),
            "Title": meta.get("Title", ""),
            "Endpoint": meta.get("Endpoint", ""),
            "Omics_type": meta.get("Omics_type", ""),
            "Ready_for_expression_scoring": meta.get("Ready_for_expression_scoring", ""),
            "Gene_count": meta.get("Gene_count", ""),
            "Unmatched_gene_count": meta.get("Unmatched_gene_count", ""),
            "Unmatched_genes": meta.get("Unmatched_genes", ""),
            "Training_cohort": meta.get("Training_cohort", ""),
            "Validation_cohort": meta.get("Validation_cohort", ""),
            "Usable_article_validation_cohorts": meta.get("Usable_article_validation_cohorts", ""),
            "Training_validation_overlap": meta.get("Training_validation_overlap", ""),
            "Overlap_note": pmeta.get("Overlap_note", ""),
            "Cutoff_type": meta.get("Cutoff_type", ""),
            "Cutoff_value": meta.get("Cutoff_value", ""),
            "Tested_cohorts_ok": tested,
            "Tested_GEOs_ok": join_sorted(ok["GSE"]) if tested else "",
            "Total_events_across_ok_tests": int(ok["Events"].sum()) if tested else 0,
            "Group_sig_p05_n": int((ok["Cox_high_vs_low_p"] < 0.05).sum()) if tested else 0,
            "Group_sig_p05_rate": rate(ok["Cox_high_vs_low_p"] < 0.05) if tested else np.nan,
            "Continuous_sig_p05_n": int((ok["Cox_continuous_p"] < 0.05).sum()) if tested else 0,
            "Continuous_sig_p05_rate": rate(ok["Cox_continuous_p"] < 0.05) if tested else np.nan,
            "Median_C_index": median(ok["C_index"]) if tested else np.nan,
            "Best_C_index": float(ok["C_index"].max()) if tested else np.nan,
            "Median_group_HR": median(ok["Cox_high_vs_low_HR"]) if tested else np.nan,
            "Group_HR_gt1_rate": rate(ok["Cox_high_vs_low_HR"] > 1) if tested else np.nan,
            "Continuous_HR_gt1_rate": rate(ok["Cox_continuous_HR_per_SD"] > 1) if tested else np.nan,
            "Reported_numeric_cutoff_tests": int((ok["Cutoff_rule_used"] == "reported_numeric_cutoff").sum()) if tested else 0,
            "Cohort_median_fallback_tests": int((ok["Cutoff_rule_used"] == "cohort_median_fallback").sum()) if tested else 0,
            "All_group_tests_use_median_fallback": "yes"
            if tested and (ok["Cutoff_rule_used"] == "cohort_median_fallback").all()
            else "no",
            "Original_validation_GSEs_tested": "; ".join(original_validation_tested),
            "Original_validation_tests_ok_n": len(original_ok) if tested else 0,
            "Original_validation_group_sig_p05_n": int((original_ok["Cox_high_vs_low_p"] < 0.05).sum()) if tested else 0,
            "Novel_external_tests_ok_n": len(novel_ok) if tested else 0,
            "Novel_external_group_sig_p05_n": int((novel_ok["Cox_high_vs_low_p"] < 0.05).sum()) if tested else 0,
            "Novel_external_continuous_sig_p05_n": int((novel_ok["Cox_continuous_p"] < 0.05).sum()) if tested else 0,
            "Skipped_training_overlap_rows": status_counts.get("skipped_training_overlap", 0),
            "Insufficient_gene_coverage_rows": status_counts.get("insufficient_gene_coverage", 0),
            "Single_group_after_cutoff_rows": status_counts.get("single_group_after_cutoff", 0),
            "Score_contains_missing_rows": status_counts.get("score_contains_missing", 0),
            "GSE50081_in_ok_tests": "yes" if tested and ok["GSE"].eq("GSE50081").any() else "no",
        }
        row["Recommendation_Tier"] = recommendation(row)
        row["Manuscript_Use_Notes"] = use_notes(row)
        sig_rows.append(row)

        if row["Single_group_after_cutoff_rows"] > 0:
            flag_rows.append(
                {
                    "Signature_ID": sig,
                    "Flag_type": "single_group_after_reported_cutoff",
                    "Severity": "high",
                    "Details": "Reported numeric cutoff produced only one risk group in one or more GEO cohorts.",
                }
            )
        if row["Score_contains_missing_rows"] > 0:
            flag_rows.append(
                {
                    "Signature_ID": sig,
                    "Flag_type": "score_contains_missing",
                    "Severity": "medium",
                    "Details": "At least one cohort had missing scores; inspect sample-level expression before manuscript use.",
                }
            )
        if row["All_group_tests_use_median_fallback"] == "yes" and row["Group_sig_p05_n"] > 0:
            flag_rows.append(
                {
                    "Signature_ID": sig,
                    "Flag_type": "median_fallback_cutoff_dependence",
                    "Severity": "medium",
                    "Details": "All significant grouped survival results use cohort-specific median fallback cutoff.",
                }
            )
        if row["Original_validation_GSEs_tested"]:
            flag_rows.append(
                {
                    "Signature_ID": sig,
                    "Flag_type": "published_validation_cohort_reuse",
                    "Severity": "medium",
                    "Details": f"Tested cohort(s) also appeared in article validation: {row['Original_validation_GSEs_tested']}.",
                }
            )
        if row["GSE50081_in_ok_tests"] == "yes":
            flag_rows.append(
                {
                    "Signature_ID": sig,
                    "Flag_type": "GSE50081_time_unit_review",
                    "Severity": "medium",
                    "Details": "GSE50081 survival_time was interpreted as years and converted to months; verify source unit.",
                }
            )

    sig_df = pd.DataFrame(sig_rows)
    tier_order = {
        "Primary_numeric_cutoff_candidate": 1,
        "Primary_median_fallback_candidate": 2,
        "Secondary_median_fallback_candidate": 3,
        "Continuous_score_only_numeric_cutoff_not_group_ready": 4,
        "Continuous_score_support_only": 5,
        "Do_not_advance_first_pass": 6,
        "Not_testable_first_pass": 7,
    }
    sig_df["_tier_order"] = sig_df["Recommendation_Tier"].map(tier_order).fillna(99)
    sig_df = sig_df.sort_values(
        [
            "_tier_order",
            "Group_sig_p05_n",
            "Novel_external_group_sig_p05_n",
            "Continuous_sig_p05_n",
            "Median_C_index",
        ],
        ascending=[True, False, False, False, False],
    ).drop(columns=["_tier_order"])

    pair_df = pd.DataFrame(pair_rows)
    if not pair_df.empty:
        pair_df = pair_df.sort_values(["Signature_ID", "Cohort_role_for_signature", "GSE"])

    flag_df = pd.DataFrame(flag_rows).sort_values(["Severity", "Flag_type", "Signature_ID"]) if flag_rows else pd.DataFrame()

    sig_df.to_csv(OUT_SIG, index=False, encoding="utf-8")
    pair_df.to_csv(OUT_PAIR, index=False, encoding="utf-8")
    flag_df.to_csv(OUT_FLAGS, index=False, encoding="utf-8")

    print(f"Wrote {OUT_SIG} ({len(sig_df)} rows)")
    print(f"Wrote {OUT_PAIR} ({len(pair_df)} rows)")
    print(f"Wrote {OUT_FLAGS} ({len(flag_df)} rows)")


if __name__ == "__main__":
    main()
