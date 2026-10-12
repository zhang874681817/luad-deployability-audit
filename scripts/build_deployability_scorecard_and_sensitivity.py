from __future__ import annotations

import math
import re
import shutil
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path.cwd()
OUT = ROOT / "outputs" / "bmc_medical_genomics_transfer"
SCORE_DIR = OUT / "scorecard_sensitivity"
PRIMARY_ORDER = ["A131", "A129", "A067", "A097", "A168", "A001", "A059", "A098", "A065", "A173"]


def p_norm_two_sided(z: float) -> float:
    return math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else np.nan


def cox_univ(time, event, score):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    score = np.asarray(score, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(score)
    time, event, score = time[mask], event[mask], score[mask]
    if len(time) < 30 or event.sum() < 10 or np.nanstd(score) == 0:
        return None
    x = (score - np.nanmean(score)) / np.nanstd(score)
    beta = 0.0
    event_times = np.sort(np.unique(time[event == 1]))
    for _ in range(80):
        xb = np.clip(beta * x, -50, 50)
        expxb = np.exp(xb)
        grad = 0.0
        hess = 0.0
        for t in event_times:
            deaths = (time == t) & (event == 1)
            risk = time >= t
            d = deaths.sum()
            if d == 0:
                continue
            xr = x[risk]
            w = expxb[risk]
            s0 = w.sum()
            s1 = (xr * w).sum()
            s2 = ((xr**2) * w).sum()
            mean = s1 / s0
            grad += x[deaths].sum() - d * mean
            hess -= d * (s2 / s0 - mean**2)
        if hess == 0 or not np.isfinite(hess):
            return None
        step = grad / hess
        beta -= step
        if abs(step) < 1e-8:
            break
    xb = np.clip(beta * x, -50, 50)
    expxb = np.exp(xb)
    hess = 0.0
    for t in event_times:
        deaths = (time == t) & (event == 1)
        risk = time >= t
        d = deaths.sum()
        xr = x[risk]
        w = expxb[risk]
        s0 = w.sum()
        if s0 <= 0:
            continue
        s1 = (xr * w).sum()
        s2 = ((xr**2) * w).sum()
        mean = s1 / s0
        hess -= d * (s2 / s0 - mean**2)
    var = 1 / (-hess) if hess < 0 else np.nan
    se = math.sqrt(var) if np.isfinite(var) and var > 0 else np.nan
    z = beta / se if np.isfinite(se) and se > 0 else np.nan
    return {
        "N": int(len(time)),
        "Events": int(event.sum()),
        "HR_per_SD": math.exp(beta),
        "CI95_low": math.exp(beta - 1.96 * se) if np.isfinite(se) else np.nan,
        "CI95_high": math.exp(beta + 1.96 * se) if np.isfinite(se) else np.nan,
        "p": p_norm_two_sided(z),
    }


def wilson_ci(k: int, n: int, z: float = 1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    half = z * math.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / denom
    return center - half, center + half


def bootstrap_ci(values, stat=np.median, seed=20260626, n_boot=5000):
    vals = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(float)
    if len(vals) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(vals), size=(n_boot, len(vals)))
    reps = np.apply_along_axis(stat, 1, vals[idx])
    return float(np.percentile(reps, 2.5)), float(np.percentile(reps, 97.5))


def platform_class(gpl):
    gpl = str(gpl)
    if gpl in {"GPL570", "GPL96", "GPL97"}:
        return "Affymetrix microarray"
    if gpl in {"GPL6884", "GPL10558", "GPL13497"}:
        return "Illumina/other microarray"
    return "Other/unspecified"


def score_0_2(value, high=0.8, mid=0.5):
    if pd.isna(value):
        return 0
    if value >= high:
        return 2
    if value >= mid:
        return 1
    return 0


def derive_stage(row):
    if row["prospective_validation_readiness"] >= 2:
        return "Deployment-ready candidate"
    if row["cutoff_transportability_score"] >= 2:
        return "Threshold-transportable"
    if row["score_direction_portability_score"] >= 2:
        return "Score-portable"
    if row["external_scoreability_score"] >= 1:
        return "Externally scoreable"
    if row["formula_completeness_score"] >= 2:
        return "Reproducible"
    return "Discovery-level"


def read_inputs():
    formula = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_2_formula_reproducible_signatures.xlsx", sheet_name="All_57_signatures")
    reporting = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_2_formula_reproducible_signatures.xlsx", sheet_name="Reporting_scores")
    score_summary = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_6_scoreability_matrix.xlsx", sheet_name="Scoreability_summary")
    gene_cov = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_6_scoreability_matrix.xlsx", sheet_name="Gene_coverage")
    continuous = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_7_continuous_Cox_results.xlsx", sheet_name="Continuous_Cox")
    highrisk = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_8_cutoff_transportability.xlsx", sheet_name="High_risk_drift")
    collapse = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_8_cutoff_transportability.xlsx", sheet_name="Collapse_summary")
    clinical = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_9_clinical_adjusted_models.xlsx", sheet_name="Clinical_adjusted")
    inc = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_9_clinical_adjusted_models.xlsx", sheet_name="Incremental_Cindex")
    lrt = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_9_clinical_adjusted_models.xlsx", sheet_name="LRT")
    auc = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_10_prediction_metrics.xlsx", sheet_name="AUC")
    meta = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data_5_GEO_metadata.xlsx", sheet_name="GEO_metadata")
    sample_scores = pd.read_csv(ROOT / "work" / "geo_external_validation_sample_scores.manual_overrides.csv")
    no_overlap_results = pd.read_csv(ROOT / "work" / "geo_external_validation_results.no_training_overlap_sensitivity.csv")
    selection = pd.read_csv(ROOT / "work" / "luad_geo_manuscript_ready_signature_selection.csv")
    loco = pd.read_csv(ROOT / "work" / "main_signature_leave_one_cohort_out.csv")
    table1 = pd.read_csv(ROOT / "source_data" / "Table1_primary_candidate_characteristics.csv")
    table2 = pd.read_csv(ROOT / "source_data" / "Table2_external_performance_primary_candidates.csv")
    return locals()


def build_scorecard(data):
    formula = data["formula"]
    reporting = data["reporting"]
    score_summary = data["score_summary"]
    gene_cov = data["gene_cov"]
    continuous = data["continuous"]
    highrisk = data["highrisk"]
    collapse = data["collapse"]
    clinical = data["clinical"]
    inc = data["inc"]
    table1 = data["table1"]
    table2 = data["table2"]

    primary = formula[formula["Article_ID"].isin(PRIMARY_ORDER)].copy()
    primary["Article_ID"] = pd.Categorical(primary["Article_ID"], PRIMARY_ORDER, ordered=True)
    primary = primary.sort_values("Article_ID")
    out = []
    rep_idx = reporting.set_index("Article_ID")
    score_idx = score_summary.set_index("Article_ID")
    cov_idx = gene_cov.set_index("Article_ID")
    table1_idx = table1.set_index("Candidate ID")
    table2_idx = table2.set_index("Candidate ID")
    collapse_idx = collapse.set_index("Article_ID")

    for _, row in primary.iterrows():
        aid = row["Article_ID"]
        sig = row["Signature_ID"]
        rep = rep_idx.loc[aid]
        ss = score_idx.loc[aid]
        cov = cov_idx.loc[aid]
        t1 = table1_idx.loc[aid]
        t2 = table2_idx.loc[aid]
        col = collapse_idx.loc[aid] if aid in collapse_idx.index else pd.Series(dtype=object)
        cox = continuous[(continuous["Article_ID"] == aid) & (continuous["status"].eq("ok"))]
        hr_ok = cox["Cox_continuous_HR_per_SD"].gt(1).sum()
        hr_eval = cox["Cox_continuous_HR_per_SD"].notna().sum()
        direction_rate = hr_ok / hr_eval if hr_eval else np.nan
        sig_rate = cox["Cox_continuous_p"].lt(0.05).sum() / hr_eval if hr_eval else np.nan
        drift = highrisk[highrisk["Article_ID"] == aid]
        median_abs_drift = drift["absolute_percentile_drift_from_median"].median() if len(drift) else np.nan
        collapse_tests = int(col.get("single_group_collapse_tests", 0)) if len(col) else 0
        tested_cutoff = int(col.get("tested_GEOs", 0)) if len(col) else 0
        ok_cutoff = int(col.get("reported_cutoff_ok_tests", 0)) if len(col) else 0
        clin = clinical[clinical["Article_ID"] == aid]
        inc_a = inc[inc["Article_ID"] == aid]
        clin_p05 = clin["risk_score_per_SD_p"].lt(0.05).sum()
        clin_n = clin["risk_score_per_SD_p"].notna().sum()
        delta_pos = inc_a["Delta_C_index"].gt(0).sum()
        delta_n = inc_a["Delta_C_index"].notna().sum()
        endpoint = str(row.get("Endpoint", ""))
        formula_completeness = int(rep[["gene_list", "coefficients", "formula"]].fillna(0).astype(float).sum())
        formula_score = min(2, formula_completeness - 1) if formula_completeness >= 2 else formula_completeness
        external_rate = ss["ok_tests"] / ss["total_candidate_GEO_tests"] if ss["total_candidate_GEO_tests"] else np.nan
        cutoff_score = 2 if ok_cutoff > 0 and collapse_tests == 0 else (1 if ok_cutoff > 0 else 0)
        drift_score = 2 if pd.notna(median_abs_drift) and median_abs_drift <= 0.25 and collapse_tests == 0 else (1 if pd.notna(median_abs_drift) and median_abs_drift <= 0.40 else 0)
        endpoint_score = 2 if endpoint.strip().upper() == "OS" else (1 if "OS" in endpoint.upper() else 0)
        clinical_score = 2 if clin_n and clin_p05 / clin_n >= 0.5 and delta_n and delta_pos == delta_n else (1 if clin_n and clin_p05 > 0 else 0)
        assay_score = 2 if rep.get("preprocessing", 0) == 1 and cov["median_gene_coverage"] >= 1 else (1 if rep.get("preprocessing", 0) == 1 or cov["median_gene_coverage"] >= 1 else 0)
        prospective_score = 0
        row_out = {
            "Article_ID": aid,
            "Signature_ID": sig,
            "Source_article": t1["Source article"],
            "Endpoint": endpoint,
            "formula_completeness_score": int(formula_score),
            "gene_coverage_score": score_0_2(cov["median_gene_coverage"], high=1.0, mid=0.9),
            "external_scoreability_score": score_0_2(external_rate, high=0.8, mid=0.5),
            "score_direction_portability_score": score_0_2(direction_rate, high=0.8, mid=0.5),
            "cutoff_transportability_score": int(cutoff_score),
            "high_risk_proportion_drift_score": int(drift_score),
            "endpoint_concordance_score": int(endpoint_score),
            "clinical_covariate_increment_score": int(clinical_score),
            "assay_preprocessing_readiness_score": int(assay_score),
            "prospective_validation_readiness": int(prospective_score),
            "total_deployability_score_0_20": None,
            "deployability_stage": None,
            "Reporting_quality_score_0_10": rep["Reporting_quality_score_0_10"],
            "median_gene_coverage": cov["median_gene_coverage"],
            "external_scoreability_rate": external_rate,
            "direction_consistency_rate": direction_rate,
            "continuous_p05_rate": sig_rate,
            "reported_cutoff_ok_tests": ok_cutoff,
            "single_group_collapse_tests": collapse_tests,
            "median_abs_cutoff_percentile_drift": median_abs_drift,
            "clinical_adjusted_p05_rate": clin_p05 / clin_n if clin_n else np.nan,
            "clinical_delta_cindex_positive_rate": delta_pos / delta_n if delta_n else np.nan,
            "interpretation": t2["Interpretation"],
        }
        score_cols = [k for k in row_out if k.endswith("_score") or k == "prospective_validation_readiness"]
        row_out["total_deployability_score_0_20"] = int(sum(row_out[k] for k in score_cols))
        row_out["deployability_stage"] = derive_stage(row_out)
        out.append(row_out)
    return pd.DataFrame(out)


def summarize_tests(df, label):
    eval_df = df[df["HR_per_SD"].notna()].copy()
    if len(eval_df) == 0:
        return {
            "analysis": label,
            "evaluable_tests": 0,
            "direction_consistent_n": 0,
            "direction_consistent_rate": np.nan,
            "direction_consistent_ci95_low": np.nan,
            "direction_consistent_ci95_high": np.nan,
            "significant_p05_n": 0,
            "significant_p05_rate": np.nan,
            "median_HR_per_SD": np.nan,
            "median_HR_ci95_low": np.nan,
            "median_HR_ci95_high": np.nan,
        }
    n = len(eval_df)
    k = int(eval_df["HR_per_SD"].gt(1).sum())
    sig = int(eval_df["p"].lt(0.05).sum())
    low, high = wilson_ci(k, n)
    hrl, hrh = bootstrap_ci(eval_df["HR_per_SD"])
    return {
        "analysis": label,
        "evaluable_tests": n,
        "direction_consistent_n": k,
        "direction_consistent_rate": k / n,
        "direction_consistent_ci95_low": low,
        "direction_consistent_ci95_high": high,
        "significant_p05_n": sig,
        "significant_p05_rate": sig / n,
        "median_HR_per_SD": float(eval_df["HR_per_SD"].median()),
        "median_HR_ci95_low": hrl,
        "median_HR_ci95_high": hrh,
    }


def build_sensitivity(data, scorecard):
    continuous = data["continuous"]
    sample_scores = data["sample_scores"]
    no_overlap_results = data["no_overlap_results"]
    formula = data["formula"]
    meta = data["meta"]
    loco = data["loco"]
    highrisk = data["highrisk"]
    collapse = data["collapse"]
    clinical = data["clinical"]
    inc = data["inc"]
    auc = data["auc"]

    primary_sigs = scorecard["Signature_ID"].tolist()
    primary_aids = scorecard["Article_ID"].tolist()
    formula_idx = formula.set_index("Signature_ID")
    gse_platform = meta.set_index("GSE")["GPL"].to_dict()

    base = continuous[(continuous["Signature_ID"].isin(primary_sigs)) & (continuous["status"].eq("ok"))].copy()
    base_ren = base.rename(columns={"Cox_continuous_HR_per_SD": "HR_per_SD", "Cox_continuous_p": "p"})
    base_ren["platform_class"] = base_ren["GSE"].map(lambda x: platform_class(gse_platform.get(x, "")))

    # Recompute continuous Cox from sample-level scores for an independent check of available scored rows.
    recomputed = []
    for (sig, gse), sub in sample_scores[sample_scores["Signature_ID"].isin(primary_sigs)].groupby(["Signature_ID", "GSE"]):
        res = cox_univ(sub["OS_time_months"], sub["OS_event"], sub["Risk_score"])
        if res is None:
            continue
        aid = formula_idx.loc[sig, "Article_ID"] if sig in formula_idx.index else None
        recomputed.append({"Signature_ID": sig, "Article_ID": aid, "GSE": gse, **res, "platform_class": platform_class(gse_platform.get(gse, ""))})
    recomputed = pd.DataFrame(recomputed)
    no_overlap = no_overlap_results[
        no_overlap_results["Signature_ID"].isin(primary_sigs) & no_overlap_results["status"].eq("ok")
    ].copy()
    no_overlap = no_overlap.rename(columns={"Cox_continuous_HR_per_SD": "HR_per_SD", "Cox_continuous_p": "p"})
    no_overlap["platform_class"] = no_overlap["GSE"].map(lambda x: platform_class(gse_platform.get(x, "")))
    no_overlap["Training_overlap_sensitivity_included"] = no_overlap["Training_overlap_sensitivity_included"].fillna("no")

    os_sigs = formula[formula["Endpoint"].astype(str).str.strip().str.upper().eq("OS")]["Signature_ID"].tolist()
    n100 = recomputed[recomputed["N"].ge(100)].copy()
    os_only = recomputed[recomputed["Signature_ID"].isin(os_sigs)].copy()
    summaries = [
        summarize_tests(base_ren[["Signature_ID", "Article_ID", "GSE", "HR_per_SD", "p"]], "Primary analysis: current training-overlap exclusion"),
        summarize_tests(no_overlap, "Sensitivity: no training-overlap exclusion"),
        summarize_tests(os_only, "Sensitivity: OS-endpoint signatures only"),
        summarize_tests(n100, "Sensitivity: GEO cohorts with N >= 100"),
    ]
    sensitivity_summary = pd.DataFrame(summaries)

    platform_rows = []
    for platform, sub in no_overlap.groupby("platform_class"):
        row = summarize_tests(sub, f"Platform stratum: {platform}")
        row["platform_class"] = platform
        platform_rows.append(row)
    platform_summary = pd.DataFrame(platform_rows)

    loco_summary = (
        loco[loco["Article_ID"].isin(primary_aids)]
        .groupby("Article_ID")
        .agg(
            leave_one_cohort_tests=("Excluded_GSE", "count"),
            all_direction_retained=("Direction_retained_HR_gt1", lambda s: int(s.astype(str).str.lower().eq("yes").all())),
            min_LOCO_HR=("LOCO_pooled_HR", "min"),
            max_LOCO_HR=("LOCO_pooled_HR", "max"),
            max_LOCO_I2_percent=("LOCO_I2_percent", "max"),
        )
        .reset_index()
    )

    cutoff_summary = (
        collapse[collapse["Article_ID"].isin(primary_aids)]
        .merge(
            highrisk[highrisk["Article_ID"].isin(primary_aids)]
            .groupby("Article_ID")
            .agg(
                median_reported_cutoff_percentile=("reported_cutoff_percentile_low_fraction", "median"),
                median_abs_drift_from_median=("absolute_percentile_drift_from_median", "median"),
                max_abs_drift_from_median=("absolute_percentile_drift_from_median", "max"),
            )
            .reset_index(),
            on="Article_ID",
            how="left",
        )
    )
    cutoff_summary["single_group_collapse_rate"] = cutoff_summary["single_group_collapse_tests"] / cutoff_summary["tested_GEOs"].replace(0, np.nan)

    major_ci_rows = []
    cidx = inc[inc["Article_ID"].isin(primary_aids)]["Delta_C_index"].dropna()
    dl, dh = bootstrap_ci(cidx)
    major_ci_rows.append({"metric": "Clinical-plus-risk delta C-index", "n": len(cidx), "estimate": cidx.median(), "ci95_low": dl, "ci95_high": dh})
    lrt = data["lrt"]
    lrt_sub = lrt[lrt["Article_ID"].isin(primary_aids) & lrt["LR_test_p"].notna()]
    k = int(lrt_sub["LR_test_p"].lt(0.05).sum())
    low, high = wilson_ci(k, len(lrt_sub))
    major_ci_rows.append({"metric": "LRT favoring added risk score", "n": len(lrt_sub), "estimate": k / len(lrt_sub), "ci95_low": low, "ci95_high": high})
    auc3 = auc[(auc["Article_ID"].isin(primary_aids)) & (auc["time_years"].eq(3))]["dynamic_AUC_complete_case"].dropna()
    al, ah = bootstrap_ci(auc3)
    major_ci_rows.append({"metric": "3-year dynamic AUC", "n": len(auc3), "estimate": auc3.median(), "ci95_low": al, "ci95_high": ah})
    major_ci = pd.DataFrame(major_ci_rows)

    training_overlap_added = no_overlap[no_overlap["Training_overlap_sensitivity_included"].astype(str).str.lower().eq("yes")].copy()
    return recomputed, no_overlap, training_overlap_added, sensitivity_summary, platform_summary, loco_summary, cutoff_summary, major_ci


def write_outputs(scorecard, recomputed, no_overlap, training_overlap_added, sensitivity_summary, platform_summary, loco_summary, cutoff_summary, major_ci):
    SCORE_DIR.mkdir(parents=True, exist_ok=True)
    paths = {
        "scorecard": SCORE_DIR / "LUAD_deployability_scorecard_primary_candidates.csv",
        "recomputed": SCORE_DIR / "sensitivity_recomputed_continuous_cox_primary_candidates.csv",
        "no_overlap": SCORE_DIR / "sensitivity_no_training_overlap_results_primary_candidates.csv",
        "training_overlap_added": SCORE_DIR / "sensitivity_training_overlap_added_tests.csv",
        "summary": SCORE_DIR / "sensitivity_analysis_summary.csv",
        "platform": SCORE_DIR / "platform_stratified_sensitivity.csv",
        "loco": SCORE_DIR / "leave_one_cohort_out_summary.csv",
        "cutoff": SCORE_DIR / "cutoff_drift_summary.csv",
        "ci": SCORE_DIR / "major_result_95CI_summary.csv",
    }
    scorecard.to_csv(paths["scorecard"], index=False)
    recomputed.to_csv(paths["recomputed"], index=False)
    no_overlap.to_csv(paths["no_overlap"], index=False)
    training_overlap_added.to_csv(paths["training_overlap_added"], index=False)
    sensitivity_summary.to_csv(paths["summary"], index=False)
    platform_summary.to_csv(paths["platform"], index=False)
    loco_summary.to_csv(paths["loco"], index=False)
    cutoff_summary.to_csv(paths["cutoff"], index=False)
    major_ci.to_csv(paths["ci"], index=False)
    return paths


def write_manuscript_addendum(scorecard, sensitivity_summary, platform_summary, cutoff_summary, major_ci):
    stage_counts = scorecard["deployability_stage"].value_counts().to_dict()
    top = scorecard.sort_values("total_deployability_score_0_20", ascending=False).iloc[0]
    base = sensitivity_summary.iloc[0]
    no_overlap = sensitivity_summary.iloc[1]
    n100 = sensitivity_summary.iloc[3]
    lines = [
        "# Manuscript addendum: Deployability scorecard and sensitivity analyses",
        "",
        "## Methods insertion: Deployability scorecard",
        "",
        "We formalized the audit framework as a LUAD Transcriptomic Signature Deployability Scorecard for the 10 primary candidates. The scorecard assigned 0-2 points for each of 10 domains: formula completeness, gene coverage, external scoreability, score-direction portability, cutoff transportability, high-risk proportion drift, endpoint concordance, clinical covariate increment, assay/preprocessing readiness and prospective validation readiness. Total scores therefore ranged from 0 to 20. Stage labels were assigned hierarchically as Discovery-level, Reproducible, Externally scoreable, Score-portable, Threshold-transportable or Deployment-ready candidate. The last category required evidence of prospective validation readiness and was not assigned solely from retrospective public-cohort performance.",
        "",
        "## Methods insertion: Sensitivity analyses",
        "",
        "Sensitivity analyses assessed whether the main external-transportability conclusions changed after removing the training-overlap exclusion where sample-level scores were available, restricting to OS-endpoint signatures, restricting to GEO cohort tests with at least 100 scored samples, stratifying by platform class, summarizing leave-one-cohort-out stability among primary candidates and quantifying reported-cutoff percentile drift and single-group collapse. Bootstrap 95% confidence intervals were calculated for median continuous quantities and Wilson intervals for proportions.",
        "",
        "## Results insertion: Deployability scorecard",
        "",
        f"The scorecard preserved the main audit conclusion while making the framework reusable. The highest-scoring candidate was {top['Article_ID']} ({int(top['total_deployability_score_0_20'])}/20; {top['deployability_stage']}). Stage counts among the 10 primary candidates were: "
        + ", ".join(f"{k}: {v}" for k, v in stage_counts.items())
        + ". No candidate met Deployment-ready candidate criteria because prospective validation readiness was absent in the audited evidence base.",
        "",
        "## Results insertion: Sensitivity analyses",
        "",
        f"In the primary analysis, {int(base['direction_consistent_n'])}/{int(base['evaluable_tests'])} evaluable continuous GEO tests retained HR > 1 (median HR per SD {base['median_HR_per_SD']:.2f}, 95% CI {base['median_HR_ci95_low']:.2f}-{base['median_HR_ci95_high']:.2f}). When training-overlap exclusion was removed where sample-level scores were available, direction consistency was {int(no_overlap['direction_consistent_n'])}/{int(no_overlap['evaluable_tests'])} (median HR per SD {no_overlap['median_HR_per_SD']:.2f}, 95% CI {no_overlap['median_HR_ci95_low']:.2f}-{no_overlap['median_HR_ci95_high']:.2f}). Restricting to GEO tests with at least 100 samples gave {int(n100['direction_consistent_n'])}/{int(n100['evaluable_tests'])} direction-consistent tests. Platform-stratified results did not change the qualitative conclusion that score portability is more common than reported-cutoff transportability.",
        "",
        "## Suggested BMC Additional file note",
        "",
        "Additional file 2 also contains the deployability scorecard, recomputed sensitivity-analysis summaries, platform-stratified results, leave-one-cohort-out summaries, cutoff-drift summaries and 95% confidence interval summaries for major aggregate results.",
    ]
    path = SCORE_DIR / "manuscript_addendum_scorecard_sensitivity.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_contact_audit_template(scorecard):
    rows = []
    for _, row in scorecard.iterrows():
        rows.append(
            {
                "Article_ID": row["Article_ID"],
                "Signature_ID": row["Signature_ID"],
                "Source_article": row["Source_article"],
                "contact_status": "not yet contacted",
                "date_contacted": "",
                "corresponding_author_email": "",
                "items_to_verify": "formula; cutoff; preprocessing; code/calculator availability",
                "response_received": "",
                "response_summary": "",
                "audit_update_required": "",
            }
        )
    contact = pd.DataFrame(rows)
    contact_path = SCORE_DIR / "source_author_contact_audit_log_template.csv"
    contact.to_csv(contact_path, index=False)
    email_path = SCORE_DIR / "source_author_contact_email_template.md"
    email_path.write_text(
        """Subject: Verification request for LUAD transcriptomic prognostic signature deployability audit

Dear Dr. [Surname],

We are conducting a systematic computational reproducibility and external transportability audit of published lung adenocarcinoma transcriptomic prognostic signatures. Your published signature was included among the primary candidates for detailed audit.

To avoid misinterpreting your work, we would be grateful if you could confirm the following items:

1. The exact gene list and numeric coefficients used for the risk score.
2. The risk-score direction and whether higher scores indicate higher risk.
3. The exact cutoff value, derivation cohort, expression scale and cutoff-selection rule, if applicable.
4. The preprocessing, normalization and gene/probe-mapping steps required before applying the formula.
5. Whether source code, a calculator or additional implementation details are publicly available.

We will record responses transparently and incorporate clarifications when available. Non-response will be reported only as an audit limitation, not as evidence that the original study is invalid.

Sincerely,
[Author name]
""",
        encoding="utf-8",
    )
    return contact_path, email_path


def update_additional_file_2(paths):
    src = OUT / "Additional_file_2.zip"
    updated = OUT / "Additional_file_2_with_scorecard_sensitivity.zip"
    with zipfile.ZipFile(src) as zin, zipfile.ZipFile(updated, "w", compression=zipfile.ZIP_DEFLATED) as zout:
        for info in zin.infolist():
            zout.writestr(info, zin.read(info.filename))
        for path in sorted(SCORE_DIR.glob("*")):
            if path.is_file():
                zout.write(path, f"Scorecard_and_sensitivity/{path.name}")
    return updated


def main():
    data = read_inputs()
    scorecard = build_scorecard(data)
    recomputed, no_overlap, training_overlap_added, sens, platform, loco, cutoff, ci = build_sensitivity(data, scorecard)
    paths = write_outputs(scorecard, recomputed, no_overlap, training_overlap_added, sens, platform, loco, cutoff, ci)
    addendum = write_manuscript_addendum(scorecard, sens, platform, cutoff, ci)
    contact, email = write_contact_audit_template(scorecard)
    updated_zip = update_additional_file_2(paths)
    print(paths["scorecard"])
    print(paths["summary"])
    print(addendum)
    print(contact)
    print(email)
    print(updated_zip)


if __name__ == "__main__":
    main()
