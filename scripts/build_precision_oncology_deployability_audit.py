import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
OUT_BASE = ROOT / "outputs" / "luad_signature_extraction"
OUT = OUT_BASE / "precision_oncology_deployability_audit"
FIG = OUT / "figure_data"

PRIMARY_ARTICLES = ["A131", "A129", "A067", "A097", "A168", "A001", "A059", "A098", "A065", "A173"]
STRONG_SCORE_PORTABILITY = {"A129", "A067", "A097", "A168", "A001", "A059"}
MODERATE_SCORE_PORTABILITY = {"A098"}
LOWER_PRIORITY_SUPPORT = {"A065", "A173"}
TIMEPOINTS_AUC = [12, 36, 60]
TIMEPOINTS_PRED = [36, 60]


def read_csv(name, **kwargs):
    return pd.read_csv(WORK / name, **kwargs)


def clean_text(value):
    if pd.isna(value):
        return ""
    return str(value).strip()


def yn(value):
    return str(value).strip().lower() == "yes"


def has_value(value):
    if pd.isna(value):
        return False
    text = str(value).strip()
    return bool(text) and text.lower() not in {"nan", "none", "na", "n/a", "unclear"}


def split_tokens(value):
    tokens = []
    for token in re.split(r"[;/|,]", clean_text(value)):
        token = re.sub(r"\s+", " ", token).strip()
        if token and token.lower() not in {"nan", "none", "na", "n/a"}:
            tokens.append(token)
    return tokens


def normalize_dataset(value):
    value = re.sub(r"\s+", " ", clean_text(value)).strip()
    if not value:
        return ""
    upper = value.upper().replace("TCGA LUAD", "TCGA-LUAD")
    upper = upper.replace("TCGA_LUAD", "TCGA-LUAD")
    upper = upper.replace("INTERNAL COHORT", "INTERNAL")
    return upper


def safe_num(series):
    return pd.to_numeric(series, errors="coerce")


def p_norm_two_sided(z):
    if not np.isfinite(z):
        return np.nan
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def p_chisq1(lr):
    if not np.isfinite(lr) or lr < 0:
        return np.nan
    return math.erfc(math.sqrt(float(lr) / 2.0))


def c_index(time, event, score):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    score = np.asarray(score, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(score)
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


def cox_ph(time, event, x, names):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(x).all(axis=1)
    time, event, x = time[mask], event[mask], x[mask]
    if len(time) < 30 or event.sum() < 10 or x.shape[1] == 0:
        return {}
    keep = [i for i in range(x.shape[1]) if np.nanstd(x[:, i]) > 0]
    if not keep:
        return {}
    x = x[:, keep]
    names = [names[i] for i in keep]
    p = x.shape[1]
    beta = np.zeros(p)
    event_times = np.sort(np.unique(time[event == 1]))
    for _ in range(80):
        grad = np.zeros(p)
        hess = np.zeros((p, p))
        xb = np.clip(x @ beta, -50, 50)
        expxb = np.exp(xb)
        for t in event_times:
            deaths = (time == t) & (event == 1)
            risk = time >= t
            d = deaths.sum()
            if d == 0:
                continue
            xr = x[risk]
            w = expxb[risk]
            s0 = w.sum()
            if s0 <= 0:
                continue
            s1 = (xr * w[:, None]).sum(axis=0)
            s2 = (xr[:, :, None] * xr[:, None, :] * w[:, None, None]).sum(axis=0)
            mean = s1 / s0
            grad += x[deaths].sum(axis=0) - d * mean
            hess -= d * (s2 / s0 - np.outer(mean, mean))
        try:
            step = np.linalg.solve(hess, grad)
        except np.linalg.LinAlgError:
            return {}
        beta -= step
        if np.max(np.abs(step)) < 1e-7:
            break
    loglik = cox_loglik(time, event, x, beta)
    xb = np.clip(x @ beta, -50, 50)
    expxb = np.exp(xb)
    hess = np.zeros((p, p))
    for t in event_times:
        deaths = (time == t) & (event == 1)
        risk = time >= t
        d = deaths.sum()
        xr = x[risk]
        w = expxb[risk]
        s0 = w.sum()
        if s0 <= 0:
            continue
        s1 = (xr * w[:, None]).sum(axis=0)
        s2 = (xr[:, :, None] * xr[:, None, :] * w[:, None, None]).sum(axis=0)
        mean = s1 / s0
        hess -= d * (s2 / s0 - np.outer(mean, mean))
    try:
        cov = np.linalg.inv(-hess)
    except np.linalg.LinAlgError:
        return {}
    out = {
        "N_model": int(len(time)),
        "Events_model": int(event.sum()),
        "loglik": float(loglik),
        "beta_vector": beta,
        "used_names": names,
        "linear_predictor": x @ beta,
    }
    for i, name in enumerate(names):
        se = math.sqrt(cov[i, i]) if cov[i, i] > 0 else np.nan
        z = beta[i] / se if np.isfinite(se) and se > 0 else np.nan
        out[f"{name}_beta"] = float(beta[i])
        out[f"{name}_HR"] = float(math.exp(beta[i]))
        out[f"{name}_CI95_low"] = float(math.exp(beta[i] - 1.96 * se)) if np.isfinite(se) else np.nan
        out[f"{name}_CI95_high"] = float(math.exp(beta[i] + 1.96 * se)) if np.isfinite(se) else np.nan
        out[f"{name}_p"] = p_norm_two_sided(z)
    return out


def cox_loglik(time, event, x, beta):
    xb = np.clip(np.asarray(x, dtype=float) @ np.asarray(beta, dtype=float), -50, 50)
    expxb = np.exp(xb)
    total = 0.0
    for t in np.sort(np.unique(time[event == 1])):
        deaths = (time == t) & (event == 1)
        risk = time >= t
        d = deaths.sum()
        if d:
            total += xb[deaths].sum() - d * math.log(expxb[risk].sum())
    return total


def zscore(values):
    v = pd.to_numeric(pd.Series(values), errors="coerce")
    sd = v.std(skipna=True)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series([np.nan] * len(v), index=v.index)
    return (v - v.mean(skipna=True)) / sd


def stage_major(value):
    s = clean_text(value).upper()
    if not s:
        return np.nan
    if re.search(r"\bIV\b|^4", s):
        return 4.0
    if re.search(r"\bIII\b|^3", s):
        return 3.0
    if re.search(r"\bII\b|^2", s):
        return 2.0
    if re.search(r"\bI\b|^1", s):
        return 1.0
    return np.nan


def sex_binary(value):
    s = clean_text(value).upper()
    if s in {"M", "MALE", "1"}:
        return 1.0
    if s in {"F", "FEMALE", "0"}:
        return 0.0
    return np.nan


def dynamic_auc(time, event, score, horizon):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    score = np.asarray(score, dtype=float)
    cases = (time <= horizon) & (event == 1)
    controls = time > horizon
    valid = np.isfinite(time) & np.isfinite(event) & np.isfinite(score) & (cases | controls)
    cases = cases & valid
    controls = controls & valid
    n_cases = int(cases.sum())
    n_controls = int(controls.sum())
    if n_cases < 5 or n_controls < 5:
        return np.nan, n_cases, n_controls
    case_scores = score[cases]
    control_scores = score[controls]
    total = conc = ties = 0
    for a in case_scores:
        for b in control_scores:
            total += 1
            if a > b:
                conc += 1
            elif a == b:
                ties += 1
    return (conc + 0.5 * ties) / total, n_cases, n_controls


def km_survival_at(time, event, horizon):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event)
    time, event = time[mask], event[mask]
    if len(time) == 0:
        return np.nan
    surv = 1.0
    for t in np.sort(np.unique(time[(event == 1) & (time <= horizon)])):
        at_risk = (time >= t).sum()
        deaths = ((time == t) & (event == 1)).sum()
        if at_risk > 0:
            surv *= 1.0 - deaths / at_risk
    return surv


def censoring_km_values(time, event):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    censor_event = 1 - event
    surv = 1.0
    rows = []
    for t in np.sort(np.unique(time[censor_event == 1])):
        at_risk = (time >= t).sum()
        cens = ((time == t) & (censor_event == 1)).sum()
        if at_risk > 0:
            before = surv
            surv *= 1.0 - cens / at_risk
            rows.append((t, before, surv))
    return rows


def censoring_surv(rows, t, before=False):
    surv = 1.0
    for event_t, surv_before, surv_after in rows:
        if before:
            if event_t >= t:
                break
            surv = surv_after
        else:
            if event_t > t:
                break
            surv = surv_after
    return max(surv, 1e-6)


def baseline_cumhaz_at(time, event, x, beta, horizon):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    x = np.asarray(x, dtype=float)
    xb = np.clip(x * beta, -50, 50)
    expxb = np.exp(xb)
    h0 = 0.0
    for t in np.sort(np.unique(time[(event == 1) & (time <= horizon)])):
        d = ((time == t) & (event == 1)).sum()
        denom = expxb[time >= t].sum()
        if denom > 0:
            h0 += d / denom
    return h0


def recalibrated_predictions(time, event, score, horizon):
    score_z = zscore(score).to_numpy(float)
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(score_z)
    if mask.sum() < 30 or event[mask].sum() < 10:
        return None
    model = cox_ph(time[mask], event[mask], score_z[mask].reshape(-1, 1), ["risk_score_per_SD"])
    if not model or "risk_score_per_SD_beta" not in model:
        return None
    beta = model["risk_score_per_SD_beta"]
    h0 = baseline_cumhaz_at(time[mask], event[mask], score_z[mask], beta, horizon)
    surv = np.full(len(time), np.nan)
    surv[mask] = np.exp(-h0 * np.exp(np.clip(beta * score_z[mask], -50, 50)))
    return {
        "survival_probability": surv,
        "risk_probability": 1.0 - surv,
        "slope_beta": beta,
        "slope_HR": model.get("risk_score_per_SD_HR", np.nan),
        "slope_p": model.get("risk_score_per_SD_p", np.nan),
        "N_model": model.get("N_model", np.nan),
        "Events_model": model.get("Events_model", np.nan),
    }


def ipcw_brier(time, event, predicted_survival, horizon):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=float)
    pred = np.asarray(predicted_survival, dtype=float)
    mask = np.isfinite(time) & np.isfinite(event) & np.isfinite(pred)
    time, event, pred = time[mask], event[mask], pred[mask]
    if len(time) < 30:
        return np.nan, 0
    crows = censoring_km_values(time, event)
    gt = censoring_surv(crows, horizon)
    total = 0.0
    used = 0
    n = len(time)
    for ti, ei, pi in zip(time, event, pred):
        if ti <= horizon and ei == 1:
            g = censoring_surv(crows, ti, before=True)
            total += ((0.0 - pi) ** 2) / g
            used += 1
        elif ti > horizon:
            total += ((1.0 - pi) ** 2) / gt
            used += 1
    return total / n if used else np.nan, used


def build_metadata():
    articles = read_csv("article_extraction_qc.csv")
    readiness = read_csv("luad_signature_readiness.csv")
    selection = read_csv("luad_geo_manuscript_ready_signature_selection.csv")
    geo = read_csv("geo_external_validation_results.csv")
    pairwise = read_csv("luad_geo_manuscript_pairwise_evidence.csv")
    clinical = read_csv("luad_geo_main_signature_clinical_adjusted_cox.csv")
    robustness = read_csv("luad_geo_main_signature_robustness_summary.csv")
    integrated = read_csv("integrated_main_signature_manuscript_evidence.csv")
    coverage = read_csv("signature_geo_gene_coverage.csv")
    cohort = read_csv("manuscript_geo_cohort_source_table.csv")
    tcga = read_csv("tcga_luad_main_signature_gdc_star_tpm_sanity_results.csv")
    loco = read_csv("main_signature_leave_one_cohort_out.csv")
    scores = read_csv("geo_external_validation_sample_scores.csv")
    clin_samples = read_csv("geo_clinical_samples.csv")
    data = {
        "articles": articles,
        "readiness": readiness,
        "selection": selection,
        "geo": geo,
        "pairwise": pairwise,
        "clinical": clinical,
        "robustness": robustness,
        "integrated": integrated,
        "coverage": coverage,
        "cohort": cohort,
        "tcga": tcga,
        "loco": loco,
        "scores": scores,
        "clin_samples": clin_samples,
    }
    apply_manual_cutoff_overrides(data)
    return data


def apply_manual_cutoff_overrides(data):
    """Carry forward completed manual QC that intentionally was not back-written to extraction CSVs."""
    a097_mask_article = data["articles"]["Article_ID"].eq("A097")
    if a097_mask_article.any():
        data["articles"].loc[a097_mask_article, "Cutoff_available"] = "yes"
        data["articles"].loc[a097_mask_article, "Cutoff_value"] = 1.66
        data["articles"].loc[a097_mask_article, "Cutoff_type"] = "median / training-derived"
        data["articles"].loc[
            a097_mask_article, "Risk_group_definition"
        ] = "Manual PDF verification: TCGA/training risk-score cutoff 1.66; Low <1.66, High >=1.66. Fixed cutoff collapses processed GEO cohorts and is not treated as transportable."
    a097_mask_ready = data["readiness"]["Article_ID"].eq("A097")
    if a097_mask_ready.any():
        data["readiness"].loc[a097_mask_ready, "Cutoff_value"] = 1.66
        data["readiness"].loc[a097_mask_ready, "Cutoff_type"] = "median / training-derived"


def ladder_for_signatures(data):
    readiness = data["readiness"].copy()
    selection = data["selection"].copy()
    geo = data["geo"].copy()
    clinical = data["clinical"].copy()
    tcga = data["tcga"].copy()
    loco = data["loco"].copy()

    scoreable_status = {"ok", "single_group_after_cutoff"}
    scoreable = geo[geo["status"].isin(scoreable_status)].groupby("Signature_ID").size().rename("Externally_scoreable_tests")
    ok = geo[geo["status"].eq("ok")].groupby("Signature_ID").agg(
        GEO_tests_ok=("GSE", "nunique"),
        External_samples=("N", "sum"),
        External_events=("Events", "sum"),
        Continuous_sig_p05_n=("Cox_continuous_p", lambda s: safe_num(s).lt(0.05).sum()),
        Continuous_HR_gt1_rate=("Cox_continuous_HR_per_SD", lambda s: safe_num(s).gt(1).mean()),
        Median_C_index=("C_index", "median"),
        Group_sig_p05_n=("Cox_high_vs_low_p", lambda s: safe_num(s).lt(0.05).sum()),
    )
    numeric_ok = geo[geo["Cutoff_rule_used"].eq("reported_numeric_cutoff") & geo["status"].eq("ok")]
    numeric_ok_count = numeric_ok.groupby("Signature_ID").size().rename("Reported_numeric_cutoff_ok_tests")
    numeric_collapse = geo[geo["status"].eq("single_group_after_cutoff")].groupby("Signature_ID").size().rename("Single_group_collapse_tests")
    sel_cols = [
        "Signature_ID",
        "Final_analysis_tier",
        "Recommendation_Tier",
        "Recommended_manuscript_role",
        "Reported_numeric_cutoff_tests",
        "Single_group_after_cutoff_rows",
    ]
    out = readiness.merge(selection[[c for c in sel_cols if c in selection.columns]], on="Signature_ID", how="left")
    out = out.merge(scoreable, on="Signature_ID", how="left")
    out = out.merge(ok, on="Signature_ID", how="left")
    out = out.merge(numeric_ok_count, on="Signature_ID", how="left")
    out = out.merge(numeric_collapse, on="Signature_ID", how="left")
    out = out.merge(
        clinical.groupby("Signature_ID").agg(
            Clinical_adjusted_models=("GSE", "nunique"),
            Clinical_adjusted_p05_n=("risk_score_per_SD_p", lambda s: safe_num(s).lt(0.05).sum()),
            Delta_C_index_positive_n=("Delta_C_index_vs_clinical", lambda s: safe_num(s).gt(0).sum()),
        ),
        on="Signature_ID",
        how="left",
    )
    out = out.merge(
        tcga[tcga["status"].eq("ok")].groupby("Signature_ID").agg(
            TCGA_sanity_ok=("status", "size"),
            TCGA_C_index=("C_index", "max"),
            TCGA_continuous_p=("Continuous_p", "min"),
        ),
        on="Signature_ID",
        how="left",
    )
    out = out.merge(
        loco.groupby("Signature_ID").agg(LOCO_all_direction_retained=("Direction_retained_HR_gt1", lambda s: (s.astype(str).str.lower() == "yes").all())),
        on="Signature_ID",
        how="left",
    )
    numeric_cols = [
        "Externally_scoreable_tests",
        "GEO_tests_ok",
        "External_samples",
        "External_events",
        "Continuous_sig_p05_n",
        "Group_sig_p05_n",
        "Reported_numeric_cutoff_ok_tests",
        "Single_group_collapse_tests",
        "Clinical_adjusted_models",
        "Clinical_adjusted_p05_n",
        "Delta_C_index_positive_n",
        "TCGA_sanity_ok",
    ]
    for col in numeric_cols:
        if col in out:
            out[col] = out[col].fillna(0)

    def assign(row):
        level = 1
        label = "L1_formula_reproducible"
        if row.get("Externally_scoreable_tests", 0) > 0:
            level, label = 2, "L2_externally_scoreable"
        score_portable = (
            row.get("GEO_tests_ok", 0) >= 3
            and row.get("Continuous_sig_p05_n", 0) >= 2
            and row.get("Continuous_HR_gt1_rate", 0) >= 0.6
            and row.get("Median_C_index", 0) > 0.55
        )
        if score_portable:
            level, label = 3, "L3_score_portable"
        cutoff_transportable = (
            row.get("Reported_numeric_cutoff_ok_tests", 0) >= 3
            and row.get("Single_group_collapse_tests", 0) == 0
            and row.get("Group_sig_p05_n", 0) >= 2
        )
        if cutoff_transportable:
            level, label = 4, "L4_cutoff_transportable_candidate"
        deployability_candidate = (
            cutoff_transportable
            and row.get("Clinical_adjusted_p05_n", 0) >= 2
            and row.get("Delta_C_index_positive_n", 0) >= 2
            and row.get("TCGA_sanity_ok", 0) > 0
            and bool(row.get("LOCO_all_direction_retained", False))
        )
        if deployability_candidate:
            level, label = 5, "L5_clinical_deployability_candidate"
        return pd.Series({"Evidence_ladder_level": level, "Evidence_ladder_label": label})

    out = pd.concat([out, out.apply(assign, axis=1)], axis=1)
    return out


def ladder_for_articles(data, sig_ladder):
    articles = data["articles"].copy()
    max_by_article = sig_ladder.groupby("Article_ID").agg(
        Max_signature_ladder_level=("Evidence_ladder_level", "max"),
        Max_signature_ladder_label=("Evidence_ladder_label", "last"),
        Formula_reproducible_signature_count=("Signature_ID", "nunique"),
    )
    articles = articles.merge(max_by_article, on="Article_ID", how="left")
    articles["Formula_reproducible_signature_count"] = articles["Formula_reproducible_signature_count"].fillna(0).astype(int)
    articles["Evidence_ladder_level"] = np.where(
        articles["Formula_reproducible_signature_count"].gt(0),
        articles["Max_signature_ladder_level"],
        0,
    )
    labels = {
        0: "L0_non_reproducible_or_insufficiently_reported",
        1: "L1_formula_reproducible",
        2: "L2_externally_scoreable",
        3: "L3_score_portable",
        4: "L4_cutoff_transportable_candidate",
        5: "L5_clinical_deployability_candidate",
    }
    articles["Evidence_ladder_label"] = articles["Evidence_ladder_level"].fillna(0).astype(int).map(labels)
    return articles


def build_evidence_ladder(data):
    sig = ladder_for_signatures(data)
    art = ladder_for_articles(data, sig)
    sig.to_csv(OUT / "evidence_ladder_all_signatures.csv", index=False, encoding="utf-8")
    art.to_csv(OUT / "evidence_ladder_all_articles.csv", index=False, encoding="utf-8")

    levels = [
        "L0_non_reproducible_or_insufficiently_reported",
        "L1_formula_reproducible",
        "L2_externally_scoreable",
        "L3_score_portable",
        "L4_cutoff_transportable_candidate",
        "L5_clinical_deployability_candidate",
    ]
    rows = []
    for entity, df in [("article", art), ("signature", sig)]:
        counts = df["Evidence_ladder_label"].value_counts()
        for level in levels:
            rows.append({"Entity": entity, "Evidence_ladder_label": level, "Count": int(counts.get(level, 0))})
    pd.DataFrame(rows).to_csv(OUT / "evidence_ladder_summary_counts.csv", index=False, encoding="utf-8")

    transition_rows = []
    sig_total = len(sig)
    art_total = len(art)
    transitions = [
        ("Article inventory", "L0_non_reproducible_or_insufficiently_reported", int((art["Evidence_ladder_level"] == 0).sum())),
        ("Article inventory", "L1_formula_reproducible", int((art["Evidence_ladder_level"] >= 1).sum())),
        ("L1_formula_reproducible", "L2_externally_scoreable", int((sig["Evidence_ladder_level"] >= 2).sum())),
        ("L2_externally_scoreable", "L3_score_portable", int((sig["Evidence_ladder_level"] >= 3).sum())),
        ("L3_score_portable", "L4_cutoff_transportable_candidate", int((sig["Evidence_ladder_level"] >= 4).sum())),
        ("L4_cutoff_transportable_candidate", "L5_clinical_deployability_candidate", int((sig["Evidence_ladder_level"] >= 5).sum())),
    ]
    for source, target, n in transitions:
        transition_rows.append({"Source": source, "Target": target, "Count": n, "Signature_denominator": sig_total, "Article_denominator": art_total})
    trans = pd.DataFrame(transition_rows)
    trans.to_csv(OUT / "evidence_ladder_transition_matrix.csv", index=False, encoding="utf-8")
    trans.rename(columns={"Source": "source", "Target": "target", "Count": "value"}).to_csv(
        OUT / "evidence_ladder_sankey_data.csv", index=False, encoding="utf-8"
    )
    return art, sig


def cutoff_derivation(cutoff_type):
    text = clean_text(cutoff_type).lower()
    parts = []
    for key in ["median", "training-derived", "optimal", "surv_cutpoint", "x-tile", "youden"]:
        if key in text:
            parts.append(key)
    return "; ".join(parts) if parts else ("not reported" if not text else "unclear")


def build_cutoff_audit(data, sig_ladder):
    readiness = data["readiness"].copy()
    geo = data["geo"].copy()
    cohort = data["cohort"].copy()
    samples = data["scores"].copy()
    a097_stress_path = WORK / "a097_fixed_cutoff_1p66_stress_test.csv"
    a097_stress = pd.read_csv(a097_stress_path) if a097_stress_path.exists() else pd.DataFrame()

    transform_summary = geo.groupby("Signature_ID").agg(
        GEO_preprocessing_observed=("Scoring_transform", lambda s: "; ".join(sorted(set(clean_text(x) for x in s if has_value(x))))),
        Reported_cutoff_GEO_ok_tests=("Cutoff_rule_used", lambda s: (s == "reported_numeric_cutoff").sum()),
        Single_group_collapse_count=("status", lambda s: (s == "single_group_after_cutoff").sum()),
    )
    status_summary = geo.groupby("Signature_ID")["status"].value_counts().unstack(fill_value=0)
    out = readiness.merge(transform_summary, on="Signature_ID", how="left")
    out = out.merge(status_summary.add_prefix("GEO_status_"), on="Signature_ID", how="left")
    out["reported_numeric_cutoff"] = safe_num(out["Cutoff_value"]).notna()
    if not a097_stress.empty:
        a097_sig = "A097_PMID36971680"
        out.loc[out["Signature_ID"].eq(a097_sig), "Single_group_collapse_count"] = len(a097_stress)
        out.loc[out["Signature_ID"].eq(a097_sig), "A097_manual_cutoff_stress_test"] = "fixed cutoff 1.66 collapses all tested GEO cohorts to one group"
    else:
        out["A097_manual_cutoff_stress_test"] = ""
    out["cutoff_rule"] = out["Cutoff_type"]
    out["cutoff_derivation"] = out["Cutoff_type"].apply(cutoff_derivation)
    out["expression_scale_reported"] = out.apply(
        lambda r: "partly inferable from formula text" if has_value(r.get("Cutoff_value")) else "not sufficiently standardized for cross-study cutoff use",
        axis=1,
    )
    out["preprocessing_reported"] = np.where(out["GEO_preprocessing_observed"].fillna("").ne(""), out["GEO_preprocessing_observed"], "not harmonized across GEO")
    out["reported_cutoff_applicable_to_GEO"] = out["Reported_cutoff_GEO_ok_tests"].fillna(0).gt(0)
    out["single_group_collapse_observed"] = out["Single_group_collapse_count"].fillna(0).gt(0)
    out.loc[out["Signature_ID"].eq("A097_PMID36971680"), "reported_cutoff_applicable_to_GEO"] = False
    out["cutoff_transportability_interpretation"] = np.where(
        out["reported_cutoff_applicable_to_GEO"] & ~out["single_group_collapse_observed"],
        "reported numeric cutoff testable in at least one GEO cohort",
        np.where(out["reported_numeric_cutoff"], "reported numeric cutoff not robustly transportable in GEO", "no numeric cutoff for strict transportability testing"),
    )
    out.to_csv(OUT / "cutoff_reporting_audit_57.csv", index=False, encoding="utf-8")

    long_rows = []
    for _, row in geo[geo["Cutoff_rule_used"].eq("reported_numeric_cutoff") | geo["status"].eq("single_group_after_cutoff")].iterrows():
        n = row.get("N")
        high = row.get("High_n")
        low = row.get("Low_n")
        if pd.notna(n) and n:
            high_prop = high / n if pd.notna(high) else np.nan
            low_prop = low / n if pd.notna(low) else np.nan
        else:
            high_prop = np.nan
            low_prop = np.nan
        long_rows.append(
            {
                "Signature_ID": row["Signature_ID"],
                "Article_ID": row["Article_ID"],
                "GSE": row["GSE"],
                "Status": row["status"],
                "Cutoff_used": row.get("Cutoff_used"),
                "N": n,
                "High_n": high,
                "Low_n": low,
                "reported_cutoff_percentile_low_fraction": low_prop,
                "absolute_percentile_drift_from_median": abs(low_prop - 0.5) if pd.notna(low_prop) else np.nan,
            }
        )
    drift = pd.DataFrame(long_rows)
    if not a097_stress.empty:
        stress_rows = []
        for _, r in a097_stress.iterrows():
            n = r.get("N_scored")
            low = r.get("Low_n_score_lt_1p66")
            high = r.get("High_n_score_ge_1p66")
            low_prop = low / n if pd.notna(n) and n else np.nan
            stress_rows.append(
                {
                    "Signature_ID": r["Signature_ID"],
                    "Article_ID": "A097",
                    "GSE": r["GSE"],
                    "Status": r["Fixed_cutoff_stress_status"],
                    "Cutoff_used": r["Fixed_cutoff_from_uploaded_pdf"],
                    "N": n,
                    "High_n": high,
                    "Low_n": low,
                    "reported_cutoff_percentile_low_fraction": low_prop,
                    "absolute_percentile_drift_from_median": abs(low_prop - 0.5) if pd.notna(low_prop) else np.nan,
                    "Manual_override_note": "A097 PDF cutoff 1.66 confirmed; fixed cutoff stress test shows single-group collapse on processed GEO scales.",
                }
            )
        drift = pd.concat([drift, pd.DataFrame(stress_rows)], ignore_index=True)
    if "Manual_override_note" not in drift.columns:
        drift["Manual_override_note"] = ""
    drift = drift.merge(cohort[["GSE", "GPL"]], on="GSE", how="left")
    drift.to_csv(OUT / "cutoff_percentile_drift_long.csv", index=False, encoding="utf-8")
    high_prop = drift.copy()
    high_prop["reported_cutoff_high_risk_proportion"] = high_prop["High_n"] / high_prop["N"]
    high_prop.to_csv(OUT / "reported_cutoff_high_risk_proportion_long.csv", index=False, encoding="utf-8")
    collapse = geo.groupby(["Signature_ID", "Article_ID"]).agg(
        tested_GEOs=("GSE", "nunique"),
        reported_cutoff_ok_tests=("Cutoff_rule_used", lambda s: (s == "reported_numeric_cutoff").sum()),
        single_group_collapse_tests=("status", lambda s: (s == "single_group_after_cutoff").sum()),
        insufficient_gene_coverage_tests=("status", lambda s: (s == "insufficient_gene_coverage").sum()),
    ).reset_index()
    if not a097_stress.empty:
        mask = collapse["Signature_ID"].eq("A097_PMID36971680")
        if mask.any():
            collapse.loc[mask, "reported_cutoff_ok_tests"] = 0
            collapse.loc[mask, "single_group_collapse_tests"] = len(a097_stress)
        else:
            collapse = pd.concat(
                [
                    collapse,
                    pd.DataFrame(
                        [
                            {
                                "Signature_ID": "A097_PMID36971680",
                                "Article_ID": "A097",
                                "tested_GEOs": len(a097_stress),
                                "reported_cutoff_ok_tests": 0,
                                "single_group_collapse_tests": len(a097_stress),
                                "insufficient_gene_coverage_tests": 0,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
    collapse.to_csv(OUT / "single_group_collapse_summary.csv", index=False, encoding="utf-8")
    summary = out.groupby(["reported_numeric_cutoff", "reported_cutoff_applicable_to_GEO", "single_group_collapse_observed"]).size().reset_index(name="Signature_count")
    summary.to_csv(OUT / "cutoff_transportability_summary.csv", index=False, encoding="utf-8")
    return out


def reporting_quality_items(df):
    item_defs = {
        "gene_list": df["Gene_list"].apply(has_value) if "Gene_list" in df else df["Gene_count"].fillna(0).gt(0),
        "coefficients": df["Coefficients_available"].apply(yn) if "Coefficients_available" in df else df["Numeric_coefficient_count"].fillna(0).gt(0),
        "formula": df["Formula_available"].apply(yn) if "Formula_available" in df else df["Ready_for_expression_scoring"].apply(yn),
        "numeric_cutoff": safe_num(df["Cutoff_value"]).notna() if "Cutoff_value" in df else pd.Series(False, index=df.index),
        "cutoff_derivation": df["Cutoff_type"].apply(has_value) if "Cutoff_type" in df else pd.Series(False, index=df.index),
        "training_cohort": df["Training_cohort"].apply(has_value),
        "validation_cohort": df["Validation_cohort"].apply(has_value),
        "endpoint": df["Endpoint"].apply(has_value),
        "preprocessing": (
            df.get("Risk_score_formula", pd.Series("", index=df.index)).astype(str).str.lower().str.contains("log|scale|normalize|fpkm|tpm|count|z-score|zscore|expression")
            | df.get("Evidence_formula_context", pd.Series("", index=df.index)).astype(str).str.lower().str.contains("log|scale|normalize|fpkm|tpm|count|z-score|zscore|expression")
        ),
        "code_data_availability": (
            df.get("Public_data_used", pd.Series("", index=df.index)).apply(has_value)
            | df.get("Source_files", pd.Series("", index=df.index)).apply(has_value)
        ),
    }
    return pd.DataFrame({k: v.astype(bool).astype(int) for k, v in item_defs.items()}, index=df.index)


def build_reporting_quality(data, art_ladder, sig_ladder):
    articles = data["articles"].copy()
    readiness = data["readiness"].copy()
    art_items = reporting_quality_items(articles)
    art_out = pd.concat([articles[["Article_ID", "PMID", "Title", "Year", "Journal", "Endpoint", "Reproducible"]], art_items], axis=1)
    art_out["Reporting_quality_score_0_10"] = art_items.sum(axis=1)
    art_out = art_out.merge(art_ladder[["Article_ID", "Evidence_ladder_level", "Evidence_ladder_label"]], on="Article_ID", how="left")
    art_out.to_csv(OUT / "reporting_quality_score_all_articles.csv", index=False, encoding="utf-8")

    sig_base = readiness.merge(data["articles"].drop_duplicates("Article_ID"), on=["Article_ID", "PMID", "Title", "Endpoint", "Training_cohort", "Validation_cohort", "Cutoff_type", "Cutoff_value"], how="left", suffixes=("", "_article"))
    sig_items = reporting_quality_items(sig_base)
    sig_out = pd.concat([readiness[["Signature_ID", "Article_ID", "PMID", "Title", "Endpoint", "Gene_count", "Cutoff_type", "Cutoff_value"]], sig_items], axis=1)
    sig_out["Reporting_quality_score_0_10"] = sig_items.sum(axis=1)
    sig_out = sig_out.merge(sig_ladder[["Signature_ID", "Evidence_ladder_level", "Evidence_ladder_label"]], on="Signature_ID", how="left")
    sig_out.to_csv(OUT / "reporting_quality_score_by_signature.csv", index=False, encoding="utf-8")

    summary = pd.DataFrame(
        [
            {"Entity": "article", "N": len(art_out), "Median_score": art_out["Reporting_quality_score_0_10"].median(), "Mean_score": art_out["Reporting_quality_score_0_10"].mean()},
            {"Entity": "signature", "N": len(sig_out), "Median_score": sig_out["Reporting_quality_score_0_10"].median(), "Mean_score": sig_out["Reporting_quality_score_0_10"].mean()},
        ]
    )
    summary.to_csv(OUT / "reporting_quality_summary.csv", index=False, encoding="utf-8")
    art_out.groupby("Year").agg(
        Article_count=("Article_ID", "nunique"),
        Median_reporting_score=("Reporting_quality_score_0_10", "median"),
        Mean_reporting_score=("Reporting_quality_score_0_10", "mean"),
        Formula_reproducible_count=("Reproducible", lambda s: (s.astype(str).str.lower() == "yes").sum()),
    ).reset_index().to_csv(OUT / "reporting_score_by_year.csv", index=False, encoding="utf-8")
    pd.concat(
        [
            art_out.assign(Entity="article")[["Entity", "Article_ID", "Reporting_quality_score_0_10", "Evidence_ladder_level", "Evidence_ladder_label"]],
            sig_out.assign(Entity="signature")[["Entity", "Signature_ID", "Reporting_quality_score_0_10", "Evidence_ladder_level", "Evidence_ladder_label"]].rename(columns={"Signature_ID": "Article_ID"}),
        ],
        ignore_index=True,
    ).to_csv(OUT / "reporting_score_vs_evidence_ladder.csv", index=False, encoding="utf-8")
    comp = art_out.assign(Reproducible_binary=art_out["Reproducible"].astype(str).str.lower().eq("yes")).groupby("Reproducible_binary").agg(
        Article_count=("Article_ID", "nunique"),
        Median_reporting_score=("Reporting_quality_score_0_10", "median"),
        Mean_reporting_score=("Reporting_quality_score_0_10", "mean"),
        Numeric_cutoff_rate=("numeric_cutoff", "mean"),
        Formula_rate=("formula", "mean"),
        Coefficients_rate=("coefficients", "mean"),
        Preprocessing_rate=("preprocessing", "mean"),
    ).reset_index()
    comp.to_csv(OUT / "reproducible_vs_nonreproducible_reporting_comparison.csv", index=False, encoding="utf-8")
    return art_out, sig_out


def build_all57_landscape(data, sig_ladder):
    geo = data["geo"].copy()
    coverage = data["coverage"].copy()
    readiness = data["readiness"].copy()
    selection = data["selection"].copy()
    primary = set(selection[selection["Final_analysis_tier"].isin(["Tier_1_strict_reported_cutoff", "Tier_2_primary_median_fallback"])]["Signature_ID"])

    scoreability = geo.groupby(["Signature_ID", "Article_ID"]).agg(
        total_candidate_GEO_tests=("GSE", "nunique"),
        ok_tests=("status", lambda s: (s == "ok").sum()),
        scoreable_or_cutoff_collapse_tests=("status", lambda s: s.isin(["ok", "single_group_after_cutoff"]).sum()),
        insufficient_gene_coverage_tests=("status", lambda s: (s == "insufficient_gene_coverage").sum()),
        skipped_training_overlap_tests=("status", lambda s: (s == "skipped_training_overlap").sum()),
        single_group_after_cutoff_tests=("status", lambda s: (s == "single_group_after_cutoff").sum()),
        score_contains_missing_tests=("status", lambda s: (s == "score_contains_missing").sum()),
    ).reset_index()
    scoreability.to_csv(OUT / "all_57_signature_scoreability.csv", index=False, encoding="utf-8")

    gene_cov = coverage.groupby(["Signature_ID", "Article_ID"]).agg(
        signature_gene_count=("Signature_gene_count", "max"),
        median_gene_coverage=("Coverage_fraction", "median"),
        min_gene_coverage=("Coverage_fraction", "min"),
        full_coverage_GEOs=("Coverage_fraction", lambda s: safe_num(s).eq(1.0).sum()),
        any_missing_gene_GEOs=("Missing_gene_count", lambda s: safe_num(s).gt(0).sum()),
        recurrent_missing_genes=("Missing_genes", lambda s: "; ".join(sorted(set("; ".join(s.dropna().astype(str)).split("; ")) - {""}))[:1000]),
    ).reset_index()
    gene_cov.to_csv(OUT / "all_57_gene_coverage_summary.csv", index=False, encoding="utf-8")

    ok = geo[geo["status"].eq("ok")].copy()
    perf = ok.groupby(["Signature_ID", "Article_ID"]).agg(
        evaluable_cohort_count=("GSE", "nunique"),
        external_samples=("N", "sum"),
        external_events=("Events", "sum"),
        grouped_significant_rate=("Cox_high_vs_low_p", lambda s: safe_num(s).lt(0.05).mean()),
        continuous_significant_rate=("Cox_continuous_p", lambda s: safe_num(s).lt(0.05).mean()),
        median_C_index=("C_index", "median"),
        median_group_HR=("Cox_high_vs_low_HR", "median"),
    ).reset_index()
    perf = perf.merge(readiness[["Signature_ID", "Cutoff_type", "Cutoff_value"]], on="Signature_ID", how="left")
    perf = perf.merge(sig_ladder[["Signature_ID", "Evidence_ladder_level", "Evidence_ladder_label"]], on="Signature_ID", how="left")
    perf["reported_cutoff_status"] = np.where(safe_num(perf["Cutoff_value"]).notna(), "numeric cutoff reported", "numeric cutoff absent or median/rule only")
    perf.to_csv(OUT / "all_57_GEO_performance_summary.csv", index=False, encoding="utf-8")

    ok[["Signature_ID", "Article_ID", "GSE", "Cox_high_vs_low_HR", "Cox_high_vs_low_p", "Cutoff_rule_used", "N", "Events"]].to_csv(
        OUT / "all_pairwise_HR_distribution.csv", index=False, encoding="utf-8"
    )
    ok[["Signature_ID", "Article_ID", "GSE", "C_index", "Cox_continuous_HR_per_SD", "Cox_continuous_p", "N", "Events"]].to_csv(
        OUT / "all_pairwise_Cindex_distribution.csv", index=False, encoding="utf-8"
    )

    reasons = selection[~selection["Signature_ID"].isin(primary)].copy()
    def reason(row):
        if row.get("Tested_cohorts_ok", 0) == 0:
            return "not externally scoreable in the first-pass GEO panel"
        if row.get("Single_group_after_cutoff_rows", 0) > 0:
            return "reported cutoff produced single-group collapse in at least one GEO test"
        if row.get("Continuous_sig_p05_n", 0) < 3:
            return "weaker continuous score-portability signal than primary candidates"
        if row.get("Median_C_index", 0) < 0.58:
            return "lower discrimination than primary-candidate set"
        return "secondary or sensitivity-tier signal retained outside primary candidate set"
    reasons["Reason_not_advanced_to_primary"] = reasons.apply(reason, axis=1)
    reasons.to_csv(OUT / "reasons_not_advanced_to_primary.csv", index=False, encoding="utf-8")

    primary_just = selection[selection["Signature_ID"].isin(primary)].copy()
    primary_just["Primary_candidate_justification"] = primary_just.apply(
        lambda r: (
            "only strict reported-cutoff primary candidate"
            if r["Article_ID"] == "A131"
            else "selected for recurrent GEO score portability; grouped analyses use cohort-median fallback and do not test the original reported threshold"
        ),
        axis=1,
    )
    primary_just.to_csv(OUT / "primary_candidate_selection_justification.csv", index=False, encoding="utf-8")


def build_clinical_samples(data):
    clin = data["clin_samples"].copy()
    clin["age_z"] = zscore(clin["age"]).values
    clin["sex_male"] = clin["sex"].apply(sex_binary)
    clin["stage_major"] = clin["stage"].apply(stage_major)
    clin["stage_z"] = zscore(clin["stage_major"]).values
    return clin


def build_prediction_performance(data, primary_sigs):
    scores = data["scores"].copy()
    pairwise = data["pairwise"].copy()
    clin = build_clinical_samples(data)
    valid_pairs = pairwise[pairwise["Signature_ID"].isin(primary_sigs) & pairwise["status"].eq("ok") & pairwise["Is_main_analysis"].eq(True)][["Signature_ID", "Article_ID", "GSE"]].drop_duplicates()
    scores = scores.merge(valid_pairs, on=["Signature_ID", "GSE"], how="inner")
    scores = scores.merge(clin[["GSE", "geo_accession", "age_z", "sex_male", "stage_z"]], on=["GSE", "geo_accession"], how="left")

    auc_rows = []
    brier_rows = []
    cal_rows = []
    curve_rows = []
    decision_rows = []
    inc_rows = []
    lrt_rows = []

    for (sig, article, gse), sub in scores.groupby(["Signature_ID", "Article_ID", "GSE"]):
        sub = sub.copy()
        time = safe_num(sub["OS_time_months"]).to_numpy(float)
        event = safe_num(sub["OS_event"]).to_numpy(float)
        score = safe_num(sub["Risk_score"]).to_numpy(float)
        score_z = zscore(score).to_numpy(float)
        n = int(np.isfinite(time + event + score).sum())
        events = int(np.nansum(event[np.isfinite(event)]))
        for horizon in TIMEPOINTS_AUC:
            auc, cases, controls = dynamic_auc(time, event, score, horizon)
            auc_rows.append(
                {
                    "Article_ID": article,
                    "Signature_ID": sig,
                    "GSE": gse,
                    "time_months": horizon,
                    "time_years": horizon / 12,
                    "dynamic_AUC_complete_case": auc,
                    "case_count_event_by_time": cases,
                    "control_count_survived_beyond_time": controls,
                    "N_scored": n,
                    "Events": events,
                    "interpretation": "cumulative/dynamic AUC using externally scored risk score",
                }
            )
        for horizon in TIMEPOINTS_PRED:
            pred = recalibrated_predictions(time, event, score, horizon)
            if pred is None:
                brier = np.nan
                used = 0
                slope_beta = slope_hr = slope_p = np.nan
                risk_prob = np.full(len(sub), np.nan)
                surv_prob = np.full(len(sub), np.nan)
            else:
                surv_prob = pred["survival_probability"]
                risk_prob = pred["risk_probability"]
                brier, used = ipcw_brier(time, event, surv_prob, horizon)
                slope_beta = pred["slope_beta"]
                slope_hr = pred["slope_HR"]
                slope_p = pred["slope_p"]
            brier_rows.append(
                {
                    "Article_ID": article,
                    "Signature_ID": sig,
                    "GSE": gse,
                    "time_months": horizon,
                    "time_years": horizon / 12,
                    "IPCW_Brier_score_external_recalibrated": brier,
                    "IPCW_complete_cases_used": used,
                    "N_scored": n,
                    "Events": events,
                    "calibration_scope": "external recalibration assessment; original absolute-risk calibration not recoverable without published baseline survival",
                }
            )
            cal_rows.append(
                {
                    "Article_ID": article,
                    "Signature_ID": sig,
                    "GSE": gse,
                    "time_months": horizon,
                    "time_years": horizon / 12,
                    "external_recalibration_slope_beta": slope_beta,
                    "external_recalibration_HR_per_SD": slope_hr,
                    "external_recalibration_p": slope_p,
                    "N_model": n,
                    "Events": events,
                    "calibration_scope": "external recalibration assessment, not original model absolute-risk calibration",
                }
            )
            curve_df = pd.DataFrame({"time": time, "event": event, "risk": risk_prob, "surv": surv_prob})
            curve_df = curve_df[np.isfinite(curve_df["time"]) & np.isfinite(curve_df["event"]) & np.isfinite(curve_df["risk"])]
            if len(curve_df) >= 30 and curve_df["risk"].nunique() >= 4:
                try:
                    curve_df["risk_quantile"] = pd.qcut(curve_df["risk"], q=min(5, curve_df["risk"].nunique()), duplicates="drop")
                    for idx, (_, g) in enumerate(curve_df.groupby("risk_quantile", observed=True), start=1):
                        obs_surv = km_survival_at(g["time"], g["event"], horizon)
                        curve_rows.append(
                            {
                                "Article_ID": article,
                                "Signature_ID": sig,
                                "GSE": gse,
                                "time_months": horizon,
                                "risk_bin": idx,
                                "N_bin": len(g),
                                "Mean_predicted_risk": g["risk"].mean(),
                                "Mean_predicted_survival": g["surv"].mean(),
                                "Observed_KM_risk": 1 - obs_surv if pd.notna(obs_surv) else np.nan,
                                "Observed_KM_survival": obs_surv,
                            }
                        )
                except ValueError:
                    pass
            dca_df = pd.DataFrame({"time": time, "event": event, "risk": risk_prob})
            usable = dca_df[(np.isfinite(dca_df["risk"])) & (((dca_df["time"] <= horizon) & (dca_df["event"] == 1)) | (dca_df["time"] > horizon))].copy()
            if len(usable) >= 30:
                y = ((usable["time"] <= horizon) & (usable["event"] == 1)).astype(int)
                prevalence = y.mean()
                for threshold in np.arange(0.05, 0.51, 0.05):
                    pred_pos = usable["risk"] >= threshold
                    tp = ((pred_pos) & (y == 1)).sum()
                    fp = ((pred_pos) & (y == 0)).sum()
                    nb = tp / len(usable) - fp / len(usable) * threshold / (1 - threshold)
                    treat_all = prevalence - (1 - prevalence) * threshold / (1 - threshold)
                    decision_rows.append(
                        {
                            "Article_ID": article,
                            "Signature_ID": sig,
                            "GSE": gse,
                            "time_months": horizon,
                            "threshold_probability": round(float(threshold), 2),
                            "net_benefit_model_complete_case": nb,
                            "net_benefit_treat_all": treat_all,
                            "net_benefit_treat_none": 0.0,
                            "N_complete_case": len(usable),
                            "event_by_time_prevalence": prevalence,
                            "decision_curve_scope": "exploratory complete-case decision curve using external recalibrated risk",
                        }
                    )
        covariates = []
        labels = []
        for col, label in [("age_z", "age_per_SD"), ("sex_male", "sex_male"), ("stage_z", "stage_per_SD")]:
            if sub[col].notna().sum() >= 50 and sub[col].nunique(dropna=True) >= 2:
                covariates.append(col)
                labels.append(label)
        if covariates:
            comp = sub[["OS_time_months", "OS_event", "Risk_score"] + covariates].copy()
            comp["risk_z"] = zscore(comp["Risk_score"]).values
            comp = comp.dropna()
            if len(comp) >= 30 and comp["OS_event"].sum() >= 10:
                clin_model = cox_ph(comp["OS_time_months"], comp["OS_event"], comp[covariates].to_numpy(float), labels)
                full_model = cox_ph(comp["OS_time_months"], comp["OS_event"], comp[["risk_z"] + covariates].to_numpy(float), ["risk_score_per_SD"] + labels)
                if clin_model and full_model:
                    clin_c = c_index(comp["OS_time_months"], comp["OS_event"], clin_model["linear_predictor"])
                    full_c = c_index(comp["OS_time_months"], comp["OS_event"], full_model["linear_predictor"])
                    inc_rows.append(
                        {
                            "Article_ID": article,
                            "Signature_ID": sig,
                            "GSE": gse,
                            "Clinical_covariates": "; ".join(labels),
                            "N_model": len(comp),
                            "Events_model": int(comp["OS_event"].sum()),
                            "Clinical_only_C_index": clin_c,
                            "Clinical_plus_risk_C_index": full_c,
                            "Delta_C_index": full_c - clin_c,
                        }
                    )
                    lr = 2 * (full_model["loglik"] - clin_model["loglik"])
                    lrt_rows.append(
                        {
                            "Article_ID": article,
                            "Signature_ID": sig,
                            "GSE": gse,
                            "Clinical_covariates": "; ".join(labels),
                            "N_model": len(comp),
                            "Events_model": int(comp["OS_event"].sum()),
                            "Clinical_only_loglik": clin_model["loglik"],
                            "Clinical_plus_risk_loglik": full_model["loglik"],
                            "LR_statistic_df1": lr,
                            "LR_test_p": p_chisq1(lr),
                        }
                    )

    pd.DataFrame(auc_rows).to_csv(OUT / "time_dependent_auc_primary_candidates.csv", index=False, encoding="utf-8")
    pd.DataFrame(brier_rows).to_csv(OUT / "brier_score_primary_candidates.csv", index=False, encoding="utf-8")
    pd.DataFrame(cal_rows).to_csv(OUT / "calibration_metrics_primary_candidates.csv", index=False, encoding="utf-8")
    pd.DataFrame(curve_rows).to_csv(OUT / "calibration_curve_data_3yr_5yr.csv", index=False, encoding="utf-8")
    pd.DataFrame(inc_rows).to_csv(OUT / "incremental_cindex_by_signature_cohort.csv", index=False, encoding="utf-8")
    pd.DataFrame(lrt_rows).to_csv(OUT / "likelihood_ratio_test_clinical_vs_clinical_plus_risk.csv", index=False, encoding="utf-8")
    pd.DataFrame(decision_rows).to_csv(OUT / "decision_curve_data_3yr_5yr.csv", index=False, encoding="utf-8")


def build_deployability_checklist(data, sig_ladder):
    selection = data["selection"].copy()
    readiness = data["readiness"].copy()
    robustness = data["robustness"].copy()
    tcga = data["tcga"].copy()
    loco = data["loco"].copy()
    clinical = data["clinical"].copy()
    primary = selection[selection["Final_analysis_tier"].isin(["Tier_1_strict_reported_cutoff", "Tier_2_primary_median_fallback"])].copy()
    primary = primary[primary["Article_ID"].isin(PRIMARY_ARTICLES)].copy()
    primary = primary.merge(readiness[["Signature_ID", "Ready_for_expression_scoring", "Gene_count", "Unmatched_gene_count", "Endpoint", "Cutoff_value", "Cutoff_type"]], on="Signature_ID", how="left", suffixes=("", "_readiness"))
    primary = primary.merge(robustness, on=["Article_ID", "Signature_ID", "Final_analysis_tier"], how="left")
    primary = primary.merge(sig_ladder[["Signature_ID", "Evidence_ladder_level", "Evidence_ladder_label"]], on="Signature_ID", how="left")
    tcga_ok = set(tcga[tcga["status"].eq("ok")]["Signature_ID"])
    loco_ok = loco.groupby("Signature_ID")["Direction_retained_HR_gt1"].apply(lambda s: (s.astype(str).str.lower() == "yes").all()).to_dict()
    auc = pd.read_csv(OUT / "time_dependent_auc_primary_candidates.csv")
    brier = pd.read_csv(OUT / "brier_score_primary_candidates.csv")
    cal = pd.read_csv(OUT / "calibration_metrics_primary_candidates.csv")
    inc = pd.read_csv(OUT / "incremental_cindex_by_signature_cohort.csv")
    lrt = pd.read_csv(OUT / "likelihood_ratio_test_clinical_vs_clinical_plus_risk.csv")
    a097_stress_path = WORK / "a097_fixed_cutoff_1p66_stress_test.csv"
    cutoff_collapse_sigs = {"A097_PMID36971680"} if a097_stress_path.exists() else set()

    rows = []
    for _, r in primary.iterrows():
        sig = r["Signature_ID"]
        article = r["Article_ID"]
        strict = r["Final_analysis_tier"] == "Tier_1_strict_reported_cutoff"
        cutoff_value = r.get("Cutoff_value_readiness", r.get("Cutoff_value"))
        endpoint_value = r.get("Endpoint_readiness", r.get("Endpoint"))
        no_collapse = float(r.get("Single_group_after_cutoff_rows", 0) or 0) == 0 and sig not in cutoff_collapse_sigs
        cont_perf = float(r.get("Continuous_sig_p05_n", 0) or 0) >= 3
        adjusted = float(r.get("Adjusted_risk_p05_n", 0) or 0) >= 2
        delta = float(r.get("Delta_C_index_gt0_n", 0) or 0) >= 1
        time_auc_ok = sig in set(auc.dropna(subset=["dynamic_AUC_complete_case"])["Signature_ID"])
        brier_ok = sig in set(brier.dropna(subset=["IPCW_Brier_score_external_recalibrated"])["Signature_ID"])
        cal_ok = sig in set(cal.dropna(subset=["external_recalibration_slope_beta"])["Signature_ID"])
        checklist = {
            "Article_ID": article,
            "Signature_ID": sig,
            "complete_formula": "yes" if yn(r.get("Ready_for_expression_scoring")) else "no",
            "gene_availability": "yes" if float(r.get("Unmatched_gene_count", 1) or 0) == 0 else "partial",
            "endpoint_compatibility_with_GEO_OS": "yes" if "OS" in clean_text(endpoint_value) else "partial",
            "reported_numeric_cutoff": "yes" if pd.notna(cutoff_value) else "no",
            "strict_cutoff_testability": "yes" if strict else "no",
            "no_single_group_collapse": "yes" if no_collapse else "no",
            "continuous_GEO_performance": "yes" if cont_perf else "partial",
            "clinical_adjusted_signal": "yes" if adjusted else "partial",
            "delta_C_index_positive": "yes" if delta else "partial",
            "time_AUC_available": "yes" if time_auc_ok else "no",
            "Brier_available": "yes" if brier_ok else "no",
            "external_recalibration_available": "yes" if cal_ok else "no",
            "TCGA_formula_direction_sanity": "yes" if sig in tcga_ok else "no",
            "leave_one_cohort_out_direction_retained": "yes" if loco_ok.get(sig, False) else "no",
            "original_code_or_calculator": "manual_check_required",
            "assay_standardization": "no",
            "prospective_validation": "no",
            "clinical_threshold_justification": "partial" if strict else "no",
            "Evidence_ladder_label": r.get("Evidence_ladder_label"),
        }
        score = 0
        for key, val in checklist.items():
            if key in {"Article_ID", "Signature_ID", "Evidence_ladder_label"}:
                continue
            score += 1 if val == "yes" else 0.5 if val == "partial" else 0
        checklist["Deployability_score_0_18"] = score
        checklist["Deployability_interpretation"] = (
            "best candidate for further clinical deployability evaluation, but not ready for deployment"
            if article == "A131"
            else "score-portability candidate; lacks original cutoff transportability for clinical thresholding"
        )
        rows.append(checklist)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "clinical_deployability_checklist_primary10.csv", index=False, encoding="utf-8")
    score_cols = [c for c in out.columns if c not in {"Article_ID", "Signature_ID", "Evidence_ladder_label", "Deployability_score_0_18", "Deployability_interpretation"}]
    summary = out[["Article_ID", "Signature_ID", "Deployability_score_0_18", "Deployability_interpretation", "Evidence_ladder_label"]].sort_values("Deployability_score_0_18", ascending=False)
    summary.to_csv(OUT / "clinical_deployability_score_summary.csv", index=False, encoding="utf-8")
    heat = out.melt(id_vars=["Article_ID", "Signature_ID"], value_vars=score_cols, var_name="Deployability_domain", value_name="Status")
    heat["Heatmap_value"] = heat["Status"].map({"yes": 1.0, "partial": 0.5, "manual_check_required": 0.25, "no": 0.0}).fillna(0)
    heat.to_csv(OUT / "clinical_deployability_heatmap_data.csv", index=False, encoding="utf-8")

    lines = [
        "# Clinical deployability gap summary",
        "",
        "This project should be framed as a precision-oncology deployability audit of published LUAD transcriptomic prognostic signatures, not as another LUAD prognostic model.",
        "",
        "Key result: prognostic signal, score portability, cutoff transportability, and clinical deployability separate sharply. A131 is the only strict reported-cutoff primary candidate, while A129, A067, A097, A168, A001, and A059 are strong score-portability candidates whose grouped GEO results rely on cohort-median fallback.",
        "",
        "No primary candidate should be described as ready for clinical deployment. Common gaps include absent assay standardization, absent prospective validation, limited original-code/calculator availability, and inadequate justification of fixed clinical thresholds across expression platforms.",
        "",
        "Calibration outputs in this package are external recalibration assessments. They must not be described as original absolute-risk calibration because most source publications did not report the baseline survival function required to reconstruct absolute risk.",
    ]
    (OUT / "clinical_deployability_gap_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def build_dataset_reuse(data):
    readiness = data["readiness"].copy()
    plan = read_csv("luad_external_validation_plan.csv")
    edge_rows = []
    for _, r in readiness.iterrows():
        for role, col in [("training", "Training_cohort"), ("validation", "Validation_cohort")]:
            for ds in split_tokens(r.get(col)):
                edge_rows.append({"Article_ID": r["Article_ID"], "Signature_ID": r["Signature_ID"], "Dataset": normalize_dataset(ds), "Role": role})
    edges = pd.DataFrame(edge_rows)
    edges.to_csv(OUT / "dataset_signature_bipartite_network_edges.csv", index=False, encoding="utf-8")
    by_article = edges.groupby(["Article_ID", "Dataset", "Role"]).agg(Signature_count=("Signature_ID", "nunique")).reset_index()
    by_article.to_csv(OUT / "dataset_reuse_by_article.csv", index=False, encoding="utf-8")
    freq = edges.groupby(["Dataset", "Role"]).agg(Signature_count=("Signature_ID", "nunique"), Article_count=("Article_ID", "nunique")).reset_index()
    freq.to_csv(OUT / "dataset_reuse_frequency.csv", index=False, encoding="utf-8")

    sigs = sorted(readiness["Signature_ID"].unique())
    datasets = sorted(edges["Dataset"].dropna().unique())
    matrix_rows = []
    for sig in sigs:
        sig_edges = edges[edges["Signature_ID"].eq(sig)]
        row = {"Signature_ID": sig, "Article_ID": readiness.loc[readiness["Signature_ID"].eq(sig), "Article_ID"].iloc[0]}
        for ds in datasets:
            roles = sorted(sig_edges.loc[sig_edges["Dataset"].eq(ds), "Role"].unique())
            row[ds] = ";".join(roles)
        matrix_rows.append(row)
    pd.DataFrame(matrix_rows).to_csv(OUT / "signature_training_validation_overlap_matrix.csv", index=False, encoding="utf-8")

    excluded = []
    for _, r in plan.iterrows():
        for ds in split_tokens(r.get("Excluded_from_external_validation")):
            excluded.append({"Article_ID": r["Article_ID"], "Signature_ID": r["Signature_ID"], "Excluded_dataset": normalize_dataset(ds), "Reason": "training or internal-cohort overlap"})
    pd.DataFrame(excluded).to_csv(OUT / "excluded_external_tests_due_to_training_overlap.csv", index=False, encoding="utf-8")

    top = freq.sort_values(["Signature_count", "Article_count"], ascending=False).head(15)
    lines = [
        "# Dataset reuse summary",
        "",
        "Public LUAD prognostic-signature studies repeatedly reuse a small set of transcriptomic cohorts. This matters for precision-oncology deployability because apparent external validation can be weakened by training/validation overlap and endpoint heterogeneity.",
        "",
        "Most frequently reused datasets in the formula-reproducible set:",
        "",
    ]
    for _, r in top.iterrows():
        lines.append(f"- {r['Dataset']} ({r['Role']}): {int(r['Signature_count'])} signatures across {int(r['Article_count'])} articles")
    lines += [
        "",
        "Current GEO testing excluded signature-cohort combinations when the GEO cohort was identified as a training or internal-overlap dataset for that signature.",
    ]
    (OUT / "dataset_reuse_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def platform_class(gpl):
    gpl = clean_text(gpl).upper()
    if gpl in {"GPL570", "GPL96"}:
        return "Affymetrix microarray"
    if gpl in {"GPL6480", "GPL6884", "GPL15048"}:
        return "Illumina/other microarray"
    if "GDC" in gpl or "STAR" in gpl:
        return "RNA-seq/GDC STAR TPM sanity"
    return "Other/unknown microarray"


def build_platform_analysis(data):
    geo = data["geo"].copy()
    cohort = data["cohort"][["GSE", "GPL"]].drop_duplicates()
    primary = set(data["selection"][data["selection"]["Final_analysis_tier"].isin(["Tier_1_strict_reported_cutoff", "Tier_2_primary_median_fallback"])]["Signature_ID"])
    geo = geo.merge(cohort, on="GSE", how="left")
    geo["Platform_class"] = geo["GPL"].apply(platform_class)
    ok = geo[geo["status"].eq("ok")].copy()
    perf = ok.groupby(["Signature_ID", "Article_ID", "GPL", "Platform_class"]).agg(
        tests=("GSE", "nunique"),
        samples=("N", "sum"),
        events=("Events", "sum"),
        median_C_index=("C_index", "median"),
        median_group_HR=("Cox_high_vs_low_HR", "median"),
        grouped_p05_rate=("Cox_high_vs_low_p", lambda s: safe_num(s).lt(0.05).mean()),
        continuous_p05_rate=("Cox_continuous_p", lambda s: safe_num(s).lt(0.05).mean()),
    ).reset_index()
    perf.to_csv(OUT / "platform_stratified_performance.csv", index=False, encoding="utf-8")
    heter = ok.groupby(["Signature_ID", "Article_ID"]).agg(
        platform_count=("GPL", "nunique"),
        C_index_min=("C_index", "min"),
        C_index_max=("C_index", "max"),
        C_index_range=("C_index", lambda s: safe_num(s).max() - safe_num(s).min()),
        HR_min=("Cox_high_vs_low_HR", "min"),
        HR_max=("Cox_high_vs_low_HR", "max"),
    ).reset_index()
    heter.to_csv(OUT / "platform_heterogeneity_summary.csv", index=False, encoding="utf-8")
    collapse = geo.groupby(["Signature_ID", "Article_ID", "GPL", "Platform_class"]).agg(
        reported_cutoff_tests=("Cutoff_rule_used", lambda s: (s == "reported_numeric_cutoff").sum()),
        single_group_collapse_tests=("status", lambda s: (s == "single_group_after_cutoff").sum()),
    ).reset_index()
    collapse.to_csv(OUT / "platform_cutoff_collapse_summary.csv", index=False, encoding="utf-8")
    heat = perf[perf["Signature_ID"].isin(primary)].copy()
    heat.to_csv(OUT / "platform_specific_heatmap_data.csv", index=False, encoding="utf-8")


def endpoint_class(value):
    text = clean_text(value).upper()
    labels = [x for x in ["OS", "DFS", "RFS", "PFS", "DSS"] if re.search(rf"\b{x}\b", text)]
    if not labels:
        return "unclear"
    if labels == ["OS"]:
        return "OS"
    return "mixed_" + "_".join(labels)


def build_endpoint_audit(data, sig_ladder):
    articles = data["articles"].copy()
    readiness = data["readiness"].copy()
    geo = data["geo"].copy()
    articles["Endpoint_class"] = articles["Endpoint"].apply(endpoint_class)
    endpoint_dist = articles.groupby("Endpoint_class").agg(Article_count=("Article_ID", "nunique")).reset_index()
    endpoint_dist.to_csv(OUT / "endpoint_distribution_all_articles.csv", index=False, encoding="utf-8")
    sig_endpoint = readiness[["Signature_ID", "Article_ID", "PMID", "Title", "Endpoint"]].copy()
    sig_endpoint["Endpoint_class"] = sig_endpoint["Endpoint"].apply(endpoint_class)
    sig_endpoint.to_csv(OUT / "endpoint_by_signature.csv", index=False, encoding="utf-8")
    mismatch = sig_endpoint.copy()
    mismatch["GEO_validation_endpoint"] = "OS"
    mismatch["Endpoint_mismatch_with_GEO_OS"] = ~mismatch["Endpoint_class"].isin(["OS"])
    mismatch.to_csv(OUT / "endpoint_mismatch_with_GEO_OS_validation.csv", index=False, encoding="utf-8")
    os_sigs = set(sig_endpoint.loc[sig_endpoint["Endpoint_class"].eq("OS"), "Signature_ID"])
    ok = geo[geo["status"].eq("ok") & geo["Signature_ID"].isin(os_sigs)].copy()
    os_results = ok.groupby(["Signature_ID", "Article_ID"]).agg(
        OS_only_evaluable_GEOs=("GSE", "nunique"),
        OS_only_events=("Events", "sum"),
        OS_only_continuous_p05_n=("Cox_continuous_p", lambda s: safe_num(s).lt(0.05).sum()),
        OS_only_group_p05_n=("Cox_high_vs_low_p", lambda s: safe_num(s).lt(0.05).sum()),
        OS_only_median_C_index=("C_index", "median"),
    ).reset_index()
    os_results = os_results.merge(sig_ladder[["Signature_ID", "Evidence_ladder_label"]], on="Signature_ID", how="left")
    os_results.to_csv(OUT / "OS_only_signature_subset_results.csv", index=False, encoding="utf-8")
    lines = [
        "# Endpoint consistency summary",
        "",
        "GEO validation in this project uses overall survival. Source signatures that reported DFS, RFS, PFS, DSS, or mixed endpoints should therefore be interpreted as endpoint-mismatched when transported to GEO OS.",
        "",
        f"OS-only articles: {int((articles['Endpoint_class'] == 'OS').sum())} of {len(articles)}.",
        f"Endpoint-mismatched formula-reproducible signatures: {int(mismatch['Endpoint_mismatch_with_GEO_OS'].sum())} of {len(mismatch)}.",
        "",
        "Endpoint mismatch is a deployability issue: a score may retain prognostic biology while not being directly qualified for a different clinical decision endpoint.",
    ]
    (OUT / "endpoint_consistency_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_figure_data():
    copies = {
        "Figure1_evidence_ladder_workflow.csv": "evidence_ladder_sankey_data.csv",
        "Figure2_reporting_reproducibility_landscape.csv": "reporting_score_vs_evidence_ladder.csv",
        "Figure3_all57_external_score_portability_landscape.csv": "all_57_GEO_performance_summary.csv",
        "Figure4_cutoff_transportability_failure_anatomy.csv": "cutoff_percentile_drift_long.csv",
        "Figure5_primary_candidate_prediction_performance_auc.csv": "time_dependent_auc_primary_candidates.csv",
        "Figure5_primary_candidate_prediction_performance_brier.csv": "brier_score_primary_candidates.csv",
        "Figure6_clinical_deployability_gap_heatmap.csv": "clinical_deployability_heatmap_data.csv",
    }
    for dest, src in copies.items():
        src_path = OUT / src
        if src_path.exists():
            pd.read_csv(src_path).to_csv(FIG / dest, index=False, encoding="utf-8")


def build_markdown_outputs(data):
    counts = pd.read_csv(OUT / "evidence_ladder_summary_counts.csv")
    deploy = pd.read_csv(OUT / "clinical_deployability_score_summary.csv")
    cutoff = pd.read_csv(OUT / "cutoff_transportability_summary.csv")
    perf = pd.read_csv(OUT / "all_57_GEO_performance_summary.csv")
    top_primary = deploy.sort_values("Deployability_score_0_18", ascending=False).head(10)
    lines = [
        "# High-impact analysis summary",
        "",
        "This is not another LUAD prognostic model. It is a precision-oncology deployability audit of published transcriptomic prognostic signatures.",
        "",
        "## 1. Published LUAD prognostic signatures frequently lack complete reproducibility elements",
        "",
        "The evidence ladder separates article-level reporting from signature-level testability. Many articles remain at L0 because gene lists, coefficients, executable formulas, numeric cutoffs, or preprocessing details are insufficient for independent scoring.",
        "",
        "## 2. External scoreability is limited by gene coverage, cohort reuse, and endpoint inconsistency",
        "",
        "Among formula-reproducible signatures, GEO scoreability is constrained by missing genes on legacy microarray platforms, training/validation cohort reuse, and the use of GEO OS as a common endpoint for signatures originally developed under mixed endpoints.",
        "",
        "## 3. Score portability is preserved in a minority of signatures",
        "",
        "A129, A067, A097, A168, A001, and A059 show strong score-portability evidence, with A098 showing moderate support and A065/A173 lower-priority support. These are score-level signals, not proof that original clinical cutoffs transport.",
        "",
        "## 4. Reported cutoff transportability is rare and represents the major clinical translation bottleneck",
        "",
        "A131 is the only strict reported-cutoff primary candidate. Median-fallback grouped analyses are cohort-specific score-portability analyses and must not be written as original reported-cutoff validation.",
        "",
        "## 5. Clinical adjustment, expression-scale sensitivity, and TCGA sanity checks support selected signals but do not establish deployability",
        "",
        "Clinical-adjusted Cox models, delta C-index, z-score sensitivity, leave-one-cohort-out checks, and TCGA GDC STAR TPM formula-direction sanity checks support selected signals. TCGA-LUAD remains a formula-direction sanity check rather than a separate external-validation source.",
        "",
        "## 6. Evidence ladder and deployability checklist reveal a gap between prognostic biology and precision-oncology readiness",
        "",
        "No signature should be described as ready for clinical deployment. The audit shows that prognostic signal, cutoff transportability, and deployability are separable evidence levels.",
    ]
    (OUT / "high_impact_analysis_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    outline = [
        "# npj Precision Oncology manuscript outline",
        "",
        "## Title concept",
        "A precision-oncology deployability audit of published LUAD transcriptomic prognostic signatures",
        "",
        "## Central question",
        "Do published LUAD transcriptomic prognostic signatures support individualized prognostic stratification in precision oncology, and where do formula reproducibility, score portability, cutoff transportability, and clinical deployability diverge?",
        "",
        "## Results structure",
        "1. Published LUAD prognostic signatures frequently lack complete reproducibility elements.",
        "2. External scoreability is limited by gene coverage, cohort reuse, and endpoint inconsistency.",
        "3. Score portability is preserved in a minority of signatures.",
        "4. Reported cutoff transportability is rare and represents the major clinical translation bottleneck.",
        "5. Clinical adjustment, expression-scale sensitivity, and TCGA sanity checks support selected signals but do not establish deployability.",
        "6. Evidence ladder and deployability checklist reveal a gap between prognostic biology and precision-oncology readiness.",
        "",
        "## Required terminology guardrails",
        "- Use reported-cutoff validation only when the original reported numeric cutoff is applied.",
        "- Use score-portability analysis for median-fallback grouped GEO analyses.",
        "- Describe TCGA-LUAD as a formula-direction sanity check, not as a separate external-validation cohort.",
        "- Describe pooled HR as descriptive pooled HR unless a full formal meta-analysis is implemented.",
        "- Avoid valid/invalid binary framing; state that prognostic signal, cutoff transportability, and clinical deployability are separable evidence levels.",
    ]
    (OUT / "npj_precision_oncology_manuscript_outline.md").write_text("\n".join(outline) + "\n", encoding="utf-8")

    ladder_sig = counts[counts["Entity"].eq("signature")][["Evidence_ladder_label", "Count"]]
    report = [
        "# npj Precision Oncology readiness report",
        "",
        "## 1. Why the study fits npj Precision Oncology",
        "",
        "The upgraded study is framed around clinical deployability of published LUAD transcriptomic prognostic signatures, not model discovery. It audits whether signatures can move from reported formulas to external scoring, cutoff transportability, clinical-adjusted prognostic value, recalibrated prediction performance, and decision-relevant readiness.",
        "",
        "## 2. Risks for this journal",
        "",
        "- Most analyses remain retrospective and public-data based.",
        "- Original absolute-risk calibration is generally unrecoverable because baseline survival functions were not reported.",
        "- Assay standardization, locked clinical thresholds, and prospective validation are absent for most signatures.",
        "- Some source signatures reused TCGA-LUAD or GEO cohorts, so TCGA sanity checks and reused GEO cohorts require cautious wording.",
        "",
        "## 3. How the new analyses support precision-oncology deployability",
        "",
        "- Evidence ladder outputs quantify attrition from 178 full-text articles to 57 formula-reproducible signatures, 27 manuscript-triage signatures, and 10 primary candidates.",
        "- Cutoff audit isolates the key translation bottleneck: score portability often survives, while fixed published cutoffs rarely transport across GEO platforms.",
        "- Reporting quality scoring links missing reproducibility elements to failure points in external scoring and cutoff testing.",
        "- Time-dependent AUC, IPCW Brier score, external recalibration metrics, delta C-index, likelihood-ratio tests, and exploratory decision-curve data shift the paper toward prediction performance and deployability rather than ordinary differential-expression bioinformatics.",
        "- Dataset reuse, platform stratification, and endpoint audits make the evidence hierarchy explicit.",
        "",
        "## 4. Manual checks still needed",
        "",
        "- Confirm original code/calculator availability for each primary candidate.",
        "- Manually verify expression-scale and preprocessing details for numeric-cutoff signatures beyond what could be inferred from extracted fields.",
        "- Decide whether A131 should be called an L5 deployability candidate or an L4 cutoff-transportable candidate with deployability gaps, depending on final manuscript conservatism.",
        "- Validate final target-journal formatting and whether decision-curve outputs should remain supplementary/exploratory.",
        "",
        "## 5. Recommended main figures and supplementary tables",
        "",
        "- Figure 1: evidence ladder workflow from article inventory to deployability candidate.",
        "- Figure 2: reporting/reproducibility landscape.",
        "- Figure 3: all-57 external score-portability landscape.",
        "- Figure 4: cutoff transportability failure anatomy.",
        "- Figure 5: primary candidate prediction performance.",
        "- Figure 6: clinical deployability gap heatmap.",
        "- Supplementary tables: full evidence ladder, cutoff audit, reporting quality scores, all-57 GEO landscape, clinical deployability checklist, prediction-performance metrics, dataset reuse, platform stratification, endpoint consistency.",
        "",
        "## Ladder counts",
        "",
    ]
    for _, r in ladder_sig.iterrows():
        report.append(f"- {r['Evidence_ladder_label']}: {int(r['Count'])} signatures")
    report += [
        "",
        "Final positioning sentence for summaries and manuscript text: This is not another LUAD prognostic model. It is a precision-oncology deployability audit of published transcriptomic prognostic signatures.",
    ]
    (OUT / "npj_precision_oncology_readiness_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    data = build_metadata()
    art_ladder, sig_ladder = build_evidence_ladder(data)
    build_cutoff_audit(data, sig_ladder)
    build_reporting_quality(data, art_ladder, sig_ladder)
    build_all57_landscape(data, sig_ladder)
    primary_sigs = list(
        data["selection"].loc[
            data["selection"]["Final_analysis_tier"].isin(["Tier_1_strict_reported_cutoff", "Tier_2_primary_median_fallback"]),
            "Signature_ID",
        ]
    )
    build_prediction_performance(data, primary_sigs)
    build_deployability_checklist(data, sig_ladder)
    build_dataset_reuse(data)
    build_platform_analysis(data)
    build_endpoint_audit(data, sig_ladder)
    build_figure_data()
    build_markdown_outputs(data)
    produced = sorted(p.name for p in OUT.glob("*") if p.is_file())
    print(f"Precision-oncology deployability audit complete: {OUT}")
    print(f"Top-level files: {len(produced)}")
    for name in produced:
        print(name)


if __name__ == "__main__":
    main()
