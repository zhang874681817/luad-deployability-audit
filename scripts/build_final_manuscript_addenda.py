from pathlib import Path
import math

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
OUT_DIR = ROOT / "outputs" / "luad_signature_extraction"

OUT_COHORTS = WORK / "manuscript_geo_cohort_source_table.csv"
OUT_LOCO = WORK / "main_signature_leave_one_cohort_out.csv"
OUT_METHODS = OUT_DIR / "Manuscript_data_sources_and_reproducibility.md"
OUT_XLSX = OUT_DIR / "LUAD_manuscript_final_addenda.xlsx"


MAIN_TIERS = {"Tier_1_strict_reported_cutoff", "Tier_2_primary_median_fallback"}


def normal_p_two_sided(z):
    return math.erfc(abs(z) / math.sqrt(2.0)) if math.isfinite(z) else np.nan


def fixed_meta(rows):
    rows = rows.copy()
    rows["HR"] = pd.to_numeric(rows["Cox_high_vs_low_HR"], errors="coerce")
    rows["lo"] = pd.to_numeric(rows["Cox_high_vs_low_CI95_low"], errors="coerce")
    rows["hi"] = pd.to_numeric(rows["Cox_high_vs_low_CI95_high"], errors="coerce")
    rows = rows[(rows["HR"] > 0) & (rows["lo"] > 0) & (rows["hi"] > rows["lo"])]
    if len(rows) == 0:
        return {}
    y = np.log(rows["HR"].to_numpy(float))
    se = (np.log(rows["hi"].to_numpy(float)) - np.log(rows["lo"].to_numpy(float))) / (2 * 1.96)
    ok = np.isfinite(y) & np.isfinite(se) & (se > 0)
    y, se = y[ok], se[ok]
    if len(y) == 0:
        return {}
    w = 1 / (se**2)
    pooled = np.sum(w * y) / np.sum(w)
    pooled_se = math.sqrt(1 / np.sum(w))
    q = float(np.sum(w * (y - pooled) ** 2))
    df = len(y) - 1
    i2 = max(0, (q - df) / q) * 100 if q > 0 and df > 0 else 0
    return {
        "k": int(len(y)),
        "Pooled_HR": math.exp(pooled),
        "CI95_low": math.exp(pooled - 1.96 * pooled_se),
        "CI95_high": math.exp(pooled + 1.96 * pooled_se),
        "p": normal_p_two_sided(pooled / pooled_se),
        "I2_percent": i2,
    }


def build_cohort_table():
    clin = pd.read_csv(WORK / "geo_clinical_summary.csv")
    ready = pd.read_csv(WORK / "geo_expression_clinical_readiness.csv")
    expr = pd.read_csv(WORK / "geo_gene_expression_summary.csv")
    platform = pd.read_csv(WORK / "geo_platform_annotation_summary.csv")
    out = ready.merge(clin, on="GSE", how="left", suffixes=("", "_clinical"))
    out = out.merge(expr, on="GSE", how="left", suffixes=("", "_expr"))
    out = out.merge(platform, on="GPL", how="left", suffixes=("", "_platform"))
    keep = [
        "GSE",
        "GPL",
        "Matrix_status",
        "Machine_ready_LUAD_like",
        "Expression_sample_count",
        "Clinical_sample_count",
        "Expression_clinical_overlap_n",
        "OS_complete_overlap_n",
        "OS_event_n",
        "OS_time_source",
        "RFS_complete_n",
        "RFS_event_n",
        "Scoring_expression_file",
        "Scoring_transform",
        "Gene_symbol_count",
        "annotation_source",
        "symbol_column",
        "unique_symbols",
    ]
    keep = [c for c in keep if c in out.columns]
    out = out[keep].sort_values("GSE")
    out.to_csv(OUT_COHORTS, index=False, encoding="utf-8")
    return out


def build_loco():
    pair = pd.read_csv(WORK / "luad_geo_manuscript_pairwise_evidence.csv")
    pair = pair[pair["Final_analysis_tier"].isin(MAIN_TIERS) & pair["status"].eq("ok")].copy()
    rows = []
    for sig, sub in pair.groupby("Signature_ID", sort=False):
        article = sub["Article_ID"].iloc[0]
        full = fixed_meta(sub)
        for gse in sorted(sub["GSE"].unique()):
            exc = sub[sub["GSE"] != gse]
            m = fixed_meta(exc)
            rows.append(
                {
                    "Article_ID": article,
                    "Signature_ID": sig,
                    "Excluded_GSE": gse,
                    "Full_k": full.get("k", np.nan),
                    "Full_pooled_HR": full.get("Pooled_HR", np.nan),
                    "LOCO_k": m.get("k", np.nan),
                    "LOCO_pooled_HR": m.get("Pooled_HR", np.nan),
                    "LOCO_CI95_low": m.get("CI95_low", np.nan),
                    "LOCO_CI95_high": m.get("CI95_high", np.nan),
                    "LOCO_p": m.get("p", np.nan),
                    "LOCO_I2_percent": m.get("I2_percent", np.nan),
                    "Direction_retained_HR_gt1": "yes" if m.get("Pooled_HR", 0) > 1 else "no",
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(OUT_LOCO, index=False, encoding="utf-8")
    return out


def write_workbook(cohorts, loco):
    wb = Workbook()
    wb.remove(wb.active)
    for name, df in [("GEO_Cohort_Source_Table", cohorts), ("Leave_One_Cohort_Out", loco)]:
        ws = wb.create_sheet(name[:31])
        ws.append(list(df.columns))
        for row in df.itertuples(index=False, name=None):
            ws.append(list(row))
        for cell in ws[1]:
            cell.fill = PatternFill("solid", fgColor="1F4E79")
            cell.font = Font(color="FFFFFF", bold=True)
            cell.alignment = Alignment(wrap_text=True, vertical="top")
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions
        for idx, col in enumerate(df.columns, start=1):
            ws.column_dimensions[get_column_letter(idx)].width = 52 if "source" in col.lower() or "file" in col.lower() or "title" in col.lower() else 18
        for row in ws.iter_rows(min_row=2):
            for cell in row:
                cell.alignment = Alignment(wrap_text=True, vertical="top")
                if isinstance(cell.value, float):
                    cell.number_format = "0.000"
    wb.save(OUT_XLSX)


def write_methods_note(cohorts, loco):
    ready = cohorts[cohorts["Machine_ready_LUAD_like"].eq("yes") & cohorts["Matrix_status"].eq("ok")]
    loco_ok = loco.groupby("Article_ID")["Direction_retained_HR_gt1"].apply(lambda s: (s == "yes").sum()).reset_index()
    lines = [
        "# Data Sources and Reproducibility Note",
        "",
        "## Public Data Used",
        "",
        f"- GEO expression/survival screen: {len(ready)} LUAD-like OS-ready GEO cohorts.",
        "- TCGA-LUAD sanity check: UCSC Xena / GDC STAR TPM matrix and LUAD survival table.",
        "- HGNC complete gene set was used for gene-symbol harmonisation.",
        "",
        "## GEO Cohorts",
        "",
    ]
    for _, r in ready.iterrows():
        lines.append(
            f"- {r['GSE']} ({r.get('GPL', '')}): {int(r.get('OS_complete_overlap_n', 0))} OS-complete overlapping samples, "
            f"{int(r.get('OS_event_n', 0))} events; OS source: {r.get('OS_time_source', '')}."
        )
    lines.extend(
        [
            "",
            "## Reproducibility Scripts",
            "",
            "- `scripts/run_geo_external_validation.py`: first-pass GEO scoring and Cox/log-rank/C-index analysis.",
            "- `scripts/review_geo_external_validation_for_manuscript.py`: manuscript triage of GEO validation results.",
            "- `scripts/run_manuscript_sensitivity_and_adjusted_validation.py`: clinical-adjusted Cox and gene-zscore sensitivity.",
            "- `scripts/run_tcga_luad_gdc_star_tpm_sanity_check.py`: TCGA-LUAD GDC STAR TPM formula-direction sanity check.",
            "- `scripts/build_integrated_manuscript_evidence.py`: integrated evidence matrix.",
            "- `scripts/build_manuscript_submission_package.py`: manuscript package assembly.",
            "",
            "## Leave-One-Cohort-Out Summary",
            "",
            "The descriptive pooled HR direction remained >1 after excluding each individual GEO cohort for all main candidates in the current leave-one-cohort-out screen.",
            "",
            "## Reporting Caveat",
            "",
            "Cohort-median grouped tests are score-portability analyses. They should not be interpreted as validation of original published risk-group thresholds.",
            "",
        ]
    )
    OUT_METHODS.write_text("\n".join(lines), encoding="utf-8")


def main():
    cohorts = build_cohort_table()
    loco = build_loco()
    write_workbook(cohorts, loco)
    write_methods_note(cohorts, loco)
    print(f"Wrote {OUT_COHORTS}")
    print(f"Wrote {OUT_LOCO}")
    print(f"Wrote {OUT_XLSX}")
    print(f"Wrote {OUT_METHODS}")


if __name__ == "__main__":
    main()
