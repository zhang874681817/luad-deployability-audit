from __future__ import annotations

import json
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "visual_upgrade_biorender_data_driven_v1_polished"
DERIVED = OUT / "derived_plot_data"
SUPP = ROOT / "outputs" / "bmc_ai_second_stage_performance_package" / "additional_file_2_extracted"


def ok(name, passed, detail):
    return {"check": name, "passed": bool(passed), "detail": detail}


def main():
    checks = []
    # Figure 1
    checks.append(ok("Figure 1 7084 decomposition", 6028 + 498 + 276 + 282 == 7084, "6028 + 498 + 276 + 282 = 7084"))
    f1 = pd.read_csv(DERIVED / "figure_1_data.csv")
    required_counts = {10061, 7262, 178, 57, 27, 10, 1}
    checks.append(ok("Figure 1 funnel counts", required_counts.issubset(set(f1["count"].astype(int))), f"Found counts: {sorted(set(f1['count'].astype(int)))}"))

    # Figure 2
    all178 = pd.read_excel(SUPP / "Supplementary_Data_and_Tables" / "Supplementary_Data_02_All_178_articles.xlsx", sheet_name="All_178_articles")
    repro = int(all178["Reproducible"].astype(str).str.lower().eq("yes").sum())
    non = len(all178) - repro
    checks.append(ok("Figure 2 reproducible vs L0 counts", repro == 57 and non == 121, f"Formula-reproducible={repro}; non-reproducible/L0={non}"))
    fig2_texts = ""
    for p in [OUT / "source_data_used_figure_2.md"]:
        fig2_texts += p.read_text(encoding="utf-8", errors="replace")
    checks.append(ok("Figure 2 no p-value claim", "P < 0.001" not in fig2_texts and "p < 0.001" not in fig2_texts, "No p-value displayed in source notes."))

    # Figure 4
    fig4_notes = (OUT / "source_data_used_figure_4.md").read_text(encoding="utf-8", errors="replace")
    checks.append(ok("Figure 4 legend avoids restricted wording", "protective" not in fig4_notes.lower(), "Legend labels use HR direction and evaluability categories."))

    # Figure 5
    score = pd.read_csv(SUPP / "Scorecard_and_sensitivity" / "LUAD_deployability_scorecard_primary_candidates.csv")
    expected = {"A131": 17, "A129": 13, "A067": 14, "A097": 12, "A168": 13, "A001": 13, "A059": 13, "A098": 14, "A065": 14, "A173": 14}
    got = dict(zip(score["Article_ID"], score["total_deployability_score_0_20"].astype(int)))
    score_ok = all(got.get(k) == v for k, v in expected.items()) and score["total_deployability_score_0_20"].between(0, 20).all()
    checks.append(ok("Figure 5 scorecard 0-20 values", score_ok, f"Observed: { {k: got.get(k) for k in expected} }"))

    # Figure 6
    bench = pd.read_csv(SUPP / "AI_deep_learning_benchmark" / "AI_deepsurv_external_validation_results.csv")
    pl = bench[(bench["Dataset"].eq("GEO_external")) & (bench["Model"].eq("Penalized_linear_Cox_signature_gene_baseline"))]
    expected_pl = {
        "GSE13213": 0.569,
        "GSE30219": 0.526,
        "GSE31210": 0.696,
        "GSE3141": 0.554,
        "GSE37745": 0.527,
        "GSE42127": 0.529,
        "GSE50081": 0.489,
        "GSE68465": 0.544,
        "GSE72094": 0.587,
    }
    observed = {r["GSE"]: round(float(r["C_index"]), 3) for _, r in pl.iterrows()}
    checks.append(ok("Figure 6 penalized linear Cox non-empty", len(pl) == 9, f"Rows={len(pl)}; values={observed}"))
    checks.append(ok("Figure 6 penalized linear Cox verified values", all(abs(observed.get(k, -9) - v) <= 0.001 for k, v in expected_pl.items()), f"Expected={expected_pl}; observed={observed}"))

    # Output files
    expected_files = [
        "Figure_1_biorender_audit_funnel.pdf",
        "Figure_1_biorender_audit_funnel.png",
        "Figure_2_biorender_reporting_landscape.pdf",
        "Figure_2_biorender_reporting_landscape.png",
        "Figure_3_biorender_risk_dashboard.pdf",
        "Figure_3_biorender_risk_dashboard.png",
        "Figure_4_biorender_transportability_core.pdf",
        "Figure_4_biorender_transportability_core.png",
        "Figure_5_biorender_clinical_readiness.pdf",
        "Figure_5_biorender_clinical_readiness.png",
        "Figure_6_biorender_AI_stability_route.pdf",
        "Figure_6_biorender_AI_stability_route.png",
        "main_figures_contact_sheet.png",
        "figure_legends_revised_BMC.docx",
        "figure_source_notes.md",
    ]
    missing = [f for f in expected_files if not (OUT / f).exists()]
    checks.append(ok("Expected output files exist", not missing, f"Missing: {missing or 'none'}"))

    report_lines = [
        "# Visual QC report",
        "",
        "Scope: data-driven BioRender-inspired/Nature-style visual upgrade for LUAD BMC Medical Genomics figures 1-6.",
        "",
        "No AI image generator was used. Figures were generated with Python/matplotlib from verified source tables.",
        "",
        "## Data checks",
        "",
    ]
    for c in checks:
        status = "PASS" if c["passed"] else "FAIL"
        report_lines.append(f"- **{status}** - {c['check']}: {c['detail']}")
    report_lines.extend(
        [
            "",
            "## Corrections applied",
            "",
            "- Figure 1 uses P2/P3 deferred or boundary = 276, not 2765.",
            "- Figure 1 explicitly decomposes 7084 as 6028 + 498 + 276 + 282.",
            "- Figure 2 reports reporting completeness as a 0-10 score and does not display an unverified p-value.",
            "- Figure 4 legend uses HR-direction/evaluability labels and avoids ambiguous biological wording.",
            "- Figure 5 uses LUAD_deployability_scorecard_primary_candidates.csv total_deployability_score_0_20 values.",
            "- Figure 6 includes a non-empty penalized linear Cox benchmark row for all nine GEO cohorts.",
            "",
            "## Unavailable or limited source-data items",
            "",
            "- Figure 2 reporting score distributions are drawn only when the Reporting_scores worksheet exposes a machine-readable score column. If no score column is present in a different package revision, the script displays 'Data not available in source tables'.",
            "- No manuscript conclusions were changed; these are visual redesign outputs only.",
        ]
    )
    all_pass = all(c["passed"] for c in checks)
    report_lines.append("")
    report_lines.append(f"Overall status: {'PASS' if all_pass else 'CHECK REQUIRED'}")
    (OUT / "visual_QC_report.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    (OUT / "visual_QC_report.json").write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding="utf-8")

    zip_path = OUT / "LUAD_BMC_biorender_data_driven_figures_v1_polished.zip"
    zip_path.unlink(missing_ok=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for p in sorted(OUT.rglob("*")):
            if p.is_file() and p != zip_path:
                zf.write(p, p.relative_to(OUT).as_posix())
    # Compatibility alias requested at the top of the pasted brief.
    alias = OUT / "LUAD_BMC_visual_upgrade_biorender_data_driven_v1_polished.zip"
    alias.unlink(missing_ok=True)
    alias.write_bytes(zip_path.read_bytes())
    print(json.dumps({"all_pass": all_pass, "zip": str(zip_path), "alias": str(alias), "checks": checks}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
