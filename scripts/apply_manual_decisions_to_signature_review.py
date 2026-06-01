from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
OUT = WORK / "geo_external_validation_manuscript_signature_review.manual_decisions.csv"


def main():
    review = pd.read_csv(WORK / "geo_external_validation_manuscript_signature_review.csv")
    review["Manual_decision"] = ""
    review["Manual_decision_reason"] = ""

    mask = review["Signature_ID"].eq("A097_PMID36971680")
    review.loc[mask, "Manual_decision"] = "Keep_as_median_fallback_score_portability"
    review.loc[mask, "Manual_decision_reason"] = (
        "Uploaded PDF confirms TCGA/training risk-score cutoff 1.66, but the paper does not provide GEO OS "
        "risk-score validation and applying 1.66 to processed GEO matrices collapses all scored cohorts to one group."
    )

    mask = review["Signature_ID"].eq("A010_PMID40591077")
    review.loc[mask, "Manual_decision"] = "Sensitivity_only_after_formula_override"
    review.loc[mask, "Manual_decision_reason"] = (
        "Original machine extraction missed negative terms from the 14-gene formula; manual override makes it testable, "
        "but grouped validation still uses cohort median fallback."
    )

    mask = review["Signature_ID"].eq("A104_PMID34760888")
    review.loc[mask, "Manual_decision"] = "Sensitivity_only_after_symbol_override"
    review.loc[mask, "Manual_decision_reason"] = (
        "Manual aliases GARS/GARS1, EFG1/GFM1, and CABC1/COQ8A make it testable; keep as sensitivity analysis because "
        "the override is manual and grouped validation uses cohort median fallback."
    )

    review.to_csv(OUT, index=False, encoding="utf-8")
    print(f"Wrote {OUT} ({len(review)} rows)")
    print(
        review[review["Manual_decision"].ne("")]
        [["Signature_ID", "Recommendation_Tier", "Manual_decision", "Manual_decision_reason"]]
        .to_string(index=False, max_colwidth=120)
    )


if __name__ == "__main__":
    main()
