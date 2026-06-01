from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
OUT = WORK / "geo_external_validation_manual_override_comparison.csv"


def load(path):
    df = pd.read_csv(path)
    for col in ["Tested_cohorts", "Significant_cohorts_p05", "Median_C_index", "Best_C_index", "Median_group_HR"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


base = load(WORK / "geo_external_validation_signature_summary.csv")
override = load(WORK / "geo_external_validation_signature_summary.manual_overrides.csv")

cols = ["Signature_ID", "Tested_cohorts", "Significant_cohorts_p05", "Median_C_index", "Best_C_index", "Median_group_HR"]
merged = override[cols].merge(base[cols], on="Signature_ID", how="left", suffixes=("_override", "_base"))

meta = pd.read_csv(WORK / "luad_signature_readiness.csv")[
    ["Signature_ID", "Article_ID", "PMID", "Title", "Cutoff_type", "Cutoff_value", "Training_cohort", "Validation_cohort"]
]
merged = meta.merge(merged, on="Signature_ID", how="right")

for col in ["Tested_cohorts", "Significant_cohorts_p05", "Median_C_index", "Best_C_index", "Median_group_HR"]:
    merged[f"Delta_{col}"] = merged[f"{col}_override"].fillna(0) - merged[f"{col}_base"].fillna(0)

merged["Manual_override_applied"] = merged["Signature_ID"].isin(["A010_PMID40591077", "A104_PMID34760888"]).map({True: "yes", False: "no"})
merged["Interpretation"] = ""
merged.loc[merged["Signature_ID"].eq("A010_PMID40591077"), "Interpretation"] = (
    "Full 14-gene formula and OGDH/TMEM5 aliases make the signature testable; use only as median-fallback sensitivity analysis."
)
merged.loc[merged["Signature_ID"].eq("A104_PMID34760888"), "Interpretation"] = (
    "GARS/GARS1 and EFG1/GFM1 aliases make the signature testable; signal is promising but relies on manual symbol override and median fallback."
)

merged = merged.sort_values(
    ["Manual_override_applied", "Significant_cohorts_p05_override", "Tested_cohorts_override", "Median_C_index_override"],
    ascending=[False, False, False, False],
)
merged.to_csv(OUT, index=False, encoding="utf-8")
print(f"Wrote {OUT} ({len(merged)} rows)")
