import csv
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
PROC_GEO = ROOT / "data" / "processed" / "geo"
CLIN = WORK / "geo_clinical_samples.csv"
CLIN_SUMMARY = WORK / "geo_clinical_summary.csv"
EXPR_SUMMARY = WORK / "geo_gene_expression_summary.csv"
OUT = WORK / "geo_expression_clinical_readiness.csv"


def read_expr_samples(gse: str) -> set[str]:
    files = list((PROC_GEO / gse).glob("*.gene_symbol_expression.for_scoring.tsv.gz"))
    if not files:
        return set()
    df = pd.read_csv(files[0], sep="\t", nrows=0)
    return set(c for c in df.columns if c != "Gene_symbol")


def is_luad_like(row) -> bool:
    text = " ".join(str(row.get(c, "")) for c in ["histology", "sample_type", "stage"])
    text_l = text.lower()
    if any(x in text_l for x in ["normal", "squamous", "large cell", "ntl"]):
        return False
    if any(x in text_l for x in ["adeno", "adenocarcinoma", "lung tumour", "lung tumor", "primary lung tumor", "nsclc"]):
        return True
    return False


def main():
    clin = pd.read_csv(CLIN)
    csum = pd.read_csv(CLIN_SUMMARY)
    esum = pd.read_csv(EXPR_SUMMARY)
    rows = []
    for _, er in esum.iterrows():
        gse = er["GSE"]
        expr_samples = read_expr_samples(gse)
        sub = clin[clin["GSE"].eq(gse)].copy()
        sub["has_expression"] = sub["geo_accession"].isin(expr_samples)
        sub["OS_complete"] = sub["OS_time_months"].notna() & sub["OS_event"].notna()
        sub["LUAD_like"] = sub.apply(is_luad_like, axis=1)
        expr_clin = sub[sub["has_expression"]]
        expr_os = expr_clin[expr_clin["OS_complete"]]
        expr_os_luad = expr_os[expr_os["LUAD_like"]]
        rows.append(
            {
                "GSE": gse,
                "Matrix_status": er.get("status", ""),
                "GPL": er.get("GPL", ""),
                "Expression_sample_count": len(expr_samples),
                "Clinical_sample_count": len(sub),
                "Expression_clinical_overlap_n": int(expr_clin.shape[0]),
                "OS_complete_overlap_n": int(expr_os.shape[0]),
                "OS_event_overlap_n": int(expr_os["OS_event"].sum()) if len(expr_os) else 0,
                "LUAD_like_OS_complete_n": int(expr_os_luad.shape[0]),
                "LUAD_like_OS_event_n": int(expr_os_luad["OS_event"].sum()) if len(expr_os_luad) else 0,
                "Machine_ready_all_histologies": "yes" if expr_os.shape[0] >= 30 and expr_os["OS_event"].sum() >= 5 and er.get("status", "") == "ok" else "no",
                "Machine_ready_LUAD_like": "yes" if expr_os_luad.shape[0] >= 30 and expr_os_luad["OS_event"].sum() >= 5 and er.get("status", "") == "ok" else "no",
                "Scale_guess": er.get("scale_guess", ""),
                "Scoring_transform": er.get("scoring_transform", ""),
                "Scoring_expression_file": er.get("scoring_expression_file", ""),
                "Clinical_OS_time_source": csum.loc[csum["GSE"].eq(gse), "OS_time_source"].iloc[0] if csum["GSE"].eq(gse).any() else "",
                "Clinical_OS_event_source": csum.loc[csum["GSE"].eq(gse), "OS_event_source"].iloc[0] if csum["GSE"].eq(gse).any() else "",
            }
        )
    pd.DataFrame(rows).to_csv(OUT, index=False, encoding="utf-8-sig")
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
