import csv
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
WORK = ROOT / "work"
PROC_GEO = DATA / "processed" / "geo"
PROC_GPL = DATA / "processed" / "geo_platforms"
ARTICLE_JSON = WORK / "article_extraction_qc.json"
PREP_JSON = WORK / "luad_reproduction_preparation.json"
HGNC_TSV = WORK / "hgnc_complete_set.txt"
OUT_GEO_SUMMARY = WORK / "geo_gene_expression_summary.csv"
OUT_COVERAGE = WORK / "signature_geo_gene_coverage.csv"


def norm(s):
    return re.sub(r"\s+", " ", str(s or "")).strip()


def signature_id(row):
    pmid = str(row.get("PMID", "") or "").replace(".0", "")
    return f"{row.get('Article_ID')}{('_PMID' + pmid) if pmid else ''}"


def split_symbols(v):
    return [x.strip() for x in str(v or "").split(";") if x.strip()]


def load_scale_guess(expr_file: Path) -> dict:
    expr = pd.read_csv(expr_file, sep="\t", nrows=5000)
    sample_cols = [c for c in expr.columns if c != "ID_REF"]
    vals = expr[sample_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return {"scale_guess": "no_numeric_values"}
    qs = np.nanpercentile(finite, [0, 50, 99, 100])
    if qs[3] <= 30 and qs[0] >= -10:
        scale = "likely log2 expression"
    elif qs[2] > 100:
        scale = "likely unlogged intensity/count-like"
    else:
        scale = "unclear"
    return {"scale_guess": scale, "preview_min": qs[0], "preview_median": qs[1], "preview_p99": qs[2], "preview_max": qs[3]}


def transform_for_scoring(df: pd.DataFrame, scale_guess: str) -> tuple[pd.DataFrame, str]:
    sample_cols = [c for c in df.columns if c != "Gene_symbol"]
    out = df.copy()
    if scale_guess == "likely unlogged intensity/count-like":
        out[sample_cols] = np.log2(out[sample_cols].clip(lower=0) + 1)
        return out, "log2(x+1) applied to unlogged-looking matrix"
    return out, "unchanged"


def collapse_one(gse_dir: Path) -> dict:
    gse = gse_dir.name
    expr_files = list(gse_dir.glob("*.expression.tsv.gz"))
    meta_files = list(gse_dir.glob("*.metadata.json"))
    if not expr_files or not meta_files:
        return {"GSE": gse, "status": "missing_expression_or_metadata"}
    expr_file = expr_files[0]
    if expr_file.stat().st_size <= 200:
        return {"GSE": gse, "status": "empty_expression", "expression_file": str(expr_file)}
    meta = json.loads(meta_files[0].read_text(encoding="utf-8"))
    platforms = sorted(set(meta.get("Series_platform_id", []) + meta.get("Sample_platform_id", [])))
    platforms = [p for p in platforms if re.fullmatch(r"GPL\d+", str(p))]
    if len(platforms) != 1:
        return {"GSE": gse, "status": "multi_or_missing_platform", "platforms": "; ".join(platforms), "expression_file": str(expr_file)}
    gpl = platforms[0]
    expr = pd.read_csv(expr_file, sep="\t")
    sample_cols = [c for c in expr.columns if c != "ID_REF"]
    for c in sample_cols:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")
    map_file = PROC_GPL / gpl / f"{gpl}.probe_to_symbol.tsv"
    mapping_note = ""
    if gpl == "GPL15048" and not map_file.exists():
        hgnc = pd.read_csv(HGNC_TSV, sep="\t", dtype=str, usecols=["symbol", "refseq_accession"]).fillna("")
        refseq = {}
        for _, r in hgnc.iterrows():
            for acc in str(r["refseq_accession"]).split("|"):
                acc = acc.split(".")[0].strip()
                if acc:
                    refseq.setdefault(acc, r["symbol"])
        dyn = []
        for probe in expr["ID_REF"].astype(str):
            m = re.search(r"((?:NM|NR|XM|XR)_\d+)", probe)
            if m and m.group(1) in refseq:
                dyn.append({"ID_REF": probe, "Gene_symbol": refseq[m.group(1)]})
        mapping = pd.DataFrame(dyn)
        mapping_note = "GPL15048 dynamic RefSeq-to-HGNC mapping from ID_REF; non-RefSeq probes unavailable"
    else:
        if not map_file.exists():
            return {"GSE": gse, "status": "missing_platform_mapping", "platforms": gpl, "expression_file": str(expr_file)}
        mapping = pd.read_csv(map_file, sep="\t", dtype=str).fillna("")
        mapping = mapping[mapping["Gene_symbol"].ne("")]
    if mapping.empty:
        return {"GSE": gse, "status": "empty_platform_mapping", "platforms": gpl, "expression_file": str(expr_file)}
    merged = expr.merge(mapping[["ID_REF", "Gene_symbol"]], on="ID_REF", how="inner")
    if merged.empty:
        return {"GSE": gse, "status": "no_probe_mapping_overlap", "platforms": gpl, "expression_file": str(expr_file)}
    gene = merged.groupby("Gene_symbol", as_index=False)[sample_cols].mean(numeric_only=True)
    scale = load_scale_guess(expr_file)
    scoring, transform = transform_for_scoring(gene, scale["scale_guess"])
    out_dir = gse_dir
    raw_gene_file = out_dir / f"{gse}.gene_symbol_expression.mean.tsv.gz"
    scoring_file = out_dir / f"{gse}.gene_symbol_expression.for_scoring.tsv.gz"
    gene.to_csv(raw_gene_file, sep="\t", index=False, compression="gzip")
    scoring.to_csv(scoring_file, sep="\t", index=False, compression="gzip")
    return {
        "GSE": gse,
        "status": "ok",
        "GPL": gpl,
        "Mapping_note": mapping_note,
        "sample_count": len(sample_cols),
        "probe_feature_count": len(expr),
        "mapped_probe_rows": len(merged),
        "gene_symbol_count": len(gene),
        "expression_file": str(expr_file),
        "gene_expression_file": str(raw_gene_file),
        "scoring_expression_file": str(scoring_file),
        "scoring_transform": transform,
        **scale,
    }


def load_signature_genes():
    payload = json.loads(PREP_JSON.read_text(encoding="utf-8"))
    by_sig = {}
    for row in payload["signature_gene_coefficients"]:
        sig = row["Signature_ID"]
        gene = row.get("Approved_symbol") or row.get("Gene_original")
        by_sig.setdefault(sig, {"metadata": row, "genes": []})
        if gene:
            by_sig[sig]["genes"].append(gene)
    for sig in by_sig:
        by_sig[sig]["genes"] = list(dict.fromkeys(by_sig[sig]["genes"]))
    return by_sig


def coverage_rows(gene_summaries):
    sigs = load_signature_genes()
    rows = []
    for gs in gene_summaries:
        if gs.get("status") != "ok":
            continue
        gse = gs["GSE"]
        gene_file = Path(gs["scoring_expression_file"])
        genes = set(pd.read_csv(gene_file, sep="\t", usecols=["Gene_symbol"])["Gene_symbol"].astype(str))
        for sig, obj in sigs.items():
            need = obj["genes"]
            present = [g for g in need if g in genes]
            missing = [g for g in need if g not in genes]
            meta = obj["metadata"]
            rows.append(
                {
                    "Signature_ID": sig,
                    "Article_ID": meta.get("Article_ID", ""),
                    "PMID": meta.get("PMID", ""),
                    "Title": meta.get("Title", ""),
                    "GSE": gse,
                    "GPL": gs.get("GPL", ""),
                    "Signature_gene_count": len(need),
                    "Present_gene_count": len(present),
                    "Missing_gene_count": len(missing),
                    "Coverage_fraction": round(len(present) / len(need), 4) if need else "",
                    "Present_genes": "; ".join(present),
                    "Missing_genes": "; ".join(missing),
                    "Scoring_expression_file": gs.get("scoring_expression_file", ""),
                    "Scoring_transform": gs.get("scoring_transform", ""),
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict]):
    fields = sorted({k for r in rows for k in r})
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main():
    summaries = []
    for gse_dir in sorted(PROC_GEO.glob("GSE*")):
        row = collapse_one(gse_dir)
        summaries.append(row)
        print(row["GSE"], row.get("status"), row.get("gene_symbol_count", ""))
    write_csv(OUT_GEO_SUMMARY, summaries)
    cov = coverage_rows(summaries)
    write_csv(OUT_COVERAGE, cov)
    print(f"wrote {OUT_GEO_SUMMARY}")
    print(f"wrote {OUT_COVERAGE}")


if __name__ == "__main__":
    main()
