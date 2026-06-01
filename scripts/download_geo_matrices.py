import csv
import gzip
import io
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
DATA = ROOT / "data"
RAW = DATA / "raw" / "geo"
PROC = DATA / "processed" / "geo"
MANIFEST = WORK / "luad_cohort_download_manifest.csv"
OUT_SUMMARY = WORK / "geo_matrix_download_summary.csv"
OUT_JSON = WORK / "geo_matrix_download_summary.json"


def geo_series_dir(gse: str) -> str:
    n = int(re.sub(r"\D", "", gse))
    prefix = f"GSE{n // 1000}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/matrix/"


def read_url_text(url: str, timeout: int = 60) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def list_matrix_urls(gse: str) -> list[str]:
    base = geo_series_dir(gse)
    try:
        html = read_url_text(base)
    except Exception:
        return [f"{base}{gse}_series_matrix.txt.gz"]
    urls = []
    for href in re.findall(r'href="([^"]+series_matrix\.txt\.gz)"', html, flags=re.I):
        if href.startswith("http"):
            urls.append(href)
        else:
            urls.append(base + href)
    if not urls:
        urls = [f"{base}{gse}_series_matrix.txt.gz"]
    return sorted(set(urls))


def download(url: str, dest: Path) -> tuple[str, int]:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return "cached", dest.stat().st_size
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=180) as resp, tmp.open("wb") as f:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)
    tmp.replace(dest)
    return "downloaded", dest.stat().st_size


def parse_series_matrix(path: Path, out_prefix: Path) -> dict:
    meta: dict[str, list[str]] = {}
    table_lines = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            if in_table:
                table_lines.append(line)
                continue
            if line.startswith("!Sample_") or line.startswith("!Series_") or line.startswith("!Platform_"):
                parts = next(csv.reader([line], delimiter="\t"))
                key = parts[0].lstrip("!")
                vals = [x.strip('"') for x in parts[1:]]
                meta.setdefault(key, vals)

    if not table_lines:
        raise ValueError("no series matrix table found")

    expr = pd.read_csv(io.StringIO("\n".join(table_lines)), sep="\t", dtype=str)
    expr = expr.rename(columns={expr.columns[0]: "ID_REF"})
    sample_cols = [c for c in expr.columns if c != "ID_REF"]
    for c in sample_cols:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    expr_path = out_prefix.with_suffix(".expression.tsv.gz")
    pheno_path = out_prefix.with_suffix(".pheno.tsv")
    meta_path = out_prefix.with_suffix(".metadata.json")
    if expr_path.exists() and pheno_path.exists() and meta_path.exists() and expr_path.stat().st_size > 0:
        # Still compute stats below from the current parse; this keeps summary fresh.
        pass
    else:
        tmp_expr = expr_path.with_suffix(expr_path.suffix + ".part")
        expr.to_csv(tmp_expr, sep="\t", index=False, compression="gzip")
        tmp_expr.replace(expr_path)

    accessions = meta.get("Sample_geo_accession", sample_cols)
    max_len = max([len(accessions), len(sample_cols)] + [len(v) for k, v in meta.items() if k.startswith("Sample_")] or [0])
    pheno = pd.DataFrame({"sample_column": sample_cols})
    if len(accessions) == len(sample_cols):
        pheno["geo_accession"] = accessions
    for key, vals in meta.items():
        if not key.startswith("Sample_") or len(vals) != len(sample_cols):
            continue
        col = key.replace("Sample_", "").lower()
        if col in {"geo_accession"}:
            continue
        pheno[col] = vals
    if not pheno_path.exists() or pheno_path.stat().st_size == 0:
        tmp_pheno = pheno_path.with_suffix(pheno_path.suffix + ".part")
        pheno.to_csv(tmp_pheno, sep="\t", index=False)
        tmp_pheno.replace(pheno_path)
    if not meta_path.exists() or meta_path.stat().st_size == 0:
        tmp_meta = meta_path.with_suffix(meta_path.suffix + ".part")
        tmp_meta.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
        tmp_meta.replace(meta_path)

    values = expr[sample_cols].to_numpy(dtype=float)
    finite = values[np.isfinite(values)]
    stats = {}
    if finite.size:
        qs = np.nanpercentile(finite, [0, 1, 25, 50, 75, 99, 100])
        stats = {
            "min": float(qs[0]),
            "p01": float(qs[1]),
            "p25": float(qs[2]),
            "median": float(qs[3]),
            "p75": float(qs[4]),
            "p99": float(qs[5]),
            "max": float(qs[6]),
            "missing_fraction": float(np.isnan(values).mean()),
        }
    max_val = stats.get("max", np.nan)
    p99 = stats.get("p99", np.nan)
    if np.isfinite(max_val) and max_val <= 30 and np.isfinite(stats.get("min", np.nan)) and stats.get("min", 0) >= -5:
        scale_guess = "likely log2 expression"
    elif np.isfinite(p99) and p99 > 100:
        scale_guess = "likely unlogged intensity/count-like"
    else:
        scale_guess = "unclear"

    platform_ids = "; ".join(sorted(set(meta.get("Series_platform_id", []) + meta.get("Sample_platform_id", []))))
    return {
        "matrix_file": str(path),
        "expression_file": str(expr_path),
        "phenotype_file": str(pheno_path),
        "metadata_file": str(meta_path),
        "platform_ids": platform_ids,
        "sample_count": len(sample_cols),
        "feature_count": len(expr),
        "scale_guess": scale_guess,
        **{f"expr_{k}": v for k, v in stats.items()},
    }


def select_cohorts(limit: int | None) -> list[str]:
    m = pd.read_csv(MANIFEST)
    geo = m[(m["Cohort_type"].eq("GEO")) & (m["Public_download_priority"].eq("yes"))].copy()
    geo["score"] = geo["Training_count"].fillna(0) * 3 + geo["Article_validation_count"].fillna(0) * 2 + geo["Eligible_external_count"].fillna(0) * 0.01
    geo = geo.sort_values(["score", "Training_count", "Article_validation_count"], ascending=False)
    cohorts = geo["Cohort"].tolist()
    if limit:
        cohorts = cohorts[:limit]
    return cohorts


def main():
    limit = 10
    explicit = None
    if len(sys.argv) > 1:
        if sys.argv[1] == "--all":
            limit = None
        elif re.fullmatch(r"\d+", sys.argv[1]):
            limit = int(sys.argv[1])
        else:
            explicit = [x.strip().upper() for x in sys.argv[1:] if x.strip()]
    RAW.mkdir(parents=True, exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)
    rows = []
    for gse in (explicit or select_cohorts(limit)):
        urls = list_matrix_urls(gse)
        for i, url in enumerate(urls, 1):
            suffix = "" if len(urls) == 1 else f"_{i}"
            raw_path = RAW / gse / Path(url).name
            row = {
                "GSE": gse,
                "matrix_url": url,
                "status": "",
                "download_status": "",
                "raw_file": str(raw_path),
            }
            try:
                dl_status, size = download(url, raw_path)
                row["download_status"] = dl_status
                row["raw_size_bytes"] = size
                parsed = parse_series_matrix(raw_path, PROC / gse / f"{gse}{suffix}")
                row.update(parsed)
                row["status"] = "ok"
            except Exception as exc:
                row["status"] = "error"
                row["error"] = repr(exc)
            rows.append(row)
            print(f"{gse}: {row['status']} {row.get('sample_count', '')} samples {row.get('feature_count', '')} features")
            time.sleep(0.2)

    fieldnames = sorted({k for row in rows for k in row})
    with OUT_SUMMARY.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    OUT_JSON.write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
