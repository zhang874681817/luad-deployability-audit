import csv
import gzip
import io
import json
import re
import urllib.request
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
WORK = ROOT / "work"
PROC_GEO = DATA / "processed" / "geo"
RAW_GPL = DATA / "raw" / "geo_platforms"
PROC_GPL = DATA / "processed" / "geo_platforms"
OUT_SUMMARY = WORK / "geo_platform_annotation_summary.csv"


def gpl_dir(gpl: str) -> str:
    n = int(re.sub(r"\D", "", gpl))
    prefix = "GPLnnn" if n < 1000 else f"GPL{n // 1000}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl}/annot/"


def list_annot_urls(gpl: str) -> list[str]:
    base = gpl_dir(gpl)
    try:
        with urllib.request.urlopen(base, timeout=60) as resp:
            html = resp.read().decode("utf-8", errors="replace")
        urls = []
        for href in re.findall(r'href="([^"]+\.annot\.gz)"', html, flags=re.I):
            urls.append(href if href.startswith("http") else base + href)
        return sorted(set(urls)) or [base + f"{gpl}.annot.gz"]
    except Exception:
        return [base + f"{gpl}.annot.gz"]


def soft_url(gpl: str) -> str:
    n = int(re.sub(r"\D", "", gpl))
    prefix = "GPLnnn" if n < 1000 else f"GPL{n // 1000}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl}/soft/{gpl}_family.soft.gz"


def download(url: str, dest: Path):
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


def find_platforms() -> list[str]:
    platforms = set()
    for meta_path in PROC_GEO.glob("GSE*/GSE*.metadata.json"):
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        for key in ["Series_platform_id", "Sample_platform_id"]:
            for val in meta.get(key, []):
                if re.fullmatch(r"GPL\d+", str(val)):
                    platforms.add(str(val))
    return sorted(platforms, key=lambda x: int(x[3:]))


def read_annot_table(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    start = None
    for i, line in enumerate(lines):
        if line.startswith("ID\t") or line.startswith("ID "):
            start = i
            break
    if start is None:
        for i, line in enumerate(lines):
            if "Gene Symbol" in line and "\t" in line:
                start = i
                break
    if start is None:
        raise ValueError("annotation table header not found")
    table = "".join(lines[start:])
    df = pd.read_csv(io.StringIO(table), sep="\t", dtype=str, comment="#").fillna("")
    return df


def read_soft_platform_table(path: Path) -> pd.DataFrame:
    rows = []
    header = None
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "!platform_table_begin":
                in_table = True
                header = None
                continue
            if line == "!platform_table_end":
                break
            if not in_table:
                continue
            if header is None:
                header = next(csv.reader([line], delimiter="\t"))
            else:
                vals = next(csv.reader([line], delimiter="\t"))
                if len(vals) < len(header):
                    vals += [""] * (len(header) - len(vals))
                rows.append(vals[: len(header)])
    if header is None:
        raise ValueError("SOFT platform table header not found")
    return pd.DataFrame(rows, columns=header).fillna("")


def pick_symbol_column(df: pd.DataFrame) -> str:
    candidates = [
        "Gene Symbol",
        "GENE_SYMBOL",
        "Symbol",
        "gene_assignment",
        "Gene symbol",
        "gene symbol",
        "ORF",
    ]
    for c in candidates:
        if c in df.columns:
            return c
    for c in df.columns:
        if "symbol" in c.lower() or "gene" in c.lower():
            return c
    return ""


def split_symbols(value: str) -> list[str]:
    value = str(value or "")
    value = re.sub(r"///", "///", value)
    parts = re.split(r"\s*///\s*|\s*//\s*|\s*;\s*|\s*,\s*|\s+\|\s+", value)
    out = []
    for p in parts:
        p = p.strip()
        if not p or p in {"---", "NA", "nan", "null"}:
            continue
        # Affymetrix annotations often put "GENE // description"; keep the symbol-ish token.
        token = p.split(" ")[0].strip()
        if re.fullmatch(r"[A-Za-z0-9_.-]+", token):
            out.append(token)
    return list(dict.fromkeys(out))


def process_platform(gpl: str) -> dict:
    urls = list_annot_urls(gpl)
    url = urls[0]
    raw = RAW_GPL / gpl / Path(url).name
    row = {"GPL": gpl, "annot_url": url, "raw_file": str(raw)}
    try:
        status, size = download(url, raw)
        row["download_status"] = status
        row["raw_size_bytes"] = size
        df = read_annot_table(raw)
        row["annotation_source"] = "annot.gz"
    except Exception as annot_exc:
        url = soft_url(gpl)
        raw = RAW_GPL / gpl / Path(url).name
        row["annot_url"] = url
        row["raw_file"] = str(raw)
        status, size = download(url, raw)
        row["download_status"] = status
        row["raw_size_bytes"] = size
        df = read_soft_platform_table(raw)
        row["annotation_source"] = f"family.soft.gz fallback; annot_error={annot_exc!r}"
    sym_col = pick_symbol_column(df)
    row["columns"] = "; ".join(df.columns[:30])
    row["symbol_column"] = sym_col
    row["probe_count"] = len(df)
    if "ID" not in df.columns or not sym_col:
        row["status"] = "no_symbol_column"
        return row
    mapping = []
    for _, r in df.iterrows():
        probe = r.get("ID", "")
        symbols = split_symbols(r.get(sym_col, ""))
        if not symbols:
            mapping.append({"ID_REF": probe, "Gene_symbol": "", "Mapping_status": "no_symbol"})
        else:
            for sym in symbols:
                mapping.append({"ID_REF": probe, "Gene_symbol": sym, "Mapping_status": "mapped"})
    out_dir = PROC_GPL / gpl
    out_dir.mkdir(parents=True, exist_ok=True)
    map_path = out_dir / f"{gpl}.probe_to_symbol.tsv"
    pd.DataFrame(mapping).to_csv(map_path, sep="\t", index=False)
    row["mapping_file"] = str(map_path)
    row["mapped_probe_rows"] = sum(1 for m in mapping if m["Mapping_status"] == "mapped")
    row["unique_symbols"] = len({m["Gene_symbol"] for m in mapping if m["Gene_symbol"]})
    row["status"] = "ok"
    return row


def main():
    rows = []
    for gpl in find_platforms():
        try:
            row = process_platform(gpl)
        except Exception as exc:
            row = {"GPL": gpl, "status": "error", "error": repr(exc)}
        rows.append(row)
        print(gpl, row.get("status"), row.get("probe_count", ""))
    fields = sorted({k for r in rows for k in r})
    with OUT_SUMMARY.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
