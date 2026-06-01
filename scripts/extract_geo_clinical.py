import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
WORK = ROOT / "work"
RAW_GEO = DATA / "raw" / "geo"
PROC_GEO = DATA / "processed" / "geo"
OUT_SAMPLE = WORK / "geo_clinical_samples.csv"
OUT_SUMMARY = WORK / "geo_clinical_summary.csv"


def norm(s):
    return re.sub(r"\s+", " ", str(s or "")).strip()


def safe_col(s):
    s = str(s or "").lower()
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return s[:80] or "field"


def parse_series_sample_metadata(raw_gz: Path) -> pd.DataFrame:
    sample_keys = defaultdict(list)
    with gzip.open(raw_gz, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                break
            if not line.startswith("!Sample_"):
                continue
            parts = next(csv.reader([line], delimiter="\t"))
            key = parts[0].lstrip("!")
            vals = [x.strip('"') for x in parts[1:]]
            sample_keys[key].append(vals)

    accessions = []
    if "Sample_geo_accession" in sample_keys:
        accessions = sample_keys["Sample_geo_accession"][0]
    else:
        max_len = max((len(v[0]) for v in sample_keys.values() if v), default=0)
        accessions = [f"sample_{i+1}" for i in range(max_len)]
    df = pd.DataFrame({"geo_accession": accessions})

    for key, rows in sample_keys.items():
        if key == "Sample_geo_accession":
            continue
        if len(rows) == 1:
            vals = rows[0]
            if len(vals) == len(df):
                df[safe_col(key.replace("Sample_", ""))] = vals
            continue
        for idx, vals in enumerate(rows, 1):
            if len(vals) != len(df):
                continue
            col = safe_col(key.replace("Sample_", f"") + f"_{idx}")
            df[col] = vals
            # Also parse key-value sample characteristics into explicit columns.
            for i, val in enumerate(vals):
                m = re.match(r"\s*([^:]+?)\s*:\s*(.+)\s*$", str(val))
                if m:
                    k = safe_col(m.group(1))
                    if k not in df.columns:
                        df[k] = ""
                    df.at[i, k] = norm(m.group(2))
    return df


def to_float(v):
    if v is None:
        return np.nan
    s = str(v).strip()
    if not s or s.lower() in {"na", "n/a", "nan", "unknown", "not available", "null"}:
        return np.nan
    m = re.search(r"[-+]?\d*\.?\d+", s)
    if not m:
        return np.nan
    try:
        return float(m.group(0))
    except Exception:
        return np.nan


def event_from_value(v):
    s = str(v or "").strip().lower()
    if not s:
        return np.nan
    if re.fullmatch(r"[01](?:\.0)?", s):
        return int(float(s))
    if s in {"d", "dead", "deceased", "y", "yes"}:
        return 1
    if s in {"a", "alive", "n", "no"}:
        return 0
    if any(x in s for x in ["dead", "deceased", "death", "yes", "relapse", "recurrence", "recurred", "progressed"]):
        if not any(x in s for x in ["alive", "no death", "no recurrence", "no relapse"]):
            return 1
    if any(x in s for x in ["alive", "living", "censored", "disease free", "non-recurr"]) or re.search(r"\bno\b", s):
        return 0
    return np.nan


def extract_clinical(df: pd.DataFrame, gse: str) -> pd.DataFrame:
    out = pd.DataFrame({"GSE": gse, "geo_accession": df["geo_accession"]})
    cols = list(df.columns)
    lower = {c: c.lower() for c in cols}

    def find_cols(patterns):
        found = []
        for c in cols:
            lc = lower[c]
            if any(re.search(p, lc) for p in patterns):
                found.append(c)
        return found

    # Time in months where possible. These heuristics prefer OS-like fields over follow-up-only fields.
    os_time_cols = find_cols([
        r"\bos\b.*time",
        r"^survival_time$",
        r"surv.*month",
        r"survival.*month",
        r"survival.*days",
        r"survival_days",
        r"overall.*survival",
        r"follow.*month",
        r"follow_up_time",
        r"time.*death",
        r"days.*death",
        r"days_before_death",
        r"last_contact_or_death",
    ])
    # Keep generic "status" last; GEO also has a non-clinical publication status field.
    os_event_cols = find_cols([r"\bos\b.*event", r"vital", r"death", r"dead", r"survival_status", r"status.*dead", r"status_0_alive"])
    status_cols = [c for c in find_cols([r"^status$"]) if c not in os_event_cols]
    os_event_cols = os_event_cols + status_cols
    rfs_time_cols = find_cols([r"rfs", r"recurr.*time", r"relapse.*time", r"days.*recurrence", r"days_before_relapse", r"months_before_relapse", r"disease.*free"])
    rfs_event_cols = find_cols([r"recurr", r"relapse", r"rfs.*event"])

    def first_numeric(candidates):
        vals = pd.Series([np.nan] * len(df), dtype=float)
        source = ""
        for c in candidates:
            parsed = df[c].map(to_float)
            if parsed.notna().sum() > vals.notna().sum():
                vals = parsed
                source = c
        return vals, source

    def first_event(candidates):
        vals = pd.Series([np.nan] * len(df), dtype=float)
        source = ""
        for c in candidates:
            parsed = df[c].map(event_from_value)
            if parsed.notna().sum() > vals.notna().sum():
                vals = parsed
                source = c
        return vals, source

    os_time, os_time_src = first_numeric(os_time_cols)
    os_event, os_event_src = first_event(os_event_cols)
    rfs_time, rfs_time_src = first_numeric(rfs_time_cols)
    rfs_event, rfs_event_src = first_event(rfs_event_cols)

    # Dataset-specific cleanup for compact characteristics strings.
    if gse == "GSE3141":
        char_cols = [c for c in cols if c.startswith("characteristics_ch1")]
        combined = df[char_cols].astype(str).agg("; ".join, axis=1)
        os_time = combined.str.extract(r"Surv\(months\):\s*([-+]?\d*\.?\d+)")[0].map(to_float)
        os_event = combined.str.extract(r"STATUS\(0=alive,\s*1=dead\):\s*([01])")[0].map(event_from_value)
        os_time_src = "characteristics_ch1 Surv(months)"
        os_event_src = "characteristics_ch1 STATUS"

    # Convert day-based OS/RFS fields to months.
    if os_time_src and re.search(r"day|days", os_time_src, re.I):
        os_time = os_time / 30.4375
        os_time_src = os_time_src + " (converted days/30.4375)"
    if rfs_time_src and re.search(r"day|days", rfs_time_src, re.I):
        rfs_time = rfs_time / 30.4375
        rfs_time_src = rfs_time_src + " (converted days/30.4375)"
    # Some GEO fields called "survival time" are in years despite lacking units.
    if gse == "GSE50081" and os_time_src == "survival_time":
        os_time = os_time * 12
        os_time_src = os_time_src + " (converted years*12; GEO field lacks explicit unit)"
    if gse == "GSE50081" and rfs_time_src == "disease_free_survival_time":
        rfs_time = rfs_time * 12
        rfs_time_src = rfs_time_src + " (converted years*12; GEO field lacks explicit unit)"

    out["OS_time_months"] = os_time
    out["OS_event"] = os_event
    out["OS_time_source"] = os_time_src
    out["OS_event_source"] = os_event_src
    out["RFS_time_months"] = rfs_time
    out["RFS_event"] = rfs_event
    out["RFS_time_source"] = rfs_time_src
    out["RFS_event_source"] = rfs_event_src

    for label, pats in [
        ("stage", [r"\bstage\b", r"pathologic.*stage", r"tumor_stage"]),
        ("histology", [r"histolog", r"subtype", r"cell_type", r"disease_state"]),
        ("sample_type", [r"tissue", r"source_name", r"disease_state"]),
        ("sex", [r"\bsex\b", r"gender"]),
        ("age", [r"\bage\b"]),
    ]:
        candidates = find_cols(pats)
        if candidates:
            best = max(candidates, key=lambda c: df[c].astype(str).replace("", np.nan).notna().sum())
            out[label] = df[best]
            out[f"{label}_source"] = best
        else:
            out[label] = ""
            out[f"{label}_source"] = ""

    return out


def summary_for(clin: pd.DataFrame, raw: pd.DataFrame) -> dict:
    gse = clin["GSE"].iloc[0] if len(clin) else ""
    os_complete = clin["OS_time_months"].notna() & clin["OS_event"].notna()
    rfs_complete = clin["RFS_time_months"].notna() & clin["RFS_event"].notna()
    return {
        "GSE": gse,
        "sample_count": len(clin),
        "OS_complete_n": int(os_complete.sum()),
        "OS_event_n": int(clin.loc[os_complete, "OS_event"].sum()) if os_complete.any() else 0,
        "RFS_complete_n": int(rfs_complete.sum()),
        "RFS_event_n": int(clin.loc[rfs_complete, "RFS_event"].sum()) if rfs_complete.any() else 0,
        "OS_time_source": clin["OS_time_source"].dropna().astype(str).replace("", np.nan).dropna().iloc[0] if clin["OS_time_source"].replace("", np.nan).notna().any() else "",
        "OS_event_source": clin["OS_event_source"].dropna().astype(str).replace("", np.nan).dropna().iloc[0] if clin["OS_event_source"].replace("", np.nan).notna().any() else "",
        "RFS_time_source": clin["RFS_time_source"].dropna().astype(str).replace("", np.nan).dropna().iloc[0] if clin["RFS_time_source"].replace("", np.nan).notna().any() else "",
        "RFS_event_source": clin["RFS_event_source"].dropna().astype(str).replace("", np.nan).dropna().iloc[0] if clin["RFS_event_source"].replace("", np.nan).notna().any() else "",
        "machine_survival_ready": "yes" if os_complete.sum() >= 30 and clin.loc[os_complete, "OS_event"].sum() >= 5 else "no",
        "raw_metadata_columns": len(raw.columns),
    }


def main():
    all_clin = []
    summaries = []
    raw_debug_dir = WORK / "geo_full_sample_metadata"
    raw_debug_dir.mkdir(exist_ok=True)
    for raw_gz in sorted(RAW_GEO.glob("GSE*/GSE*_series_matrix.txt.gz")):
        gse = raw_gz.parent.name
        try:
            raw = parse_series_sample_metadata(raw_gz)
            raw.to_csv(raw_debug_dir / f"{gse}.sample_metadata_all.tsv", sep="\t", index=False)
            clin = extract_clinical(raw, gse)
            all_clin.append(clin)
            summaries.append(summary_for(clin, raw))
            print(gse, summaries[-1]["machine_survival_ready"], "OS", summaries[-1]["OS_complete_n"], "events", summaries[-1]["OS_event_n"])
        except Exception as exc:
            summaries.append({"GSE": gse, "sample_count": 0, "machine_survival_ready": "error", "error": repr(exc)})
            print(gse, "error", exc)

    sample = pd.concat(all_clin, ignore_index=True) if all_clin else pd.DataFrame()
    sample.to_csv(OUT_SAMPLE, index=False, encoding="utf-8-sig")
    pd.DataFrame(summaries).to_csv(OUT_SUMMARY, index=False, encoding="utf-8-sig")
    print(f"wrote {OUT_SAMPLE}")
    print(f"wrote {OUT_SUMMARY}")


if __name__ == "__main__":
    main()
