import csv
import json
import re
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
ARTICLE_JSON = WORK / "article_extraction_qc.json"
HGNC_TSV = WORK / "hgnc_complete_set.txt"
OUT_JSON = WORK / "luad_reproduction_preparation.json"

HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"


def norm(s):
    return re.sub(r"\s+", " ", str(s or "")).strip()


def clean_gene(g):
    g = norm(g).strip(" ;,")
    g = g.replace("EXP-", "")
    g = re.sub(r"\s+", "", g)
    return g


def signature_id(row):
    pmid = str(row.get("PMID", "") or "").replace(".0", "")
    return f"{row.get('Article_ID')}{('_PMID' + pmid) if pmid else ''}"


def split_semi(v):
    return [x.strip() for x in str(v or "").split(";") if x.strip() and x.strip().lower() != "nan"]


def parse_coefficients(row):
    genes = split_semi(row.get("Gene_list", ""))
    coef_pairs = []
    for item in split_semi(row.get("Coefficients", "")):
        if ":" not in item:
            continue
        gene, coef = item.split(":", 1)
        gene = clean_gene(gene)
        coef_raw = coef.strip()
        coef_num = None
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", coef_raw)
        if m:
            try:
                coef_num = float(m.group(0))
            except Exception:
                coef_num = None
        coef_pairs.append((gene, coef_raw, coef_num))
    if coef_pairs:
        # Preserve coefficient order as authoritative.
        return coef_pairs
    return [(clean_gene(g), "", None) for g in genes]


def download_hgnc():
    if HGNC_TSV.exists() and HGNC_TSV.stat().st_size > 1_000_000:
        return
    WORK.mkdir(exist_ok=True)
    with urllib.request.urlopen(HGNC_URL, timeout=120) as resp:
        HGNC_TSV.write_bytes(resp.read())


def load_hgnc():
    download_hgnc()
    df = pd.read_csv(HGNC_TSV, sep="\t", dtype=str).fillna("")
    by_symbol = {}
    by_alias = {}
    by_prev = {}
    for _, r in df.iterrows():
        symbol = r.get("symbol", "")
        if not symbol:
            continue
        rec = {
            "approved_symbol": symbol,
            "hgnc_id": r.get("hgnc_id", ""),
            "name": r.get("name", ""),
            "locus_group": r.get("locus_group", ""),
            "locus_type": r.get("locus_type", ""),
            "status": r.get("status", ""),
            "prev_symbol": r.get("prev_symbol", ""),
            "alias_symbol": r.get("alias_symbol", ""),
        }
        by_symbol[symbol.upper()] = rec
        for field, target in [("alias_symbol", by_alias), ("prev_symbol", by_prev)]:
            for token in re.split(r"\|", r.get(field, "") or ""):
                token = token.strip()
                if token:
                    target.setdefault(token.upper(), rec)
    return by_symbol, by_alias, by_prev


def map_gene(gene, by_symbol, by_alias, by_prev):
    q = clean_gene(gene)
    key = q.upper()
    if key in by_symbol:
        rec = by_symbol[key]
        return {**rec, "query": q, "match_status": "approved_symbol", "match_source": "symbol"}
    if key in by_prev:
        rec = by_prev[key]
        return {**rec, "query": q, "match_status": "previous_symbol", "match_source": "prev_symbol"}
    if key in by_alias:
        rec = by_alias[key]
        return {**rec, "query": q, "match_status": "alias_symbol", "match_source": "alias_symbol"}
    # Common PDF/expression suffix cleanups.
    stripped = re.sub(r"(\.|\-)$", "", q).upper()
    if stripped in by_symbol:
        rec = by_symbol[stripped]
        return {**rec, "query": q, "match_status": "approved_symbol_after_cleanup", "match_source": "cleanup"}
    return {
        "query": q,
        "approved_symbol": "",
        "hgnc_id": "",
        "name": "",
        "locus_group": "",
        "locus_type": "",
        "status": "",
        "prev_symbol": "",
        "alias_symbol": "",
        "match_status": "no_hgnc_match",
        "match_source": "",
    }


def is_non_hgnc_but_expression_feature(gene):
    g = clean_gene(gene).upper()
    return bool(
        re.fullmatch(r"(AC|AL|AP)\d+\.\d+", g)
        or re.fullmatch(r"RP\d{1,2}-[A-Z0-9]+(?:\.\d+)?", g)
        or re.fullmatch(r"CTD-\d+[A-Z]\d+(?:\.\d+)?", g)
        or re.fullmatch(r"[A-Z0-9-]+-AS\d*", g)
    )


def split_cohorts(v):
    out = []
    for token in re.split(r"[;/|,]", str(v or "")):
        t = norm(token).upper()
        if not t or t == "NAN":
            continue
        t = t.replace("TCGA LUAD", "TCGA-LUAD")
        t = t.replace("TCGA_LUAD", "TCGA-LUAD")
        t = re.sub(r"\s+", " ", t)
        out.append(t)
    return list(dict.fromkeys(out))


def public_cohort_type(cohort):
    if cohort == "TCGA-LUAD":
        return "TCGA RNA-seq"
    if re.fullmatch(r"GSE\d{4,8}", cohort):
        return "GEO"
    if cohort == "ICGC":
        return "ICGC"
    if cohort == "CPTAC":
        return "CPTAC"
    if "INTERNAL" in cohort:
        return "internal/not public"
    return "other/unclear"


def is_public_external_candidate(cohort):
    return public_cohort_type(cohort) in {"TCGA RNA-seq", "GEO", "ICGC", "CPTAC"}


def main():
    rows = json.loads(ARTICLE_JSON.read_text(encoding="utf-8"))
    included = [r for r in rows if r.get("Screening_label") == "Included_reproducible"]
    by_symbol, by_alias, by_prev = load_hgnc()

    long_rows = []
    for row in included:
        sig = signature_id(row)
        for idx, (gene, coef_raw, coef_num) in enumerate(parse_coefficients(row), 1):
            mapped = map_gene(gene, by_symbol, by_alias, by_prev)
            long_rows.append(
                {
                    "Signature_ID": sig,
                    "Article_ID": row.get("Article_ID", ""),
                    "PMID": str(row.get("PMID", "") or "").replace(".0", ""),
                    "Title": row.get("Title", ""),
                    "Signature_name": row.get("Signature_name", ""),
                    "Omics_type": row.get("Omics_type", ""),
                    "Endpoint": row.get("Endpoint", ""),
                    "Gene_order": idx,
                    "Gene_original": gene,
                    "Coefficient_raw": coef_raw,
                    "Coefficient_numeric": coef_num,
                    "Coefficient_direction": "risk" if coef_num is not None and coef_num > 0 else ("protective" if coef_num is not None and coef_num < 0 else "zero_or_missing"),
                    "HGNC_match_status": mapped["match_status"],
                    "Approved_symbol": mapped["approved_symbol"],
                    "HGNC_ID": mapped["hgnc_id"],
                    "Locus_group": mapped["locus_group"],
                    "Locus_type": mapped["locus_type"],
                    "Match_source": mapped["match_source"],
                    "Training_cohort": row.get("Training_cohort", ""),
                    "Validation_cohort": row.get("Validation_cohort", ""),
                    "Cutoff_type": row.get("Cutoff_type", ""),
                    "Cutoff_value": row.get("Cutoff_value", ""),
                    "Risk_score_formula": row.get("Risk_score_formula", ""),
                }
            )

    all_public_pool = sorted(
        {
            c
            for row in included
            for c in split_cohorts(row.get("Training_cohort", "")) + split_cohorts(row.get("Validation_cohort", ""))
            if is_public_external_candidate(c)
        }
    )

    sig_rows = []
    validation_plan = []
    for row in included:
        sig = signature_id(row)
        genes = [r for r in long_rows if r["Signature_ID"] == sig]
        total = len(genes)
        numeric = sum(1 for g in genes if g["Coefficient_numeric"] is not None)
        approved = sum(1 for g in genes if g["HGNC_match_status"].startswith("approved"))
        mapped = sum(1 for g in genes if g["Approved_symbol"])
        unmatched = [g["Gene_original"] for g in genes if not g["Approved_symbol"]]
        previous_alias = [g["Gene_original"] for g in genes if g["HGNC_match_status"] in {"previous_symbol", "alias_symbol"}]
        train = split_cohorts(row.get("Training_cohort", ""))
        valid = split_cohorts(row.get("Validation_cohort", ""))
        train_set = set(train)
        usable_article_validation = [c for c in valid if c not in train_set and is_public_external_candidate(c)]
        overlapping_validation = [c for c in valid if c in train_set]
        candidate_pool = [c for c in all_public_pool if c not in train_set]
        acceptable_non_hgnc = sum(1 for g in genes if not g["Approved_symbol"] and is_non_hgnc_but_expression_feature(g["Gene_original"]))
        needs = "yes" if total and numeric == total and (mapped + acceptable_non_hgnc) == total else "check"
        sig_rows.append(
            {
                "Signature_ID": sig,
                "Article_ID": row.get("Article_ID", ""),
                "PMID": str(row.get("PMID", "") or "").replace(".0", ""),
                "Title": row.get("Title", ""),
                "Omics_type": row.get("Omics_type", ""),
                "Endpoint": row.get("Endpoint", ""),
                "Gene_count": total,
                "Numeric_coefficient_count": numeric,
                "HGNC_mapped_count": mapped,
                "HGNC_approved_symbol_count": approved,
                "Unmatched_gene_count": len(unmatched),
                "Unmatched_genes": "; ".join(unmatched),
                "Non_HGNC_expression_feature_count": acceptable_non_hgnc,
                "Alias_or_previous_symbol_genes": "; ".join(previous_alias),
                "Ready_for_expression_scoring": needs,
                "Training_cohort": row.get("Training_cohort", ""),
                "Validation_cohort": row.get("Validation_cohort", ""),
                "Usable_article_validation_cohorts": "; ".join(usable_article_validation),
                "Training_validation_overlap": "; ".join(overlapping_validation),
                "Candidate_external_pool_minus_training": "; ".join(candidate_pool),
                "Cutoff_type": row.get("Cutoff_type", ""),
                "Cutoff_value": row.get("Cutoff_value", ""),
            }
        )
        validation_plan.append(
            {
                "Signature_ID": sig,
                "Article_ID": row.get("Article_ID", ""),
                "PMID": str(row.get("PMID", "") or "").replace(".0", ""),
                "Title": row.get("Title", ""),
                "Training_cohorts_parsed": "; ".join(train),
                "Article_validation_cohorts_parsed": "; ".join(valid),
                "Excluded_from_external_validation": "; ".join(train),
                "Usable_article_validation_cohorts": "; ".join(usable_article_validation),
                "All_candidate_external_cohorts_minus_training": "; ".join(candidate_pool),
                "Recommended_first_pass_external_cohorts": "; ".join(usable_article_validation or candidate_pool[:5]),
                "Overlap_note": "Article validation includes training-overlap cohorts: " + "; ".join(overlapping_validation) if overlapping_validation else "",
            }
        )

    cohort_counts = defaultdict(lambda: Counter())
    for row in included:
        train = split_cohorts(row.get("Training_cohort", ""))
        valid = split_cohorts(row.get("Validation_cohort", ""))
        for c in train:
            cohort_counts[c]["training_count"] += 1
        for c in valid:
            cohort_counts[c]["article_validation_count"] += 1
        for c in set(all_public_pool) - set(train):
            cohort_counts[c]["eligible_external_count"] += 1
    cohort_rows = []
    for c in sorted(cohort_counts):
        ctype = public_cohort_type(c)
        cohort_rows.append(
            {
                "Cohort": c,
                "Cohort_type": ctype,
                "Public_download_priority": "yes" if ctype in {"TCGA RNA-seq", "GEO", "ICGC", "CPTAC"} else "no",
                "Training_count": cohort_counts[c]["training_count"],
                "Article_validation_count": cohort_counts[c]["article_validation_count"],
                "Eligible_external_count": cohort_counts[c]["eligible_external_count"],
                "Download_hint": "TCGA GDC/Xena; harmonize Ensembl-to-symbol and log2(TPM+1) carefully"
                if c == "TCGA-LUAD"
                else ("GEOquery/GEOmetadb; inspect platform GPL and expression scale" if ctype == "GEO" else ""),
            }
        )

    summary = [
        {"Metric": "Included_reproducible_signatures", "Value": len(included)},
        {"Metric": "Signature_gene_rows", "Value": len(long_rows)},
        {"Metric": "Unique_gene_original", "Value": len({r["Gene_original"] for r in long_rows})},
        {"Metric": "Unique_approved_symbols", "Value": len({r["Approved_symbol"] for r in long_rows if r["Approved_symbol"]})},
        {"Metric": "Rows_without_HGNC_match", "Value": sum(1 for r in long_rows if not r["Approved_symbol"])},
        {"Metric": "Public_candidate_cohorts", "Value": len(all_public_pool)},
    ]

    payload = {
        "summary": summary,
        "signature_gene_coefficients": long_rows,
        "signature_readiness": sig_rows,
        "external_validation_plan": validation_plan,
        "cohort_download_manifest": cohort_rows,
        "hgnc_source_url": HGNC_URL,
        "hgnc_local_file": str(HGNC_TSV),
    }
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    for key, records in payload.items():
        if isinstance(records, list):
            out = WORK / f"luad_{key}.csv"
            if records:
                with out.open("w", encoding="utf-8-sig", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=list(records[0].keys()))
                    writer.writeheader()
                    writer.writerows(records)
    print(f"wrote {OUT_JSON}")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
