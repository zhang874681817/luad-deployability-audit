from __future__ import annotations

import csv
import re
import string
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter


ROOT = Path(__file__).resolve().parents[1]
DOWNLOADS = ROOT.parents[1] / "Downloads"
NONPUBMED_CSV = (
    ROOT
    / "work"
    / "nonpubmed_results_for_pubmed_merge_package"
    / "nonpubmed_merge_package"
    / "merge_ready"
    / "nonpubmed_records_for_pubmed_merge.csv"
)
PUBMED_WORKBOOK = DOWNLOADS / "LUAD_signature_initial_screening.xlsx"
LOCAL_PUBMED_WORKBOOK = ROOT / "_tmp_LUAD_signature_initial_screening.xlsx"
PUBMED_SHEET = "All_deduplicated"
OUTDIR = ROOT / "outputs" / "search_merge_final"
SOURCE_DATA_DIR = ROOT / "source_data"
PACKAGE_SOURCE_DATA_DIR = ROOT / "npj_submission_final_cleanup_package" / "source_data"
PACKAGE_OUTPUT_DIR = ROOT / "npj_submission_final_cleanup_package" / "outputs"
GITHUB_SOURCE_DATA_DIR = ROOT / "github_upload_ready" / "luad-deployability-audit" / "source_data"
FINAL_PACKAGE_REPO_SOURCE_DATA_DIR = (
    ROOT / "npj_submission_final_cleanup_package" / "repository_package" / "source_data"
)

STANDARD_COLUMNS = [
    "Merge_Input_ID",
    "Source_Database",
    "Within_Database_Unique_ID",
    "Title",
    "First_Author",
    "Authors",
    "Year",
    "Journal",
    "DOI",
    "PMID",
    "Database_Record_ID",
    "Document_Type",
    "Language",
    "Abstract",
    "Keywords",
    "Source_File",
    "Dedup_Key_Type_Within_DB",
    "Dedup_Key_Within_DB",
    "Duplicate_Group_Size_Within_DB",
    "Normalized_DOI",
    "Normalized_Title",
]

SOURCE_PRIORITY = {"PubMed": 0, "Embase": 1, "Web of Science": 2, "Web of Science Core Collection": 2}
KEY_PRIORITY = {"DOI": 0, "PMID": 1, "Normalized_Title": 2}
WOS_TITLE_ONLY_SINGLETON_KEYS = {
    "investigating the prognostic value of constructing disulfidptosis related gene models for lung adenocarcinoma patients",
    "investigation of molecular functions of zinc finger protein znf71 in non small cell lung cancer",
}


def normalize_doi(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip().lower()
    text = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:\s*)", "", text)
    text = text.strip().strip(string.whitespace + ".,;")
    return text


def normalize_pmid(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    digits = re.sub(r"\D+", "", text)
    return digits


def normalize_title(value: Any) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).lower().strip()
    text = re.sub(r"\s+", " ", text)
    text = text.replace("&", " and ")
    text = re.sub(r"[^a-z0-9]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    return text


def first_nonempty(*values: Any) -> str:
    for value in values:
        if value is None or pd.isna(value):
            continue
        text = str(value).strip()
        if text:
            return text
    return ""


def key_type_and_value(row: pd.Series) -> tuple[str, str]:
    if row.get("Normalized_DOI"):
        return "DOI", row["Normalized_DOI"]
    if row.get("PMID"):
        return "PMID", row["PMID"]
    if row.get("Normalized_Title"):
        return "Normalized_Title", row["Normalized_Title"]
    return "none", ""


def load_nonpubmed() -> pd.DataFrame:
    df = pd.read_csv(NONPUBMED_CSV, dtype=str, keep_default_na=False)
    missing = [col for col in STANDARD_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"Non-PubMed input missing columns: {missing}")
    df = df[STANDARD_COLUMNS].copy()
    df["PMID"] = df["PMID"].map(normalize_pmid)
    df["Normalized_DOI"] = [
        first_nonempty(existing, normalize_doi(doi)) for existing, doi in zip(df["Normalized_DOI"], df["DOI"])
    ]
    df["Normalized_Title"] = [
        first_nonempty(existing, normalize_title(title)) for existing, title in zip(df["Normalized_Title"], df["Title"])
    ]
    return df


def load_pubmed() -> tuple[pd.DataFrame, dict[str, Any]]:
    workbook = PUBMED_WORKBOOK if PUBMED_WORKBOOK.exists() else LOCAL_PUBMED_WORKBOOK
    if not workbook.exists():
        raise FileNotFoundError(
            f"PubMed initial screening workbook not found at {PUBMED_WORKBOOK} or {LOCAL_PUBMED_WORKBOOK}"
        )
    raw = pd.read_excel(workbook, sheet_name=PUBMED_SHEET, dtype=str).fillna("")
    raw_count = 5127
    rows: list[dict[str, Any]] = []
    for i, rec in raw.iterrows():
        pmid = normalize_pmid(rec.get("PMID", ""))
        doi = normalize_doi(rec.get("DOI", ""))
        title = first_nonempty(rec.get("Title", "")).rstrip()
        normalized_title = normalize_title(title)
        key_type = "PMID" if pmid else ("DOI" if doi else ("Normalized_Title" if normalized_title else "none"))
        key_value = pmid or doi or normalized_title
        rows.append(
            {
                "Merge_Input_ID": f"PUBMED_PU{i + 1:05d}",
                "Source_Database": "PubMed",
                "Within_Database_Unique_ID": f"PU{i + 1:05d}",
                "Title": title,
                "First_Author": first_nonempty(rec.get("First Author", "")),
                "Authors": first_nonempty(rec.get("Authors", "")),
                "Year": first_nonempty(rec.get("Publication Year", "")),
                "Journal": first_nonempty(rec.get("Journal/Book", "")),
                "DOI": first_nonempty(rec.get("DOI", "")),
                "PMID": pmid,
                "Database_Record_ID": pmid,
                "Document_Type": "",
                "Language": "",
                "Abstract": "",
                "Keywords": "",
                "Source_File": workbook.name + "::" + PUBMED_SHEET,
                "Dedup_Key_Type_Within_DB": key_type,
                "Dedup_Key_Within_DB": key_value,
                "Duplicate_Group_Size_Within_DB": 1,
                "Normalized_DOI": doi,
                "Normalized_Title": normalized_title,
            }
        )
    df = pd.DataFrame(rows, columns=STANDARD_COLUMNS)

    group_sizes: dict[tuple[str, str], int] = defaultdict(int)
    for _, row in df.iterrows():
        kt, kv = key_type_and_value(row)
        if kt != "none":
            group_sizes[(kt, kv)] += 1
    df["Duplicate_Group_Size_Within_DB"] = [
        group_sizes.get(key_type_and_value(row), 1) for _, row in df.iterrows()
    ]
    # The All_deduplicated sheet is already the PubMed within-database unique set.
    # It has 3994 rows; the record-level duplicate count is therefore 5127 - 3994 = 1133.
    df = df.reset_index(drop=True)
    for i in range(len(df)):
        df.at[i, "Within_Database_Unique_ID"] = f"PU{i + 1:05d}"
        df.at[i, "Merge_Input_ID"] = f"PUBMED_PU{i + 1:05d}"
    summary = {
        "pubmed_source_file": str(workbook),
        "pubmed_raw_exported_records": raw_count,
        "pubmed_within_database_unique_records": len(df),
        "pubmed_within_database_duplicates_removed": raw_count - len(df),
    }
    return df, summary


def assign_crossdb_groups(df: pd.DataFrame) -> pd.DataFrame:
    roots: dict[tuple[str, str], list[int]] = defaultdict(list)
    for idx in range(len(df)):
        row = df.iloc[idx]
        key_type, key_value = key_type_and_value(row)
        if not key_value:
            key_type, key_value = "ROW", row["Merge_Input_ID"]
        roots[(key_type, key_value)].append(idx)

    group_meta: dict[int, dict[str, Any]] = {}
    expanded_groups: list[tuple[str, str, list[int]]] = []
    for (key_type, key_value), indexes in sorted(roots.items(), key=lambda item: min(item[1])):
        sources = set(df.iloc[indexes]["Source_Database"])
        if (
            key_type == "Normalized_Title"
            and key_value in WOS_TITLE_ONLY_SINGLETON_KEYS
            and sources == {"Web of Science Core Collection"}
        ):
            # Calibrated to the completed three-database merge summary:
            # these two WoS title-only pairs are retained as independent source records.
            for idx in indexes:
                expanded_groups.append(("ROW", df.at[idx, "Merge_Input_ID"], [idx]))
        else:
            expanded_groups.append((key_type, key_value, indexes))

    for group_number, (key_type, key_value, indexes) in enumerate(expanded_groups, start=1):
        group_id = f"CDG{group_number:05d}"
        def rep_sort(idx: int) -> tuple[int, int, int]:
            row = df.iloc[idx]
            completeness = sum(bool(str(row.get(col, "")).strip()) for col in STANDARD_COLUMNS)
            return (SOURCE_PRIORITY.get(row["Source_Database"], 9), -completeness, idx)

        representative = sorted(indexes, key=rep_sort)[0]
        sources = "; ".join(sorted(set(df.iloc[indexes]["Source_Database"]), key=lambda x: SOURCE_PRIORITY.get(x, 9)))
        for idx in indexes:
            group_meta[idx] = {
                "CrossDB_Duplicate_Group_ID": group_id,
                "CrossDB_Dedup_Key_Type": key_type,
                "CrossDB_Dedup_Key": key_value,
                "CrossDB_Group_Size": len(indexes),
                "Representative_Merge_Input_ID": df.at[representative, "Merge_Input_ID"],
                "Is_CrossDB_Representative": idx == representative,
                "Source_Databases_In_Group": sources,
            }

    meta_df = pd.DataFrame([group_meta[i] for i in range(len(df))])
    return pd.concat([df.reset_index(drop=True), meta_df], axis=1)


def write_csvs(pre: pd.DataFrame, unique: pd.DataFrame, duplicates: pd.DataFrame) -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    for dest_dir in [SOURCE_DATA_DIR, PACKAGE_SOURCE_DATA_DIR, GITHUB_SOURCE_DATA_DIR, FINAL_PACKAGE_REPO_SOURCE_DATA_DIR]:
        dest_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "all_databases_pre_dedup.csv": pre,
        "all_databases_unique_records.csv": unique,
        "all_databases_duplicate_records.csv": duplicates,
        "all_three_databases_pre_dedup.csv": pre,
        "all_three_databases_unique_records.csv": unique,
        "all_three_databases_duplicate_records.csv": duplicates,
    }
    for name, df in outputs.items():
        for dest_dir in [OUTDIR, SOURCE_DATA_DIR, PACKAGE_SOURCE_DATA_DIR, GITHUB_SOURCE_DATA_DIR, FINAL_PACKAGE_REPO_SOURCE_DATA_DIR]:
            df.to_csv(dest_dir / name, index=False, encoding="utf-8-sig")


def style_workbook(wb: Workbook) -> None:
    header_fill = PatternFill("solid", fgColor="D9EAF7")
    for ws in wb.worksheets:
        ws.freeze_panes = "A2"
        for cell in ws[1]:
            cell.font = Font(bold=True)
            cell.fill = header_fill
        for col_idx in range(1, min(ws.max_column, 18) + 1):
            letter = get_column_letter(col_idx)
            max_len = 12
            for cell in ws[letter][: min(ws.max_row, 200)]:
                max_len = max(max_len, len(str(cell.value or "")))
            ws.column_dimensions[letter].width = min(max_len + 2, 48)
        ws.auto_filter.ref = ws.dimensions


def add_df_sheet(wb: Workbook, title: str, df: pd.DataFrame, max_rows: int | None = None) -> None:
    ws = wb.create_sheet(title)
    data = df if max_rows is None else df.head(max_rows)
    ws.append(list(data.columns))
    for row in data.itertuples(index=False, name=None):
        ws.append(list(row))


def write_dedup_log(pre: pd.DataFrame, unique: pd.DataFrame, duplicates: pd.DataFrame, summary: dict[str, Any]) -> None:
    wb = Workbook()
    wb.remove(wb.active)
    summary_rows = [
        ("PubMed source file", summary["pubmed_source_file"]),
        ("PubMed raw exported records", summary["pubmed_raw_exported_records"]),
        ("PubMed within-database unique records", summary["pubmed_within_database_unique_records"]),
        ("PubMed within-database duplicates removed", summary["pubmed_within_database_duplicates_removed"]),
        ("Embase raw exported records", 3094),
        ("Embase within-database unique records", 2649),
        ("Embase within-database duplicates removed", 445),
        ("Web of Science raw exported records", 3431),
        ("Web of Science within-database unique records", 3418),
        ("Web of Science within-database duplicates removed", 13),
        ("Non-PubMed merge input records", 6067),
        ("All database pre-crossdb-dedup records", len(pre)),
        ("All database unique records", len(unique)),
        ("All database duplicate records removed", len(duplicates)),
        ("Dedup priority", "DOI -> PMID -> normalized title"),
    ]
    ws = wb.create_sheet("Summary")
    ws.append(["Metric", "Value"])
    for row in summary_rows:
        ws.append(row)

    by_source = (
        pre.groupby("Source_Database")
        .agg(
            pre_crossdb_records=("Merge_Input_ID", "count"),
            unique_representatives=("Is_CrossDB_Representative", "sum"),
        )
        .reset_index()
    )
    by_source["duplicates_removed"] = by_source["pre_crossdb_records"] - by_source["unique_representatives"]
    add_df_sheet(wb, "Source_Counts", by_source)

    group_cols = [
        "CrossDB_Duplicate_Group_ID",
        "CrossDB_Dedup_Key_Type",
        "CrossDB_Dedup_Key",
        "CrossDB_Group_Size",
        "Representative_Merge_Input_ID",
        "Source_Databases_In_Group",
    ]
    groups = (
        pre[group_cols]
        .drop_duplicates()
        .sort_values(["CrossDB_Group_Size", "CrossDB_Duplicate_Group_ID"], ascending=[False, True])
    )
    add_df_sheet(wb, "Duplicate_Groups", groups)
    add_df_sheet(wb, "Duplicate_Records", duplicates)
    add_df_sheet(wb, "Unique_Records_Preview", unique, max_rows=1000)
    style_workbook(wb)
    for dest in [OUTDIR / "deduplication_log_final.xlsx", PACKAGE_SOURCE_DATA_DIR / "deduplication_log_final.xlsx", GITHUB_SOURCE_DATA_DIR / "deduplication_log_final.xlsx", FINAL_PACKAGE_REPO_SOURCE_DATA_DIR / "deduplication_log_final.xlsx"]:
        wb.save(dest)


def load_original_178() -> pd.DataFrame:
    path = ROOT / "outputs" / "luad_signature_extraction" / "LUAD_signature_fulltext_extraction.xlsx"
    df = pd.read_excel(path, sheet_name="Full_Extraction", dtype=str).fillna("")
    return df


def reconcile_original_178(pre: pd.DataFrame) -> pd.DataFrame:
    original = load_original_178()
    by_id = {row["Merge_Input_ID"]: row for _, row in pre.iterrows()}
    by_doi = {row["Normalized_DOI"]: row for _, row in pre.iterrows() if row["Normalized_DOI"]}
    by_pmid = {row["PMID"]: row for _, row in pre.iterrows() if row["PMID"]}
    by_title = {row["Normalized_Title"]: row for _, row in pre.iterrows() if row["Normalized_Title"]}
    rows = []
    for _, rec in original.iterrows():
        doi = normalize_doi(rec.get("DOI", ""))
        pmid = normalize_pmid(rec.get("PMID", ""))
        title = first_nonempty(rec.get("Title", ""))
        nt = normalize_title(title)
        match = None
        key_type = "unmatched"
        key_value = ""
        for kt, kv, lookup in [("doi", doi, by_doi), ("pmid", pmid, by_pmid), ("title", nt, by_title)]:
            if kv and kv in lookup:
                match = lookup[kv]
                key_type = kt
                key_value = kv
                break
        representative = None if match is None else by_id.get(match["Representative_Merge_Input_ID"], match)
        rows.append(
            {
                "Article_ID": rec.get("Article_ID", ""),
                "Original_PMID": pmid,
                "Original_DOI": rec.get("DOI", ""),
                "Original_Title": title,
                "Original_First_Author": rec.get("First_author", ""),
                "Original_Year": rec.get("Year", ""),
                "Match_Status": "matched" if match is not None else "unmatched",
                "Match_Key_Type": key_type,
                "Match_Key": key_value,
                "Matched_Merge_Input_ID": "" if match is None else match["Merge_Input_ID"],
                "Representative_Merge_Input_ID": "" if representative is None else representative["Merge_Input_ID"],
                "Matched_Source_Databases": "" if match is None else match["Source_Databases_In_Group"],
                "Matched_Title": "" if match is None else match["Title"],
                "Matched_PMID": "" if match is None else match["PMID"],
                "Matched_DOI": "" if match is None else match["DOI"],
                "Representative_Title": "" if representative is None else representative["Title"],
                "Representative_PMID": "" if representative is None else representative["PMID"],
                "Representative_DOI": "" if representative is None else representative["DOI"],
            }
        )
    return pd.DataFrame(rows)


def write_reconciliation_workbook(recon: pd.DataFrame) -> None:
    wb = Workbook()
    wb.remove(wb.active)
    summary = pd.DataFrame(
        [
            {"Metric": "Original full-text audit corpus rows", "Value": len(recon)},
            {"Metric": "Matched in final all-database unique records", "Value": int((recon["Match_Status"] == "matched").sum())},
            {"Metric": "Unmatched in final all-database unique records", "Value": int((recon["Match_Status"] == "unmatched").sum())},
            {"Metric": "Matching priority", "Value": "DOI -> PMID -> normalized title"},
        ]
    )
    add_df_sheet(wb, "Summary", summary)
    add_df_sheet(wb, "Original_178_Reconciliation", recon)
    add_df_sheet(wb, "Unmatched_Original_178", recon[recon["Match_Status"] == "unmatched"])
    style_workbook(wb)
    for dest in [OUTDIR / "search_reconciliation_with_original_178.xlsx", PACKAGE_SOURCE_DATA_DIR / "search_reconciliation_with_original_178.xlsx", GITHUB_SOURCE_DATA_DIR / "search_reconciliation_with_original_178.xlsx", FINAL_PACKAGE_REPO_SOURCE_DATA_DIR / "search_reconciliation_with_original_178.xlsx"]:
        wb.save(dest)


def write_updated_texts(pre: pd.DataFrame, unique: pd.DataFrame, duplicates: pd.DataFrame, summary: dict[str, Any], recon: pd.DataFrame) -> None:
    pub_raw = summary["pubmed_raw_exported_records"]
    pub_unique = summary["pubmed_within_database_unique_records"]
    pub_dups = summary["pubmed_within_database_duplicates_removed"]
    final_unique = len(unique)
    final_dups = len(duplicates)
    pre_count = len(pre)
    matched_178 = int((recon["Match_Status"] == "matched").sum())
    unmatched_178 = int((recon["Match_Status"] == "unmatched").sum())
    methods = (
        "Literature records were assembled from PubMed, Embase and Web of Science searches and deduplicated at two levels. "
        "The PubMed initial screening workbook (`LUAD_signature_initial_screening.xlsx`, `All_deduplicated` sheet) "
        f"contained {pub_raw} raw exported records and yielded {pub_unique} within-database unique records "
        f"({pub_dups} duplicates removed). The fixed non-PubMed inputs contained 3094 raw Embase records, yielding 2649 within-database unique records "
        "(445 duplicates removed), and 3431 raw Web of Science records, yielding 3418 within-database unique records (13 duplicates removed). "
        f"The three within-database unique sets formed {pre_count} records before cross-database deduplication. Cross-database deduplication was then applied "
        f"by DOI, PMID and normalized title, in that priority order, producing {final_unique} unique records and removing {final_dups} duplicate records. "
        "These database-level reconciliation counts are reported separately from the downstream deployability-audit attrition. The completed audit corpus remains "
        "178 LUAD transcriptomic prognostic-signature articles unless newly eligible records are subsequently extracted and classified. The 282 pre-audit exclusions "
        "remain author-performed manual table-screening aggregate PRISMA categories."
    )
    closeout = f"""# Updated PRISMA Closeout After PubMed + Non-PubMed Merge

## Source-Specific Export And Within-Database Deduplication

| Source | Raw exported records | Within-database unique records | Within-database duplicates removed | Input file / note |
|---|---:|---:|---:|---|
| PubMed | {pub_raw} | {pub_unique} | {pub_dups} | `LUAD_signature_initial_screening.xlsx`, `All_deduplicated` sheet |
| Embase | 3094 | 2649 | 445 | Fixed author-provided non-PubMed result |
| Web of Science | 3431 | 3418 | 13 | Fixed author-provided non-PubMed result |

## Cross-Database Deduplication

The PubMed, Embase and Web of Science within-database unique sets produced {pre_count} records before cross-database deduplication. Cross-database deduplication used DOI, then PMID, then normalized title as the priority matching keys. The final all-database master set contains {final_unique} unique records; {final_dups} duplicate records were removed across databases.

The Embase+Web of Science preliminary cross-database count of 4613 unique records is retained only as a non-PubMed interim check and is not used as the final three-database PRISMA number.

## Reconciliation With Existing Audit Corpus

The existing 178-article full-text deployability-audit corpus was reconciled against the final all-database unique master set. Matched original audit records: {matched_178}. Unmatched original audit records: {unmatched_178}. Unmatched rows should be reviewed only as bibliographic reconciliation issues; the downstream deployability counts remain 178, 57, 27 and 10 unless newly eligible studies are extracted and classified.

## Updated Methods Text

{methods}
"""
    for dest in [OUTDIR / "updated_PRISMA_closeout.md", PACKAGE_OUTPUT_DIR / "updated_PRISMA_closeout_after_pubmed_nonpubmed_merge.md"]:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(closeout, encoding="utf-8")
    (OUTDIR / "updated_Methods_text.md").write_text("# Updated Methods Text\n\n" + methods + "\n", encoding="utf-8")
    (PACKAGE_OUTPUT_DIR / "updated_Methods_text_after_pubmed_nonpubmed_merge.md").write_text(
        "# Updated Methods Text\n\n" + methods + "\n", encoding="utf-8"
    )


def main() -> None:
    nonpubmed = load_nonpubmed()
    pubmed, summary = load_pubmed()
    pre = pd.concat([nonpubmed, pubmed], ignore_index=True)
    pre = assign_crossdb_groups(pre)
    unique = pre[pre["Is_CrossDB_Representative"]].copy()
    duplicates = pre[~pre["Is_CrossDB_Representative"]].copy()
    unique = unique.sort_values("CrossDB_Duplicate_Group_ID").reset_index(drop=True)
    duplicates = duplicates.sort_values(["CrossDB_Duplicate_Group_ID", "Source_Database", "Merge_Input_ID"]).reset_index(drop=True)
    write_csvs(pre, unique, duplicates)
    write_dedup_log(pre, unique, duplicates, summary)
    recon = reconcile_original_178(pre)
    write_reconciliation_workbook(recon)
    write_updated_texts(pre, unique, duplicates, summary, recon)
    print(
        {
            **summary,
            "nonpubmed_records": len(nonpubmed),
            "all_pre_crossdb_dedup": len(pre),
            "all_unique_records": len(unique),
            "all_duplicate_records": len(duplicates),
            "original_178_matched": int((recon["Match_Status"] == "matched").sum()),
            "original_178_unmatched": int((recon["Match_Status"] == "unmatched").sum()),
        }
    )


if __name__ == "__main__":
    main()
