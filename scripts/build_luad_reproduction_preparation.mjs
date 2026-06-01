import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = process.cwd();
const inputJson = path.join(root, "work", "luad_reproduction_preparation.json");
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_reproduction_preparation_workbook.xlsx");

function sanitize(value) {
  if (value === null || value === undefined) return "";
  if (typeof value === "number") return Number.isFinite(value) ? value : "";
  const s = String(value);
  return s.length > 32000 ? s.slice(0, 31950) + " [truncated]" : s;
}

function colLetter(n) {
  let s = "";
  while (n > 0) {
    const m = (n - 1) % 26;
    s = String.fromCharCode(65 + m) + s;
    n = Math.floor((n - 1) / 26);
  }
  return s;
}

function widthFor(field) {
  if (/Title|Formula|Context|Notes|Recommended|Candidate|Excluded|Action|Hint|Unmatched/.test(field)) return 48;
  if (/Gene|Cohort|Validation|Training|Coefficient/.test(field)) return 34;
  if (/Count|Value|PMID|Year|Order|Priority|Ready/.test(field)) return 14;
  return 20;
}

async function addSheet(wb, name, rows, fields) {
  const sheet = wb.worksheets.add(name);
  const data = [fields, ...rows.map((r) => fields.map((f) => sanitize(r[f])))];
  const endCol = colLetter(fields.length);
  const endRow = data.length;
  sheet.getRange(`A1:${endCol}${endRow}`).values = data;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${endCol}1`).format = {
    fill: "#375623",
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
  };
  sheet.getRange(`A1:${endCol}${endRow}`).format = { wrapText: true };
  fields.forEach((field, idx) => {
    sheet.getRange(`${colLetter(idx + 1)}:${colLetter(idx + 1)}`).format.columnWidth = widthFor(field);
  });
  try {
    sheet.tables.add(`A1:${endCol}${endRow}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {
    // Keep plain range if a table cannot be added.
  }
}

await fs.mkdir(outputDir, { recursive: true });
const payload = JSON.parse(await fs.readFile(inputJson, "utf8"));

const wb = Workbook.create();
await addSheet(wb, "Summary", [
  ...payload.summary,
  { Metric: "HGNC_source_url", Value: payload.hgnc_source_url },
  { Metric: "HGNC_local_file", Value: payload.hgnc_local_file },
], ["Metric", "Value"]);

await addSheet(wb, "Signature_Readiness", payload.signature_readiness, [
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Title",
  "Omics_type",
  "Endpoint",
  "Gene_count",
  "Numeric_coefficient_count",
  "HGNC_mapped_count",
  "HGNC_approved_symbol_count",
  "Unmatched_gene_count",
  "Non_HGNC_expression_feature_count",
  "Unmatched_genes",
  "Alias_or_previous_symbol_genes",
  "Ready_for_expression_scoring",
  "Training_cohort",
  "Validation_cohort",
  "Usable_article_validation_cohorts",
  "Training_validation_overlap",
  "Candidate_external_pool_minus_training",
  "Cutoff_type",
  "Cutoff_value",
]);

await addSheet(wb, "Gene_Coefficients_Long", payload.signature_gene_coefficients, [
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Title",
  "Signature_name",
  "Omics_type",
  "Endpoint",
  "Gene_order",
  "Gene_original",
  "Coefficient_raw",
  "Coefficient_numeric",
  "Coefficient_direction",
  "HGNC_match_status",
  "Approved_symbol",
  "HGNC_ID",
  "Locus_group",
  "Locus_type",
  "Match_source",
  "Training_cohort",
  "Validation_cohort",
  "Cutoff_type",
  "Cutoff_value",
  "Risk_score_formula",
]);

await addSheet(wb, "External_Validation_Plan", payload.external_validation_plan, [
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Title",
  "Training_cohorts_parsed",
  "Article_validation_cohorts_parsed",
  "Excluded_from_external_validation",
  "Usable_article_validation_cohorts",
  "All_candidate_external_cohorts_minus_training",
  "Recommended_first_pass_external_cohorts",
  "Overlap_note",
]);

await addSheet(wb, "Cohort_Download_Manifest", payload.cohort_download_manifest, [
  "Cohort",
  "Cohort_type",
  "Public_download_priority",
  "Training_count",
  "Article_validation_count",
  "Eligible_external_count",
  "Download_hint",
]);

await wb.inspect({
  kind: "table",
  range: "Summary!A1:B8",
  tableMaxRows: 10,
  tableMaxCols: 2,
  maxChars: 1200,
});
try {
  await wb.render({ sheetName: "Signature_Readiness", range: "A1:O15", scale: 1, format: "png" });
} catch {
  // Best-effort visual check.
}
const output = await SpreadsheetFile.exportXlsx(wb);
await output.save(outputPath);
console.log(outputPath);
