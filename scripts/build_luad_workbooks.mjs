import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = process.cwd();
const qcJson = path.join(root, "work", "article_extraction_qc.json");
const inputJson = (await fs.stat(qcJson).then(() => true).catch(() => false))
  ? qcJson
  : path.join(root, "work", "article_extraction.json");
const outputDir = path.join(root, "outputs", "luad_signature_extraction");

const mainFields = [
  "PMID",
  "Title",
  "First_author",
  "Year",
  "Journal",
  "Cancer_type",
  "LUAD_specific",
  "Omics_type",
  "Endpoint",
  "Model_type",
  "Training_cohort",
  "Validation_cohort",
  "Public_data_used",
  "GEO_accession",
  "Signature_name",
  "Gene_number",
  "Gene_list",
  "Coefficients_available",
  "Coefficients",
  "Formula_available",
  "Risk_score_formula",
  "Cutoff_available",
  "Cutoff_value",
  "Cutoff_type",
  "Risk_group_definition",
  "Reported_HR",
  "Reported_C_index",
  "Reported_AUC",
  "Independent_prognostic",
  "Nomogram_included",
  "External_validation",
  "Supplement_needed",
  "Reproducible",
  "Exclusion_reason",
  "Notes",
  "Screening_level",
  "Screening_label",
  "Article_ID",
  "DOI",
  "Source_files",
  "Manual_check_needed",
  "Evidence_formula_context",
];

const formulaFields = [
  "Signature_ID",
  "PMID",
  "Genes",
  "Coefficients",
  "Formula",
  "Cutoff",
  "Cutoff_type",
  "Training_cohort",
  "Validation_cohort",
  "Endpoint",
  "Screening_level",
  "Screening_label",
  "Title",
];

function sanitize(value) {
  if (value === null || value === undefined) return "";
  if (typeof value === "number" && !Number.isFinite(value)) return "";
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

function tableRows(rows, fields) {
  return [fields, ...rows.map((r) => fields.map((f) => sanitize(r[f])))];
}

function widthFor(field) {
  if (["Title", "Risk_score_formula", "Evidence_formula_context", "Risk_group_definition", "Exclusion_reason", "Notes", "Source_files"].includes(field)) return 48;
  if (["Gene_list", "Coefficients", "Training_cohort", "Validation_cohort", "GEO_accession"].includes(field)) return 40;
  if (["PMID", "Year", "Gene_number", "Reproducible"].includes(field)) return 13;
  return 20;
}

async function addDataSheet(wb, name, rows, fields) {
  const sheet = wb.worksheets.add(name);
  const data = tableRows(rows, fields);
  const endCol = colLetter(fields.length);
  const endRow = data.length;
  sheet.getRange(`A1:${endCol}${endRow}`).values = data;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${endCol}1`).format = {
    fill: "#1F4E79",
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
    // Table styling is optional; raw ranges remain fully usable if the table add fails.
  }
  return sheet;
}

function signatureId(row) {
  const pmid = sanitize(row.PMID);
  return `${row.Article_ID || "SIG"}${pmid ? `_PMID${pmid.replace(/\.0$/, "")}` : ""}`;
}

function splitCohorts(value) {
  return String(value || "")
    .split(/[;/|,]/)
    .map((x) => x.trim())
    .filter(Boolean)
    .map((x) => x.toUpperCase().replace(/\s+/g, " "))
    .filter((x) => x && x !== "NAN");
}

function candidateRows(rows) {
  return rows.filter((r) => ["Included_reproducible", "Pending_supplement"].includes(r.Screening_label));
}

async function saveWorkbook(fileName, sheets) {
  const wb = Workbook.create();
  for (const spec of sheets) {
    await addDataSheet(wb, spec.name, spec.rows, spec.fields);
  }
  const previewSheet = sheets[0].name;
  await wb.inspect({
    kind: "table",
    range: `${previewSheet}!A1:${colLetter(Math.min(sheets[0].fields.length, 8))}${Math.min(sheets[0].rows.length + 1, 8)}`,
    tableMaxRows: 8,
    tableMaxCols: 8,
    maxChars: 2000,
  });
  try {
    await wb.render({ sheetName: previewSheet, range: `A1:${colLetter(Math.min(sheets[0].fields.length, 8))}${Math.min(sheets[0].rows.length + 1, 12)}`, scale: 1, format: "png" });
  } catch {
    // Rendering can be unavailable in some minimal runtimes; export still validates workbook structure.
  }
  const output = await SpreadsheetFile.exportXlsx(wb);
  const outPath = path.join(outputDir, fileName);
  await output.save(outPath);
  return outPath;
}

await fs.mkdir(outputDir, { recursive: true });
const rows = JSON.parse(await fs.readFile(inputJson, "utf8"));

const reproducible = candidateRows(rows);
const excluded = rows.filter((r) => ["Relevant_not_reproducible", "Exclude"].includes(r.Screening_label));
const formulaRows = rows
  .filter((r) => ["Included_reproducible", "Pending_supplement", "Partially_reproducible"].includes(r.Screening_label))
  .map((r) => ({
    Signature_ID: signatureId(r),
    PMID: r.PMID,
    Genes: r.Gene_list,
    Coefficients: r.Coefficients,
    Formula: r.Risk_score_formula,
    Cutoff: r.Cutoff_value || r.Risk_group_definition,
    Cutoff_type: r.Cutoff_type,
    Training_cohort: r.Training_cohort,
    Validation_cohort: r.Validation_cohort,
    Endpoint: r.Endpoint,
    Screening_level: r.Screening_level,
    Screening_label: r.Screening_label,
    Title: r.Title,
  }));

const cohortCandidates = reproducible.map((r) => {
  const training = splitCohorts(r.Training_cohort);
  const validation = splitCohorts(r.Validation_cohort);
  const all = [...new Set([...training, ...validation])];
  return {
    Signature_ID: signatureId(r),
    PMID: r.PMID,
    Title: r.Title,
    Training_cohort: r.Training_cohort,
    Validation_cohort: r.Validation_cohort,
    Training_cohorts_parsed: training.join("; "),
    Validation_cohorts_parsed: validation.join("; "),
    All_cohorts_parsed: all.join("; "),
    Screening_label: r.Screening_label,
  };
});

const allCohorts = [...new Set(cohortCandidates.flatMap((r) => splitCohorts(r.All_cohorts_parsed)))].sort();
const matrixFields = ["Signature_ID", "PMID", "Title", ...allCohorts];
const matrixRows = cohortCandidates.map((r) => {
  const present = new Set(splitCohorts(r.All_cohorts_parsed));
  const out = { Signature_ID: r.Signature_ID, PMID: r.PMID, Title: r.Title };
  for (const cohort of allCohorts) out[cohort] = present.has(cohort) ? 1 : 0;
  return out;
});

const outputs = [];
outputs.push(await saveWorkbook("LUAD_signature_fulltext_extraction.xlsx", [{ name: "Full_Extraction", rows, fields: mainFields }]));
outputs.push(await saveWorkbook("LUAD_reproducible_signature_inventory.xlsx", [{ name: "Reproducible", rows: reproducible, fields: mainFields }]));
outputs.push(await saveWorkbook("LUAD_excluded_articles_with_reasons.xlsx", [{ name: "Excluded_Level4_5", rows: excluded, fields: mainFields }]));
outputs.push(await saveWorkbook("LUAD_signature_formula_table.xlsx", [{ name: "Formula_Table", rows: formulaRows, fields: formulaFields }]));
outputs.push(
  await saveWorkbook("LUAD_signature_cohort_overlap_matrix.xlsx", [
    {
      name: "Candidate_Cohorts",
      rows: cohortCandidates,
      fields: [
        "Signature_ID",
        "PMID",
        "Title",
        "Training_cohort",
        "Validation_cohort",
        "Training_cohorts_parsed",
        "Validation_cohorts_parsed",
        "All_cohorts_parsed",
        "Screening_label",
      ],
    },
    { name: "Overlap_Matrix", rows: matrixRows, fields: matrixFields },
  ]),
);

console.log(JSON.stringify(outputs, null, 2));
