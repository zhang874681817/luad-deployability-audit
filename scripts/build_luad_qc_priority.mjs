import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = process.cwd();
const qcJson = path.join(root, "work", "article_extraction_qc.json");
const inputJson = (await fs.stat(qcJson).then(() => true).catch(() => false))
  ? qcJson
  : path.join(root, "work", "article_extraction.json");
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_reproducible_QC_priority.xlsx");

const artifactTerms = new Set([
  "AND",
  "AT",
  "BASED",
  "BLOOD",
  "CUTOFF",
  "EXPRESSION",
  "FIG",
  "FIG.",
  "FOR",
  "GENDER",
  "HAZARD",
  "OF",
  "OVERALL",
  "PERFORMANCEEVALUATIONAND",
  "RESSION",
  "RESPECTIVELY",
  "RISK",
  "SCORE",
  "SENSITIVITY",
  "STAGE",
  "TABLE",
  "TABLE1",
  "TIME",
  "TO",
]);

function sanitize(value) {
  if (value === null || value === undefined) return "";
  const s = String(value);
  return s.length > 32000 ? s.slice(0, 31950) + " [truncated]" : s;
}

function splitSemi(value) {
  return String(value || "")
    .split(";")
    .map((x) => x.trim())
    .filter(Boolean);
}

function splitCohorts(value) {
  return String(value || "")
    .split(/[;/|,]/)
    .map((x) => x.trim().toUpperCase())
    .filter(Boolean)
    .map((x) => x.replace(/\s+/g, " "));
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

function extractCoefGenes(coefficients) {
  return splitSemi(coefficients)
    .map((x) => x.split(":")[0]?.trim())
    .filter(Boolean);
}

function suspiciousGenes(geneList) {
  return splitSemi(geneList).filter((gene) => {
    const g = gene.toUpperCase().trim();
    if (artifactTerms.has(g)) return true;
    if (/^(FIG|TABLE|SUPPLE|PERFORMANCE|RESSION|SION|EXPRESSION|HAZARD|GENDER|TIME|SENSITIVITY)/.test(g)) return true;
    if (/[-.]$/.test(g)) return true;
    if (/^[A-Z]$/.test(g)) return true;
    if (/^\d/.test(g)) return true;
    return false;
  });
}

function invalidGseTokens(...values) {
  const joined = values.join("; ");
  const tokens = joined.match(/\bGSE\d+\b/gi) || [];
  return [...new Set(tokens.map((x) => x.toUpperCase()).filter((x) => !/^GSE\d{4,7}$/.test(x)))];
}

function cohortOverlap(training, validation) {
  const train = new Set(splitCohorts(training));
  const valid = splitCohorts(validation);
  return [...new Set(valid.filter((x) => train.has(x)))];
}

function titleLooksSuspicious(title) {
  const t = String(title || "").trim();
  return !t || t.length < 35 || /^(background|abstract|lung cancer remains|introduction)\b/i.test(t);
}

function qcRow(row) {
  const flags = [];
  const actions = [];
  const genes = splitSemi(row.Gene_list);
  const coefGenes = extractCoefGenes(row.Coefficients);
  const badGenes = suspiciousGenes(row.Gene_list);
  const overlap = cohortOverlap(row.Training_cohort, row.Validation_cohort);
  const invalidGse = invalidGseTokens(row.Training_cohort, row.Validation_cohort, row.GEO_accession);

  if (!row.PMID) {
    flags.push("missing_pmid");
    actions.push("Confirm PMID from PDF title/DOI.");
  }
  if (titleLooksSuspicious(row.Title)) {
    flags.push("suspicious_title");
    actions.push("Replace extracted title with article title from PDF/PubMed.");
  }
  if (row.Screening_label === "Pending_supplement") {
    flags.push("pending_supplement");
    actions.push("Open supplementary table and extract missing coefficients/cutoff.");
  }
  if (row.Formula_available !== "yes") {
    flags.push("formula_missing");
    actions.push("Locate formula in main text, figure, or supplement.");
  }
  if (row.Coefficients_available !== "yes") {
    flags.push("coefficients_missing");
    actions.push("Extract beta/coefficient values before risk-score reproduction.");
  }
  if (!row.Cutoff_type || row.Cutoff_type === "unclear") {
    flags.push("cutoff_unclear");
    actions.push("Confirm whether cutoff is median, optimal, training-derived, Youden, X-tile, or surv_cutpoint.");
  }
  if (!row.Cutoff_value && !/median|training-derived|optimal|youden|x-tile|surv_cutpoint/i.test(String(row.Cutoff_type || ""))) {
    flags.push("cutoff_value_or_rule_missing");
    actions.push("Record numeric cutoff or grouping rule.");
  }
  if (badGenes.length) {
    flags.push(`suspicious_gene_tokens:${badGenes.join("|")}`);
    actions.push("Manually correct gene list and coefficient parsing from the formula.");
  }
  if (genes.length !== coefGenes.length && row.Coefficients_available === "yes") {
    flags.push(`gene_coef_count_mismatch:${genes.length}_genes_vs_${coefGenes.length}_coefficients`);
    actions.push("Reconcile gene list with coefficient list.");
  }
  if (invalidGse.length) {
    flags.push(`invalid_gse_token:${invalidGse.join("|")}`);
    actions.push("Correct malformed GEO accession.");
  }
  if (overlap.length) {
    flags.push(`training_validation_overlap:${overlap.join("|")}`);
    actions.push("Do not use overlapping training/validation cohorts as external validation.");
  }
  if (row.Manual_check_needed === "yes") {
    flags.push("manual_check_needed");
    actions.push("Inspect original PDF/supplement evidence context.");
  }

  const severe = flags.some((f) =>
    /formula_missing|coefficients_missing|suspicious_gene|gene_coef_count|pending_supplement|missing_pmid|suspicious_title|invalid_gse/.test(f),
  );
  const priority = severe ? "High" : flags.length ? "Medium" : "Low";
  return {
    QC_priority: priority,
    QC_flags: flags.join("; "),
    Suggested_action: [...new Set(actions)].join(" "),
    Signature_ID: `${row.Article_ID}${row.PMID ? `_PMID${String(row.PMID).replace(/\.0$/, "")}` : ""}`,
    Article_ID: row.Article_ID,
    PMID: row.PMID,
    Title: row.Title,
    Screening_level: row.Screening_level,
    Screening_label: row.Screening_label,
    Gene_number: row.Gene_number,
    Gene_list: row.Gene_list,
    Coefficients_available: row.Coefficients_available,
    Coefficients: row.Coefficients,
    Formula_available: row.Formula_available,
    Risk_score_formula: row.Risk_score_formula,
    Cutoff_available: row.Cutoff_available,
    Cutoff_value: row.Cutoff_value,
    Cutoff_type: row.Cutoff_type,
    Training_cohort: row.Training_cohort,
    Validation_cohort: row.Validation_cohort,
    GEO_accession: row.GEO_accession,
    Manual_check_needed: row.Manual_check_needed,
    Evidence_formula_context: row.Evidence_formula_context,
    Source_files: row.Source_files,
  };
}

function sortPriority(a, b) {
  const score = { High: 0, Medium: 1, Low: 2 };
  return score[a.QC_priority] - score[b.QC_priority] || String(a.Article_ID).localeCompare(String(b.Article_ID));
}

async function addSheet(wb, name, rows, fields) {
  const sheet = wb.worksheets.add(name);
  const data = [fields, ...rows.map((r) => fields.map((f) => sanitize(r[f])))];
  const endCol = colLetter(fields.length);
  const endRow = data.length;
  sheet.getRange(`A1:${endCol}${endRow}`).values = data;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${endCol}1`).format = {
    fill: "#244062",
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
  };
  sheet.getRange(`A1:${endCol}${endRow}`).format = { wrapText: true };
  fields.forEach((field, idx) => {
    let width = 18;
    if (["Title", "Suggested_action", "QC_flags", "Risk_score_formula", "Evidence_formula_context", "Source_files"].includes(field)) width = 48;
    if (["Gene_list", "Coefficients", "Training_cohort", "Validation_cohort", "GEO_accession"].includes(field)) width = 38;
    if (["QC_priority", "PMID", "Article_ID", "Gene_number"].includes(field)) width = 14;
    sheet.getRange(`${colLetter(idx + 1)}:${colLetter(idx + 1)}`).format.columnWidth = width;
  });
  try {
    sheet.tables.add(`A1:${endCol}${endRow}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {
    // The range still exports cleanly if table styling is skipped.
  }
  return sheet;
}

await fs.mkdir(outputDir, { recursive: true });
const rows = JSON.parse(await fs.readFile(inputJson, "utf8"));
const candidates = rows
  .filter((r) => ["Included_reproducible", "Pending_supplement"].includes(r.Screening_label))
  .map(qcRow)
  .sort(sortPriority);

const fields = [
  "QC_priority",
  "QC_flags",
  "Suggested_action",
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Title",
  "Screening_level",
  "Screening_label",
  "Gene_number",
  "Gene_list",
  "Coefficients_available",
  "Coefficients",
  "Formula_available",
  "Risk_score_formula",
  "Cutoff_available",
  "Cutoff_value",
  "Cutoff_type",
  "Training_cohort",
  "Validation_cohort",
  "GEO_accession",
  "Manual_check_needed",
  "Evidence_formula_context",
  "Source_files",
];

const summaryRows = ["High", "Medium", "Low"].map((p) => ({
  QC_priority: p,
  Count: candidates.filter((r) => r.QC_priority === p).length,
}));

const wb = Workbook.create();
await addSheet(wb, "QC_Priority", candidates, fields);
await addSheet(wb, "QC_Summary", summaryRows, ["QC_priority", "Count"]);
await wb.inspect({
  kind: "table",
  range: `QC_Priority!A1:H${Math.min(candidates.length + 1, 12)}`,
  tableMaxRows: 12,
  tableMaxCols: 8,
  maxChars: 2000,
});
try {
  await wb.render({ sheetName: "QC_Priority", range: "A1:H12", scale: 1, format: "png" });
} catch {
  // Export is the required artifact; render is a best-effort visual pass.
}
const output = await SpreadsheetFile.exportXlsx(wb);
await output.save(outputPath);
console.log(
  JSON.stringify(
    {
      outputPath,
      rows: candidates.length,
      summary: summaryRows,
    },
    null,
    2,
  ),
);
