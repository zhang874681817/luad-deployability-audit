import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_GEO_validation_manual_decisions.xlsx");

function parseCsv(text) {
  const lines = text.replace(/^\uFEFF/, "").trimEnd().split(/\r?\n/);
  if (!lines.length || !lines[0]) return [];
  const parseLine = (line) => {
    const out = [];
    let cur = "";
    let q = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (ch === '"') {
        if (q && line[i + 1] === '"') {
          cur += '"';
          i++;
        } else q = !q;
      } else if (ch === "," && !q) {
        out.push(cur);
        cur = "";
      } else cur += ch;
    }
    out.push(cur);
    return out;
  };
  const header = parseLine(lines[0]);
  return lines.slice(1).filter(Boolean).map((line) => {
    const vals = parseLine(line);
    const row = {};
    header.forEach((h, i) => (row[h] = vals[i] ?? ""));
    return row;
  });
}

async function readCsv(file) {
  return parseCsv(await fs.readFile(file, "utf8"));
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

function value(field, v) {
  if (v === null || v === undefined || v === "") return "";
  if (/PMID|cohorts|p05|C_index|HR|count|N_|_n$|Score_|Cutoff|Delta|Events|Gene_|High_|Low_|p75|p25|min|max|median/i.test(field)) {
    const n = Number(v);
    if (Number.isFinite(n)) return n;
  }
  return String(v).slice(0, 32000);
}

async function addSheet(wb, name, rows, fields, fill) {
  const sheet = wb.worksheets.add(name);
  const data = [fields, ...rows.map((r) => fields.map((f) => value(f, r[f])))];
  const endCol = colLetter(fields.length);
  const endRow = data.length;
  sheet.getRange(`A1:${endCol}${endRow}`).values = data;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${endCol}1`).format = {
    fill,
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
  };
  sheet.getRange(`A1:${endCol}${endRow}`).format = { wrapText: true };
  fields.forEach((field, i) => {
    const width = /Title|Reason|Interpretation|cohort|Validation|Notes|GEOs|genes/i.test(field) ? 52 : 17;
    sheet.getRange(`${colLetter(i + 1)}:${colLetter(i + 1)}`).format.columnWidth = width;
    if (/p05|C_index|HR|Cutoff|Score_|Delta|min|max|median|p75|p25/i.test(field)) {
      sheet.getRange(`${colLetter(i + 1)}2:${colLetter(i + 1)}${endRow}`).numberFormat = "0.000";
    }
  });
  try {
    sheet.tables.add(`A1:${endCol}${endRow}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {}
}

await fs.mkdir(outputDir, { recursive: true });

const decisions = await readCsv(path.join(root, "work", "geo_external_validation_manuscript_signature_review.manual_decisions.csv"));
const a097 = await readCsv(path.join(root, "work", "a097_fixed_cutoff_1p66_stress_test.csv"));
const overrides = await readCsv(path.join(root, "work", "geo_external_validation_manual_override_comparison.csv"));

const manualRows = decisions.filter((r) => r.Manual_decision);
const primaryStrict = decisions.filter((r) => r.Recommendation_Tier === "Primary_numeric_cutoff_candidate");
const primaryFallback = decisions.filter((r) => r.Recommendation_Tier === "Primary_median_fallback_candidate");

const workbook = Workbook.create();
await addSheet(
  workbook,
  "Decision_Summary",
  [
    { Metric: "Strict reported-cutoff primary candidate", Value: primaryStrict.map((r) => r.Signature_ID).join("; ") },
    {
      Metric: "Primary median-fallback candidates",
      Value: primaryFallback.map((r) => r.Signature_ID).join("; "),
    },
    {
      Metric: "A097 decision",
      Value:
        "PDF confirms TCGA/training cutoff 1.66, but external GEO fixed-cutoff stress test collapses every scored cohort to one group; keep as median-fallback score portability.",
    },
    {
      Metric: "A010/A104 decision",
      Value:
        "Manual formula/symbol overrides make them testable, but both remain sensitivity/secondary median-fallback candidates.",
    },
  ],
  ["Metric", "Value"],
  "#0B5A6B",
);

await addSheet(
  workbook,
  "Signature_Manual_Decisions",
  decisions,
  [
    "Recommendation_Tier",
    "Signature_ID",
    "Article_ID",
    "PMID",
    "Tested_cohorts_ok",
    "Group_sig_p05_n",
    "Novel_external_group_sig_p05_n",
    "Continuous_sig_p05_n",
    "Median_C_index",
    "Reported_numeric_cutoff_tests",
    "Cohort_median_fallback_tests",
    "Cutoff_type",
    "Cutoff_value",
    "Manual_decision",
    "Manual_decision_reason",
    "Title",
  ],
  "#1F4E79",
);

await addSheet(
  workbook,
  "A097_Fixed_Cutoff_Test",
  a097,
  [
    "Signature_ID",
    "Fixed_cutoff_from_uploaded_pdf",
    "GSE",
    "N_scored",
    "Score_min",
    "Score_p25",
    "Score_median",
    "Score_p75",
    "Score_max",
    "Low_n_score_lt_1p66",
    "High_n_score_ge_1p66",
    "Fixed_cutoff_stress_status",
    "Interpretation",
  ],
  "#7A5195",
);

await addSheet(
  workbook,
  "A010_A104_Override_Compare",
  overrides.filter((r) => r.Manual_override_applied === "yes"),
  [
    "Manual_override_applied",
    "Signature_ID",
    "Article_ID",
    "PMID",
    "Tested_cohorts_override",
    "Significant_cohorts_p05_override",
    "Median_C_index_override",
    "Median_group_HR_override",
    "Delta_Tested_cohorts",
    "Delta_Significant_cohorts_p05",
    "Cutoff_type",
    "Cutoff_value",
    "Interpretation",
    "Title",
  ],
  "#355C7D",
);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
});
console.log(errors.ndjson);

const preview = await workbook.render({
  sheetName: "A097_Fixed_Cutoff_Test",
  range: "A1:M9",
  scale: 1,
  format: "png",
});
await fs.writeFile(
  path.join(outputDir, "LUAD_GEO_validation_manual_decisions_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
console.log(`Wrote ${outputPath}`);
