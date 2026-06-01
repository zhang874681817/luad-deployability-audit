import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_GEO_manuscript_ready_signature_selection.xlsx");

function parseCsv(text) {
  const lines = text.replace(/^\uFEFF/, "").trimEnd().split(/\r?\n/);
  if (!lines.length || !lines[0]) return [];
  const parseLine = (line) => {
    const out = [];
    let cur = "";
    let quoted = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (ch === '"') {
        if (quoted && line[i + 1] === '"') {
          cur += '"';
          i++;
        } else quoted = !quoted;
      } else if (ch === "," && !quoted) {
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

function numericField(field) {
  return /PMID|cohorts|p05|rate|C_index|HR|count|N$|_n$|Events|Cutoff|Delta|Median_|Best_|Gene_/i.test(field);
}

function value(field, v) {
  if (v === null || v === undefined || v === "") return "";
  if (numericField(field)) {
    const n = Number(v);
    if (Number.isFinite(n)) return n;
  }
  return String(v).slice(0, 32000);
}

async function addSheet(wb, name, rows, fields, fill = "#1F4E79") {
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
    const width = /Title|Role|caveat|Reason|cohort|Validation|Notes|GEOs|genes|decision/i.test(field) ? 52 : 17;
    sheet.getRange(`${colLetter(i + 1)}:${colLetter(i + 1)}`).format.columnWidth = width;
    if (/p05|C_index|HR|Cutoff|Median|Best|Delta|rate/i.test(field)) {
      sheet.getRange(`${colLetter(i + 1)}2:${colLetter(i + 1)}${endRow}`).numberFormat = "0.000";
    }
  });
  try {
    sheet.tables.add(`A1:${endCol}${endRow}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {}
}

function tierRows(rows, tier) {
  return rows.filter((r) => r.Final_analysis_tier === tier);
}

await fs.mkdir(outputDir, { recursive: true });

const selection = await readCsv(path.join(root, "work", "luad_geo_manuscript_ready_signature_selection.csv"));
const a097 = await readCsv(path.join(root, "work", "a097_fixed_cutoff_1p66_stress_test.csv"));
const override = await readCsv(path.join(root, "work", "geo_external_validation_manual_override_comparison.csv"));

const tierCounts = Object.values(
  selection.reduce((acc, row) => {
    const tier = row.Final_analysis_tier;
    acc[tier] ||= { Final_analysis_tier: tier, Signature_count: 0 };
    acc[tier].Signature_count++;
    return acc;
  }, {}),
);

const fields = [
  "Final_analysis_tier",
  "Recommended_manuscript_role",
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Tested_cohorts_ok",
  "Group_sig_p05_n",
  "Novel_external_group_sig_p05_n",
  "Continuous_sig_p05_n",
  "Median_C_index",
  "Median_group_HR",
  "Reported_numeric_cutoff_tests",
  "Cohort_median_fallback_tests",
  "Cutoff_type",
  "Cutoff_value",
  "Manual_decision",
  "Key_caveats_for_manuscript",
  "Title",
];

const workbook = Workbook.create();

await addSheet(
  workbook,
  "Selection_Summary",
  [
    { Metric: "Strict reported-cutoff primary", Value: "A131_PMID39735553" },
    { Metric: "Primary median-fallback count", Value: tierRows(selection, "Tier_2_primary_median_fallback").length },
    { Metric: "Secondary median-fallback count", Value: tierRows(selection, "Tier_3_secondary_median_fallback").length },
    { Metric: "Continuous-score support count", Value: tierRows(selection, "Tier_4_continuous_score_support").length },
    { Metric: "Manual override sensitivity count", Value: tierRows(selection, "Tier_5_manual_override_sensitivity").length },
    {
      Metric: "A097 decision",
      Value:
        "PDF confirms cutoff 1.66, but fixed-cutoff GEO stress test collapses every scored cohort to one group; keep as median-fallback score portability.",
    },
    {
      Metric: "A010/A104 decision",
      Value:
        "Manual formula/symbol overrides make them testable; keep as sensitivity-only, not strict primary analysis.",
    },
  ],
  ["Metric", "Value"],
  "#0B5A6B",
);

await addSheet(workbook, "Tier_Counts", tierCounts, ["Final_analysis_tier", "Signature_count"], "#0B5A6B");
await addSheet(workbook, "All_Signatures", selection, fields, "#1F4E79");
await addSheet(workbook, "Strict_Cutoff", tierRows(selection, "Tier_1_strict_reported_cutoff"), fields, "#2E7D32");
await addSheet(workbook, "Primary_Median_Fallback", tierRows(selection, "Tier_2_primary_median_fallback"), fields, "#7A5195");
await addSheet(workbook, "Secondary_Sensitivity", [
  ...tierRows(selection, "Tier_3_secondary_median_fallback"),
  ...tierRows(selection, "Tier_4_continuous_score_support"),
  ...tierRows(selection, "Tier_5_manual_override_sensitivity"),
], fields, "#355C7D");
await addSheet(workbook, "Do_Not_Advance", [
  ...tierRows(selection, "Do_not_advance_first_pass"),
  ...tierRows(selection, "Not_testable_first_pass"),
], fields, "#8F3A3A");

await addSheet(
  workbook,
  "A097_Fixed_Cutoff",
  a097,
  [
    "Signature_ID",
    "Fixed_cutoff_from_uploaded_pdf",
    "GSE",
    "N_scored",
    "Score_min",
    "Score_median",
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
  "Manual_Overrides",
  override.filter((r) => r.Manual_override_applied === "yes"),
  [
    "Signature_ID",
    "Tested_cohorts_override",
    "Significant_cohorts_p05_override",
    "Median_C_index_override",
    "Median_group_HR_override",
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
  sheetName: "All_Signatures",
  range: "A1:R25",
  scale: 1,
  format: "png",
});
await fs.writeFile(
  path.join(outputDir, "LUAD_GEO_manuscript_ready_signature_selection_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
console.log(`Wrote ${outputPath}`);
