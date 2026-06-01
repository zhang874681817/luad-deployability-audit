import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_GEO_external_validation_manuscript_review.xlsx");

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
        } else {
          quoted = !quoted;
        }
      } else if (ch === "," && !quoted) {
        out.push(cur);
        cur = "";
      } else {
        cur += ch;
      }
    }
    out.push(cur);
    return out;
  };
  const header = parseLine(lines[0]);
  return lines.slice(1).filter(Boolean).map((line) => {
    const vals = parseLine(line);
    const row = {};
    header.forEach((h, i) => {
      row[h] = vals[i] ?? "";
    });
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

function num(v) {
  if (v === null || v === undefined || v === "") return "";
  const n = Number(v);
  return Number.isFinite(n) ? n : v;
}

function value(field, v) {
  if (v === null || v === undefined) return "";
  if (
    /n$|count|rate|C_index|HR|p05|tests|rows|events|PMID|Gene_count|Cutoff_value|Cox_|Logrank|N$|High_|Low_|Score_|Coverage|Cutoff_used/i.test(
      field,
    )
  ) {
    return num(v);
  }
  const s = String(v);
  return s.length > 32000 ? `${s.slice(0, 31950)} [truncated]` : s;
}

function widthFor(field) {
  if (/Title|Notes|Details|cohort|GEOs|genes|Overlap|Validation/i.test(field)) return 48;
  if (/Recommendation|Flag|Cutoff|status|role|Signature_ID/i.test(field)) return 26;
  if (/rate|C_index|HR|p05|tests|rows|n$|PMID|Events/i.test(field)) return 14;
  return 18;
}

async function addSheet(wb, name, rows, fields, fill = "#1F4E79") {
  const sheet = wb.worksheets.add(name);
  const data = [fields, ...rows.map((r) => fields.map((f) => value(f, r[f])))];
  const endCol = colLetter(fields.length);
  const endRow = Math.max(1, data.length);
  sheet.getRange(`A1:${endCol}${endRow}`).values = data;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${endCol}1`).format = {
    fill,
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
  };
  sheet.getRange(`A1:${endCol}${endRow}`).format = { wrapText: true };
  fields.forEach((field, i) => {
    sheet.getRange(`${colLetter(i + 1)}:${colLetter(i + 1)}`).format.columnWidth = widthFor(field);
  });
  fields.forEach((field, i) => {
    if (/rate|C_index|HR|p05|p$|Cutoff_value|Cutoff_used|Score_|Coverage/i.test(field)) {
      sheet.getRange(`${colLetter(i + 1)}2:${colLetter(i + 1)}${endRow}`).numberFormat = "0.000";
    }
  });
  try {
    sheet.tables.add(`A1:${endCol}${endRow}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {}
  return sheet;
}

function countBy(rows, field) {
  const map = new Map();
  for (const row of rows) {
    const key = row[field] || "";
    map.set(key, (map.get(key) || 0) + 1);
  }
  return [...map.entries()].map(([key, count]) => ({ [field]: key, Count: count }));
}

await fs.mkdir(outputDir, { recursive: true });

const sigRows = await readCsv(path.join(root, "work", "geo_external_validation_manuscript_signature_review.csv"));
const pairRows = await readCsv(path.join(root, "work", "geo_external_validation_manuscript_pairwise_review.csv"));
const flagRows = await readCsv(path.join(root, "work", "geo_external_validation_manuscript_manual_audit_flags.csv"));

const reviewSummaryRows = [
  { Metric: "Signatures reviewed", Value: sigRows.length },
  { Metric: "Successful pairwise GEO tests", Value: pairRows.length },
  { Metric: "Primary numeric-cutoff candidates", Value: sigRows.filter((r) => r.Recommendation_Tier === "Primary_numeric_cutoff_candidate").length },
  { Metric: "Primary median-fallback candidates", Value: sigRows.filter((r) => r.Recommendation_Tier === "Primary_median_fallback_candidate").length },
  { Metric: "Secondary median-fallback candidates", Value: sigRows.filter((r) => r.Recommendation_Tier === "Secondary_median_fallback_candidate").length },
  { Metric: "Continuous-score-only numeric-cutoff candidates", Value: sigRows.filter((r) => r.Recommendation_Tier === "Continuous_score_only_numeric_cutoff_not_group_ready").length },
  { Metric: "Signatures with reported numeric cutoff tests", Value: sigRows.filter((r) => Number(r.Reported_numeric_cutoff_tests) > 0).length },
  { Metric: "Signatures where all grouped tests use median fallback", Value: sigRows.filter((r) => r.All_group_tests_use_median_fallback === "yes").length },
  { Metric: "Manual audit flags", Value: flagRows.length },
  {
    Metric: "Recommended main manuscript split",
    Value:
      "Use reported numeric cutoffs as the strict grouped-risk analysis; present median-fallback signatures separately as score portability or sensitivity analyses.",
  },
  {
    Metric: "Mandatory caveat",
    Value:
      "GSE50081 survival_time was converted from years to months because the GEO field lacks explicit unit; verify before final tables.",
  },
];

const workbook = Workbook.create();
await addSheet(workbook, "Review_Summary", reviewSummaryRows, ["Metric", "Value"], "#0B5A6B");
await addSheet(workbook, "Tier_Counts", countBy(sigRows, "Recommendation_Tier"), ["Recommendation_Tier", "Count"], "#0B5A6B");

const candidateFields = [
  "Recommendation_Tier",
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Tested_cohorts_ok",
  "Group_sig_p05_n",
  "Novel_external_group_sig_p05_n",
  "Continuous_sig_p05_n",
  "Median_C_index",
  "Median_group_HR",
  "Group_HR_gt1_rate",
  "Continuous_HR_gt1_rate",
  "Reported_numeric_cutoff_tests",
  "Cohort_median_fallback_tests",
  "All_group_tests_use_median_fallback",
  "Original_validation_GSEs_tested",
  "GSE50081_in_ok_tests",
  "Cutoff_type",
  "Cutoff_value",
  "Manuscript_Use_Notes",
  "Title",
];
await addSheet(workbook, "Signature_Review", sigRows, candidateFields, "#1F4E79");

const pairFields = [
  "Signature_ID",
  "GSE",
  "Cohort_role_for_signature",
  "Cutoff_rule_used",
  "Cutoff_type_original",
  "N",
  "Events",
  "C_index",
  "Cox_continuous_HR_per_SD",
  "Cox_continuous_p",
  "Cox_high_vs_low_HR",
  "Cox_high_vs_low_p",
  "Logrank_p",
  "Group_p_lt_0_05",
  "Continuous_p_lt_0_05",
  "GSE50081_time_unit_flag",
  "Title",
];
await addSheet(workbook, "Pairwise_OK_Reviewed", pairRows, pairFields, "#355C7D");

const fallbackRows = sigRows.filter((r) => r.All_group_tests_use_median_fallback === "yes" && Number(r.Group_sig_p05_n) > 0);
await addSheet(workbook, "Median_Fallback_Dependent", fallbackRows, candidateFields, "#7A5195");

await addSheet(
  workbook,
  "Manual_Audit_Flags",
  flagRows,
  ["Severity", "Flag_type", "Signature_ID", "Details"],
  "#8F3A3A",
);

const inspect = await workbook.inspect({
  kind: "sheet,table",
  tableMaxRows: 4,
  tableMaxCols: 8,
  maxChars: 5000,
});
console.log(inspect.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "formula error scan",
});
console.log(errors.ndjson);

const preview = await workbook.render({
  sheetName: "Signature_Review",
  range: "A1:U20",
  scale: 1,
  format: "png",
});
await fs.writeFile(
  path.join(outputDir, "LUAD_GEO_external_validation_manuscript_review_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
console.log(`Wrote ${outputPath}`);
