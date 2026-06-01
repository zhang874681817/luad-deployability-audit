import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_GEO_external_validation_results.xlsx");

async function readCsv(file) {
  const text = await fs.readFile(file, "utf8");
  const lines = text.replace(/^\uFEFF/, "").trimEnd().split(/\r?\n/);
  if (!lines.length || !lines[0]) return [];
  const parse = (line) => {
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
  const header = parse(lines[0]);
  return lines.slice(1).filter(Boolean).map((line) => {
    const vals = parse(line);
    const row = {};
    header.forEach((h, i) => (row[h] = vals[i] ?? ""));
    return row;
  });
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
  return Number.isFinite(n) ? n : "";
}

function sanitize(v) {
  if (v === null || v === undefined) return "";
  if (typeof v === "number") return Number.isFinite(v) ? v : "";
  const s = String(v);
  return s.length > 32000 ? `${s.slice(0, 31950)} [truncated]` : s;
}

const numericFields = new Set([
  "PMID",
  "Gene_count",
  "Present_gene_count",
  "Missing_gene_count",
  "Gene_coverage",
  "N",
  "Events",
  "Cutoff_used",
  "High_n",
  "Low_n",
  "High_events",
  "Low_events",
  "Score_mean",
  "Score_sd",
  "Score_median",
  "C_index",
  "Cox_continuous_HR_per_SD",
  "Cox_continuous_CI95_low",
  "Cox_continuous_CI95_high",
  "Cox_continuous_p",
  "Cox_high_vs_low_HR",
  "Cox_high_vs_low_CI95_low",
  "Cox_high_vs_low_CI95_high",
  "Cox_high_vs_low_p",
  "Logrank_p",
  "Tested_cohorts",
  "Significant_cohorts_p05",
  "Median_C_index",
  "Best_C_index",
  "Median_group_HR",
  "OS_time_months",
  "OS_event",
  "Risk_score",
]);

function cellValue(field, value) {
  return numericFields.has(field) ? num(value) : sanitize(value);
}

function widthFor(field) {
  if (/Title|Training|Validation|Missing_genes|Formula|Notes/i.test(field)) return 52;
  if (/Signature_ID|Article_ID|Cutoff_rule|status|transform/i.test(field)) return 24;
  if (/p$|HR|C_index|CI95|score|cutoff|events|cohort|coverage|count|PMID|GSE|N$/i.test(field)) return 16;
  return 20;
}

async function addSheet(wb, name, rows, fields, fill = "#1F4E79") {
  const sheet = wb.worksheets.add(name);
  const data = [fields, ...rows.map((r) => fields.map((f) => cellValue(f, r[f])))];
  const endCol = colLetter(fields.length);
  const endRow = Math.max(data.length, 1);
  sheet.getRange(`A1:${endCol}${endRow}`).values = data;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${endCol}1`).format = {
    fill,
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
  };
  sheet.getRange(`A1:${endCol}${endRow}`).format = { wrapText: true };
  fields.forEach((f, i) => {
    sheet.getRange(`${colLetter(i + 1)}:${colLetter(i + 1)}`).format.columnWidth = widthFor(f);
  });
  for (const [i, f] of fields.entries()) {
    if (/p$|C_index|coverage|HR|score|cutoff/i.test(f)) {
      sheet.getRange(`${colLetter(i + 1)}2:${colLetter(i + 1)}${endRow}`).numberFormat = "0.000";
    }
  }
  try {
    sheet.tables.add(`A1:${endCol}${endRow}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {}
}

function uniqueCount(rows, field) {
  return new Set(rows.map((r) => r[field]).filter(Boolean)).size;
}

await fs.mkdir(outputDir, { recursive: true });

const results = await readCsv(path.join(root, "work", "geo_external_validation_results.csv"));
const scores = await readCsv(path.join(root, "work", "geo_external_validation_sample_scores.csv"));
const signatureSummaryRaw = await readCsv(path.join(root, "work", "geo_external_validation_signature_summary.csv"));

const ok = results.filter((r) => r.status === "ok");
const sigSummary = signatureSummaryRaw
  .map((r) => ({
    ...r,
    Significant_rate: num(r.Tested_cohorts) ? num(r.Significant_cohorts_p05) / num(r.Tested_cohorts) : "",
  }))
  .sort((a, b) => {
    const d1 = num(b.Significant_cohorts_p05) - num(a.Significant_cohorts_p05);
    if (d1) return d1;
    return num(b.Best_C_index) - num(a.Best_C_index);
  });

const statusCounts = Object.values(
  results.reduce((acc, r) => {
    acc[r.status] ||= { Status: r.status, Pairwise_rows: 0 };
    acc[r.status].Pairwise_rows++;
    return acc;
  }, {}),
).sort((a, b) => b.Pairwise_rows - a.Pairwise_rows);

const cohortSummary = Object.values(
  ok.reduce((acc, r) => {
    const g = r.GSE;
    acc[g] ||= {
      GSE: g,
      Tested_signatures: 0,
      Significant_high_vs_low_p05: 0,
      Median_C_index: [],
      Median_HR: [],
      Total_scored_samples: 0,
      Total_events_in_tests: 0,
    };
    const o = acc[g];
    o.Tested_signatures++;
    if (num(r.Cox_high_vs_low_p) !== "" && num(r.Cox_high_vs_low_p) < 0.05) o.Significant_high_vs_low_p05++;
    if (num(r.C_index) !== "") o.Median_C_index.push(num(r.C_index));
    if (num(r.Cox_high_vs_low_HR) !== "") o.Median_HR.push(num(r.Cox_high_vs_low_HR));
    o.Total_scored_samples += num(r.N) || 0;
    o.Total_events_in_tests += num(r.Events) || 0;
    return acc;
  }, {}),
).map((r) => {
  const median = (xs) => {
    const s = xs.filter((x) => Number.isFinite(x)).sort((a, b) => a - b);
    if (!s.length) return "";
    const mid = Math.floor(s.length / 2);
    return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
  };
  return {
    ...r,
    Median_C_index: median(r.Median_C_index),
    Median_HR: median(r.Median_HR),
    Significant_rate: r.Tested_signatures ? r.Significant_high_vs_low_p05 / r.Tested_signatures : "",
  };
}).sort((a, b) => b.Tested_signatures - a.Tested_signatures || b.Significant_high_vs_low_p05 - a.Significant_high_vs_low_p05);

const summaryRows = [
  { Metric: "Pairwise_signature_GEO_rows", Value: results.length },
  { Metric: "Successful_external_tests", Value: ok.length },
  { Metric: "Signatures_with_successful_GEO_tests", Value: uniqueCount(ok, "Signature_ID") },
  { Metric: "GEO_cohorts_used", Value: uniqueCount(ok, "GSE") },
  { Metric: "Sample_level_scores", Value: scores.length },
  { Metric: "High_vs_low_tests_p_lt_0.05", Value: ok.filter((r) => num(r.Cox_high_vs_low_p) !== "" && num(r.Cox_high_vs_low_p) < 0.05).length },
  { Metric: "Cutoff_policy", Value: "Original numeric cutoff is used when available; otherwise cohort median fallback is used and flagged in Cutoff_rule_used." },
  { Metric: "Overlap_policy", Value: "GEO cohorts named in the original training cohort are skipped for that signature." },
  { Metric: "Gene_coverage_policy", Value: "Only full gene coverage is scored; incomplete coverage rows remain in Pairwise_Results with status." },
  { Metric: "Caveat", Value: "Automated first-pass validation; verify cutoff interpretation, endpoint compatibility, and GSE50081 time-unit conversion before manuscript use." },
];

const methodsRows = [
  { Step: "Input signatures", Detail: "Included_reproducible signatures with parsed gene coefficients from the full-text extraction workflow." },
  { Step: "Expression", Detail: "Gene-level GEO matrices collapsed by gene symbol; unlogged-looking matrices use log2(x+1) scoring copies." },
  { Step: "Survival endpoint", Detail: "LUAD-like samples with complete overall survival time and event fields." },
  { Step: "Risk score", Detail: "Sum of coefficient times gene expression for all genes in the signature." },
  { Step: "Cutoff", Detail: "Use numeric extracted cutoff when available; otherwise use the current external cohort median to create high/low groups." },
  { Step: "Statistics", Detail: "Univariate Cox model for continuous score per SD and high vs low risk group; log-rank p and Harrell-style C-index are computed in script." },
  { Step: "Training overlap", Detail: "A GEO cohort is skipped for a signature when that GSE appears in the original training cohort field." },
];

const wb = Workbook.create();
await addSheet(wb, "Summary", summaryRows, ["Metric", "Value"], "#2F5597");
await addSheet(wb, "Signature_Summary", sigSummary, [
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Title",
  "Tested_cohorts",
  "Significant_cohorts_p05",
  "Significant_rate",
  "Median_C_index",
  "Best_C_index",
  "Median_group_HR",
], "#548235");
await addSheet(wb, "Cohort_Summary", cohortSummary, [
  "GSE",
  "Tested_signatures",
  "Significant_high_vs_low_p05",
  "Significant_rate",
  "Median_C_index",
  "Median_HR",
  "Total_scored_samples",
  "Total_events_in_tests",
], "#7F6000");
await addSheet(wb, "Status_Counts", statusCounts, ["Status", "Pairwise_rows"], "#8064A2");
await addSheet(wb, "Pairwise_Results", results, [
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Title",
  "GSE",
  "Training_cohort",
  "Validation_cohort",
  "Gene_count",
  "Present_gene_count",
  "Missing_gene_count",
  "Gene_coverage",
  "Missing_genes",
  "N",
  "Events",
  "Scoring_transform",
  "status",
  "Cutoff_type_original",
  "Cutoff_used",
  "Cutoff_rule_used",
  "High_n",
  "Low_n",
  "High_events",
  "Low_events",
  "Score_mean",
  "Score_sd",
  "Score_median",
  "C_index",
  "Cox_continuous_HR_per_SD",
  "Cox_continuous_CI95_low",
  "Cox_continuous_CI95_high",
  "Cox_continuous_p",
  "Cox_high_vs_low_HR",
  "Cox_high_vs_low_CI95_low",
  "Cox_high_vs_low_CI95_high",
  "Cox_high_vs_low_p",
  "Logrank_p",
], "#1F4E79");
await addSheet(wb, "Sample_Scores", scores, [
  "Signature_ID",
  "GSE",
  "geo_accession",
  "OS_time_months",
  "OS_event",
  "Risk_score",
  "Risk_group",
  "Cutoff_used",
], "#375623");
await addSheet(wb, "Methods_Notes", methodsRows, ["Step", "Detail"], "#7030A0");

await wb.inspect({ kind: "table", range: "Summary!A1:B11", tableMaxRows: 11, tableMaxCols: 2, maxChars: 2000 });
await wb.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 50 }, summary: "formula error scan" });
try {
  await wb.render({ sheetName: "Signature_Summary", range: "A1:J16", scale: 1, format: "png" });
  await wb.render({ sheetName: "Pairwise_Results", range: "A1:L16", scale: 1, format: "png" });
  await wb.render({ sheetName: "Sample_Scores", range: "A1:H16", scale: 1, format: "png" });
} catch {}

const out = await SpreadsheetFile.exportXlsx(wb);
await out.save(outputPath);
console.log(outputPath);
