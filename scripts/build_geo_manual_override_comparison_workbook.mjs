import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_GEO_manual_override_comparison.xlsx");

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
  if (/PMID|cohorts|p05|C_index|HR|Delta|count|N$|Events|Gene_|High_|Low_|Score_|Cutoff/i.test(field)) {
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
    const width = /Title|Interpretation|cohort|Validation|Resolved|Missing/i.test(field) ? 52 : 17;
    sheet.getRange(`${colLetter(i + 1)}:${colLetter(i + 1)}`).format.columnWidth = width;
    if (/p05|C_index|HR|Delta|Cutoff/i.test(field)) {
      sheet.getRange(`${colLetter(i + 1)}2:${colLetter(i + 1)}${endRow}`).numberFormat = "0.000";
    }
  });
  try {
    sheet.tables.add(`A1:${endCol}${endRow}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {}
}

await fs.mkdir(outputDir, { recursive: true });

const comparison = await readCsv(path.join(root, "work", "geo_external_validation_manual_override_comparison.csv"));
const results = await readCsv(path.join(root, "work", "geo_external_validation_results.manual_overrides.csv"));
const overrideRows = results.filter((r) => r.Override_note === "manual_formula_symbol_override");

const workbook = Workbook.create();
await addSheet(
  workbook,
  "Override_Summary",
  [
    {
      Metric: "Purpose",
      Value:
        "Sensitivity validation after manual formula/symbol overrides for A010 and A104. Original first-pass outputs are unchanged.",
    },
    { Metric: "A010 result", Value: "5 tested non-training GEO cohorts; 3 significant; median C-index 0.589." },
    { Metric: "A104 result", Value: "5 tested non-training GEO cohorts; 4 significant; median C-index 0.575." },
    {
      Metric: "Interpretation",
      Value:
        "Both are sensitivity/secondary median-fallback candidates, not strict reported-cutoff manuscript candidates.",
    },
  ],
  ["Metric", "Value"],
  "#0B5A6B",
);

await addSheet(
  workbook,
  "Signature_Comparison",
  comparison,
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
    "Delta_Median_C_index",
    "Cutoff_type",
    "Cutoff_value",
    "Interpretation",
    "Title",
  ],
  "#1F4E79",
);

await addSheet(
  workbook,
  "Override_Pairwise",
  overrideRows,
  [
    "Signature_ID",
    "GSE",
    "status",
    "Gene_count",
    "Present_gene_count",
    "Missing_gene_count",
    "Missing_genes",
    "Cutoff_rule_used",
    "N",
    "Events",
    "C_index",
    "Cox_high_vs_low_HR",
    "Cox_high_vs_low_p",
    "Cox_continuous_HR_per_SD",
    "Cox_continuous_p",
    "Resolved_gene_symbols",
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
  sheetName: "Signature_Comparison",
  range: "A1:O15",
  scale: 1,
  format: "png",
});
await fs.writeFile(
  path.join(outputDir, "LUAD_GEO_manual_override_comparison_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outputPath);
console.log(`Wrote ${outputPath}`);
