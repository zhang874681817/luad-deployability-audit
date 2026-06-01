import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_expression_matrix_inventory.xlsx");

async function readCsv(file) {
  const text = await fs.readFile(file, "utf8");
  const lines = text.replace(/^\uFEFF/, "").trimEnd().split(/\r?\n/);
  if (!lines.length) return [];
  const parseLine = (line) => {
    const out = [];
    let cur = "";
    let inQ = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (ch === '"') {
        if (inQ && line[i + 1] === '"') {
          cur += '"';
          i++;
        } else {
          inQ = !inQ;
        }
      } else if (ch === "," && !inQ) {
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

function sanitize(v) {
  if (v === null || v === undefined) return "";
  if (typeof v === "number") return Number.isFinite(v) ? v : "";
  const s = String(v);
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
  if (/file|path|Title|genes|note|transform|url/i.test(field)) return 48;
  if (/count|fraction|sample|feature|GPL|GSE|status/i.test(field)) return 16;
  return 22;
}

async function addSheet(wb, name, rows, fields) {
  const sheet = wb.worksheets.add(name);
  const data = [fields, ...rows.map((r) => fields.map((f) => sanitize(r[f])))];
  const endCol = colLetter(fields.length);
  const endRow = Math.max(data.length, 1);
  sheet.getRange(`A1:${endCol}${endRow}`).values = data;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${endCol}1`).format = { fill: "#5F497A", font: { bold: true, color: "#FFFFFF" }, wrapText: true };
  sheet.getRange(`A1:${endCol}${endRow}`).format = { wrapText: true };
  fields.forEach((f, i) => {
    sheet.getRange(`${colLetter(i + 1)}:${colLetter(i + 1)}`).format.columnWidth = widthFor(f);
  });
  try {
    sheet.tables.add(`A1:${endCol}${endRow}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {}
}

await fs.mkdir(outputDir, { recursive: true });
const geoSummary = await readCsv(path.join(root, "work", "geo_gene_expression_summary.csv"));
const gplSummary = await readCsv(path.join(root, "work", "geo_platform_annotation_summary.csv"));
const coverage = await readCsv(path.join(root, "work", "signature_geo_gene_coverage.csv"));

const usable = geoSummary.filter((r) => r.status === "ok");
const summary = [
  { Metric: "Downloaded_GEO_series", Value: geoSummary.length },
  { Metric: "Usable_gene_level_GEO_matrices", Value: usable.length },
  { Metric: "Empty_or_unusable_GEO_series", Value: geoSummary.filter((r) => r.status !== "ok").length },
  { Metric: "Total_GEO_samples_in_usable_matrices", Value: usable.reduce((a, r) => a + Number(r.sample_count || 0), 0) },
  { Metric: "Coverage_rows", Value: coverage.length },
  { Metric: "Note", Value: "Gene-level matrices are probe-collapsed by mean. Unlogged-looking matrices have log2(x+1) scoring copies." },
];

const covNums = coverage.map((r) => ({ ...r, Coverage_fraction_num: Number(r.Coverage_fraction || 0) }));
const coverageSummaryBySignature = Object.values(
  covNums.reduce((acc, r) => {
    const id = r.Signature_ID;
    acc[id] ||= {
      Signature_ID: id,
      Article_ID: r.Article_ID,
      PMID: r.PMID,
      Title: r.Title,
      Cohorts_tested: 0,
      Cohorts_full_coverage: 0,
      Cohorts_80pct_coverage: 0,
      Best_coverage_fraction: 0,
      Best_cohorts: "",
    };
    const o = acc[id];
    o.Cohorts_tested++;
    if (r.Coverage_fraction_num >= 1) o.Cohorts_full_coverage++;
    if (r.Coverage_fraction_num >= 0.8) o.Cohorts_80pct_coverage++;
    if (r.Coverage_fraction_num > o.Best_coverage_fraction) {
      o.Best_coverage_fraction = r.Coverage_fraction_num;
      o.Best_cohorts = r.GSE;
    } else if (r.Coverage_fraction_num === o.Best_coverage_fraction) {
      o.Best_cohorts += `; ${r.GSE}`;
    }
    return acc;
  }, {}),
);

const wb = Workbook.create();
await addSheet(wb, "Summary", summary, ["Metric", "Value"]);
await addSheet(wb, "GEO_Matrix_Inventory", geoSummary, [
  "GSE",
  "status",
  "GPL",
  "Mapping_note",
  "sample_count",
  "probe_feature_count",
  "gene_symbol_count",
  "scale_guess",
  "scoring_transform",
  "expression_file",
  "gene_expression_file",
  "scoring_expression_file",
]);
await addSheet(wb, "Platform_Annotations", gplSummary, [
  "GPL",
  "status",
  "annotation_source",
  "symbol_column",
  "probe_count",
  "mapped_probe_rows",
  "unique_symbols",
  "mapping_file",
  "annot_url",
  "error",
]);
await addSheet(wb, "Coverage_By_Signature", coverageSummaryBySignature, [
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Title",
  "Cohorts_tested",
  "Cohorts_full_coverage",
  "Cohorts_80pct_coverage",
  "Best_coverage_fraction",
  "Best_cohorts",
]);
await addSheet(wb, "Coverage_Long", coverage, [
  "Signature_ID",
  "Article_ID",
  "PMID",
  "Title",
  "GSE",
  "GPL",
  "Signature_gene_count",
  "Present_gene_count",
  "Missing_gene_count",
  "Coverage_fraction",
  "Missing_genes",
  "Scoring_expression_file",
  "Scoring_transform",
]);

await wb.inspect({ kind: "table", range: "Summary!A1:B8", tableMaxRows: 8, tableMaxCols: 2, maxChars: 1200 });
try {
  await wb.render({ sheetName: "GEO_Matrix_Inventory", range: "A1:H12", scale: 1, format: "png" });
} catch {}
const out = await SpreadsheetFile.exportXlsx(wb);
await out.save(outputPath);
console.log(outputPath);
