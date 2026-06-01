import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "@oai/artifact-tool";

const root = process.cwd();
const outputDir = path.join(root, "outputs", "luad_signature_extraction");
const outputPath = path.join(outputDir, "LUAD_GEO_clinical_survival_inventory.xlsx");

async function readCsv(file) {
  const text = await fs.readFile(file, "utf8");
  const lines = text.replace(/^\uFEFF/, "").trimEnd().split(/\r?\n/);
  if (!lines.length) return [];
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

function sanitize(v) {
  if (v === null || v === undefined) return "";
  const s = String(v);
  return s.length > 32000 ? s.slice(0, 31950) + " [truncated]" : s;
}

function widthFor(field) {
  if (/file|source|sample_type|histology|stage|transform/i.test(field)) return 36;
  if (/GSE|GPL|event|count|ready|time/i.test(field)) return 16;
  return 22;
}

async function addSheet(wb, name, rows, fields) {
  const sheet = wb.worksheets.add(name);
  const data = [fields, ...rows.map((r) => fields.map((f) => sanitize(r[f])))];
  const endCol = colLetter(fields.length);
  sheet.getRange(`A1:${endCol}${data.length}`).values = data;
  sheet.freezePanes.freezeRows(1);
  sheet.getRange(`A1:${endCol}1`).format = { fill: "#7F6000", font: { bold: true, color: "#FFFFFF" }, wrapText: true };
  sheet.getRange(`A1:${endCol}${data.length}`).format = { wrapText: true };
  fields.forEach((f, i) => (sheet.getRange(`${colLetter(i + 1)}:${colLetter(i + 1)}`).format.columnWidth = widthFor(f)));
  try {
    sheet.tables.add(`A1:${endCol}${data.length}`, true, `${name.replace(/[^A-Za-z0-9]/g, "")}Table`);
  } catch {}
}

await fs.mkdir(outputDir, { recursive: true });
const readiness = await readCsv(path.join(root, "work", "geo_expression_clinical_readiness.csv"));
const summary = await readCsv(path.join(root, "work", "geo_clinical_summary.csv"));
const samples = await readCsv(path.join(root, "work", "geo_clinical_samples.csv"));

const summaryRows = [
  { Metric: "GEO_series_with_gene_expression", Value: readiness.filter((r) => r.Matrix_status === "ok").length },
  { Metric: "Machine_ready_LUAD_like_OS", Value: readiness.filter((r) => r.Machine_ready_LUAD_like === "yes").length },
  { Metric: "Total_LUAD_like_OS_complete_samples", Value: readiness.reduce((a, r) => a + Number(r.LUAD_like_OS_complete_n || 0), 0) },
  { Metric: "Total_LUAD_like_OS_events", Value: readiness.reduce((a, r) => a + Number(r.LUAD_like_OS_event_n || 0), 0) },
  { Metric: "Note", Value: "GSE50081 survival_time is converted from years to months based on field values and dataset convention; verify before final manuscript." },
];

const wb = Workbook.create();
await addSheet(wb, "Summary", summaryRows, ["Metric", "Value"]);
await addSheet(wb, "Readiness_By_Cohort", readiness, [
  "GSE",
  "Matrix_status",
  "GPL",
  "Expression_sample_count",
  "Clinical_sample_count",
  "Expression_clinical_overlap_n",
  "OS_complete_overlap_n",
  "OS_event_overlap_n",
  "LUAD_like_OS_complete_n",
  "LUAD_like_OS_event_n",
  "Machine_ready_all_histologies",
  "Machine_ready_LUAD_like",
  "Scale_guess",
  "Scoring_transform",
  "Clinical_OS_time_source",
  "Clinical_OS_event_source",
  "Scoring_expression_file",
]);
await addSheet(wb, "Clinical_Summary_Raw", summary, [
  "GSE",
  "sample_count",
  "OS_complete_n",
  "OS_event_n",
  "RFS_complete_n",
  "RFS_event_n",
  "OS_time_source",
  "OS_event_source",
  "RFS_time_source",
  "RFS_event_source",
  "machine_survival_ready",
  "raw_metadata_columns",
]);
await addSheet(wb, "Clinical_Samples", samples, [
  "GSE",
  "geo_accession",
  "OS_time_months",
  "OS_event",
  "RFS_time_months",
  "RFS_event",
  "stage",
  "histology",
  "sample_type",
  "sex",
  "age",
  "OS_time_source",
  "OS_event_source",
  "RFS_time_source",
  "RFS_event_source",
]);

await wb.inspect({ kind: "table", range: "Readiness_By_Cohort!A1:L12", tableMaxRows: 12, tableMaxCols: 12, maxChars: 1800 });
try {
  await wb.render({ sheetName: "Readiness_By_Cohort", range: "A1:L12", scale: 1, format: "png" });
} catch {}
const out = await SpreadsheetFile.exportXlsx(wb);
await out.save(outputPath);
console.log(outputPath);
