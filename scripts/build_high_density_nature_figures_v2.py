import html
import math
import shutil
import zipfile
from pathlib import Path

import fitz
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "outputs" / "luad_signature_extraction" / "precision_oncology_deployability_audit"
WORK = ROOT / "work"
OUT = AUDIT / "publication_figures_nature_style_v2"
SVG = OUT / "svg"
PDF = OUT / "pdf"
PNG = OUT / "png"
TIFF = OUT / "tiff"
SRC = OUT / "source_data"
PKG = ROOT / "outputs" / "luad_signature_extraction" / "precision_oncology_deployability_package"
ZIP = ROOT / "outputs" / "luad_signature_extraction" / "LUAD_precision_oncology_deployability_package.zip"

W = 183.0
INK = "#222222"
MUTED = "#707070"
GRID = "#d9d9d9"
BG = "#ffffff"
COL = {
    "blue": "#2f6fb3",
    "teal": "#1b9e91",
    "green": "#4c9f70",
    "orange": "#d98c29",
    "red": "#c84d4d",
    "purple": "#8067b7",
    "grey": "#999999",
    "light_grey": "#efefef",
    "missing": "#f3f3f3",
    "overlap": "#f6e6cc",
    "coverage": "#f4cccc",
    "collapse": "#eadcf8",
    "score_missing": "#dddddd",
}
PRIMARY = ["A131", "A129", "A067", "A097", "A168", "A001", "A059", "A098", "A065", "A173"]
STRONG = {"A129", "A067", "A097", "A168", "A001", "A059"}


def read(path):
    p = AUDIT / path
    if not p.exists():
        p = WORK / path
    return pd.read_csv(p)


def esc(x):
    return html.escape(str(x), quote=True)


def tag(name, attrs=None, content=None):
    attrs = attrs or {}
    a = " ".join(f'{k}="{esc(v)}"' for k, v in attrs.items() if v is not None)
    return f"<{name} {a}/>" if content is None else f"<{name} {a}>{content}</{name}>"


def text(x, y, s, size=2.0, weight="normal", fill=INK, anchor="start", rotate=None):
    attrs = {
        "x": round(x, 3),
        "y": round(y, 3),
        "font-family": "Arial, Helvetica, DejaVu Sans, sans-serif",
        "font-size": size,
        "font-weight": weight,
        "fill": fill,
        "text-anchor": anchor,
    }
    if rotate is not None:
        attrs["transform"] = f"rotate({rotate} {round(x,3)} {round(y,3)})"
    return tag("text", attrs, esc(s))


def rect(x, y, w, h, fill, stroke=None, sw=0.22, rx=0):
    return tag("rect", {"x": round(x, 3), "y": round(y, 3), "width": round(w, 3), "height": round(h, 3), "fill": fill, "stroke": stroke, "stroke-width": sw if stroke else None, "rx": rx})


def line(x1, y1, x2, y2, stroke=INK, sw=0.35, dash=None):
    return tag("line", {"x1": round(x1, 3), "y1": round(y1, 3), "x2": round(x2, 3), "y2": round(y2, 3), "stroke": stroke, "stroke-width": sw, "stroke-dasharray": dash})


def circle(cx, cy, r, fill, stroke=None, sw=0.22, opacity=1):
    return tag("circle", {"cx": round(cx, 3), "cy": round(cy, 3), "r": round(r, 3), "fill": fill, "stroke": stroke, "stroke-width": sw if stroke else None, "opacity": opacity})


def polyline(points, stroke, sw=0.45, fill="none"):
    pts = " ".join(f"{round(x,3)},{round(y,3)}" for x, y in points)
    return tag("polyline", {"points": pts, "fill": fill, "stroke": stroke, "stroke-width": sw, "stroke-linecap": "round", "stroke-linejoin": "round"})


def panel(x, y, letter):
    return text(x, y, letter, size=3.0, weight="bold")


def save(name, width, height, parts):
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}mm" height="{height}mm" viewBox="0 0 {width} {height}">
<metadata>{{"generator":"build_high_density_nature_figures_v2.py","style":"Nature-style evidence atlas, editable SVG text"}}</metadata>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{"\n".join(parts)}
</svg>
'''
    (SVG / name).write_text(svg, encoding="utf-8")


def scale(v, lo, hi, a, b):
    if pd.isna(v):
        return (a + b) / 2
    if hi == lo:
        return (a + b) / 2
    return a + (float(v) - lo) / (hi - lo) * (b - a)


def diverge(v, lo=0.48, mid=0.58, hi=0.70):
    if pd.isna(v):
        return "#f5f5f5"
    v = max(lo, min(hi, float(v)))
    if v < mid:
        t = (v - lo) / (mid - lo)
        return rgb_interp((244, 204, 204), (245, 245, 245), t)
    t = (v - mid) / (hi - mid)
    return rgb_interp((245, 245, 245), (47, 111, 179), t)


def hr_color(hr):
    if pd.isna(hr):
        return "#f5f5f5"
    val = math.log(max(0.2, min(5, float(hr))))
    lo, hi = math.log(0.6), math.log(2.8)
    t = max(0, min(1, (val - lo) / (hi - lo)))
    return rgb_interp((210, 80, 80), (45, 125, 175), t)


def rgb_interp(a, b, t):
    return "#%02x%02x%02x" % tuple(int(round(a[i] + (b[i] - a[i]) * t)) for i in range(3))


def axis(parts, x, y, w, h, xticks=None, yticks=None, xlim=(0, 1), ylim=(0, 1), xlabel="", ylabel="", xfmt=None, yfmt=None):
    parts.append(line(x, y + h, x + w, y + h, sw=0.45))
    parts.append(line(x, y, x, y + h, sw=0.45))
    for val in yticks or []:
        yy = scale(val, ylim[1], ylim[0], y, y + h)
        parts.append(line(x - 1, yy, x, yy, sw=0.25))
        parts.append(line(x, yy, x + w, yy, stroke="#eeeeee", sw=0.20))
        parts.append(text(x - 1.8, yy + 0.7, yfmt(val) if yfmt else str(val), size=1.75, fill=MUTED, anchor="end"))
    for val in xticks or []:
        xx = scale(val, xlim[0], xlim[1], x, x + w)
        parts.append(line(xx, y + h, xx, y + h + 1, sw=0.25))
        parts.append(text(xx, y + h + 3.5, xfmt(val) if xfmt else str(val), size=1.75, fill=MUTED, anchor="middle"))
    if xlabel:
        parts.append(text(x + w / 2, y + h + 7.2, xlabel, size=2.0, anchor="middle"))
    if ylabel:
        parts.append(text(x - 9, y + h / 2, ylabel, size=2.0, anchor="middle", rotate=-90))


def status_color(status):
    return {
        "ok": None,
        "insufficient_gene_coverage": COL["coverage"],
        "skipped_training_overlap": COL["overlap"],
        "single_group_after_cutoff": COL["collapse"],
        "score_contains_missing": COL["score_missing"],
    }.get(status, COL["missing"])


def prep_ordered_signatures():
    perf = read("all_57_GEO_performance_summary.csv")
    score = read("all_57_signature_scoreability.csv")
    all_sigs = score[["Article_ID", "Signature_ID"]].drop_duplicates()
    perf2 = all_sigs.merge(perf[["Signature_ID", "median_C_index", "continuous_significant_rate", "Evidence_ladder_level"]], on="Signature_ID", how="left")
    perf2["is_primary"] = perf2["Article_ID"].isin(PRIMARY)
    perf2["primary_order"] = perf2["Article_ID"].map({a: i for i, a in enumerate(PRIMARY)}).fillna(99)
    perf2 = perf2.sort_values(["is_primary", "Evidence_ladder_level", "median_C_index", "continuous_significant_rate", "Article_ID"], ascending=[False, False, False, False, True])
    return perf2


def figure1_v2():
    parts = [panel(5, 8, "a")]
    flow = [("Full-text\narticles", 178, COL["blue"]), ("Formula-\nreproducible", 57, COL["blue"]), ("Manuscript\ntriage", 27, COL["teal"]), ("Primary\ncandidates", 10, COL["teal"])]
    x0, y0, bw, bh, gap = 17, 12, 31, 17, 7
    for i, (lab, n, c) in enumerate(flow):
        x = x0 + i * (bw + gap)
        parts.append(rect(x, y0, bw, bh, "#eef3f8", stroke="#b8c6d5", rx=1.5))
        parts.append(text(x + bw / 2, y0 + 6, lab.split("\n")[0], size=2.1, weight="bold", anchor="middle"))
        parts.append(text(x + bw / 2, y0 + 10, lab.split("\n")[1], size=2.1, weight="bold", anchor="middle"))
        parts.append(text(x + bw / 2, y0 + 15, str(n), size=3.1, weight="bold", fill=c, anchor="middle"))
        if i < len(flow) - 1:
            parts.append(line(x + bw + 1, y0 + bh / 2, x + bw + gap - 1, y0 + bh / 2, stroke=MUTED, sw=0.45))
    parts.append(text(17, 36, "Audit question: can published signatures move from formula reproducibility to deployability?", size=2.05, weight="bold"))

    parts.append(panel(5, 49, "b"))
    ladder = read("evidence_ladder_summary_counts.csv")
    sig = ladder[ladder.Entity.eq("signature")].copy()
    labels = ["L1\nformula", "L2\nscoreable", "L3\nscore\nportable", "L4\ncutoff", "L5\ndeployability\neval."]
    vals = [19, 10, 27, 0, 1]
    x, y, w, h = 21, 55, 57, 37
    axis(parts, x, y, w, h, yticks=[0, 10, 20, 30], ylim=(0, 30), ylabel="Signatures")
    barw = w / len(vals) * 0.6
    for i, v in enumerate(vals):
        xx = x + (i + 0.5) * w / len(vals) - barw / 2
        bh = v / 30 * h
        parts.append(rect(xx, y + h - bh, barw, bh, COL["blue"] if i < 2 else COL["teal"], rx=0.4))
        parts.append(text(xx + barw / 2, y + h - bh - 1, str(v), size=1.75, anchor="middle"))
        for k, l in enumerate(labels[i].split("\n")):
            parts.append(text(xx + barw / 2, y + h + 3.5 + k * 2.3, l, size=1.45, anchor="middle", fill=MUTED))

    parts.append(panel(95, 49, "c"))
    cohort = read("manuscript_geo_cohort_source_table.csv")
    cohort = cohort[cohort["Machine_ready_LUAD_like"].eq("yes")].copy()
    cohort = cohort.sort_values("GSE")
    x, y, w, h = 105, 55, 58, 37
    maxn = cohort["OS_complete_overlap_n"].max()
    for i, r in enumerate(cohort.itertuples()):
        yy = y + i * 3.8
        parts.append(text(x - 1, yy + 1.7, r.GSE, size=1.45, anchor="end"))
        parts.append(rect(x, yy, scale(r.OS_complete_overlap_n, 0, maxn, 0, w * 0.75), 2.4, COL["blue"]))
        parts.append(rect(x + w * 0.78, yy, scale(r.OS_event_n, 0, maxn, 0, w * 0.22), 2.4, COL["orange"]))
        parts.append(text(x + w + 1, yy + 1.8, str(r.GPL), size=1.4, fill=MUTED))
    parts.append(text(x, y - 3, "OS samples", size=1.7, fill=COL["blue"], weight="bold"))
    parts.append(text(x + w * 0.78, y - 3, "events", size=1.7, fill=COL["orange"], weight="bold"))

    parts.append(panel(5, 113, "d"))
    edges = read("dataset_signature_bipartite_network_edges.csv")
    freq = edges.groupby(["Dataset", "Role"]).size().reset_index(name="n").sort_values("n", ascending=False).head(16)
    x, y, maxw = 50, 114, 88
    maxv = freq["n"].max()
    for i, r in enumerate(freq.itertuples()):
        yy = y + i * 3.2
        parts.append(text(x - 2, yy + 1.6, f"{r.Dataset} ({r.Role[0]})", size=1.45, anchor="end"))
        parts.append(rect(x, yy, scale(r.n, 0, maxv, 0, maxw), 2.1, COL["purple"] if r.Role == "training" else COL["green"]))
        parts.append(text(x + scale(r.n, 0, maxv, 0, maxw) + 1, yy + 1.6, str(r.n), size=1.45))
    parts.append(text(50, 170, "t: training; v: validation", size=1.55, fill=MUTED))
    save("Figure1_study_design_evidence_ladder_v2.svg", W, 176, parts)


def figure2_v2():
    rq = read("reporting_quality_score_all_articles.csv")
    items = ["gene_list", "coefficients", "formula", "numeric_cutoff", "cutoff_derivation", "training_cohort", "validation_cohort", "endpoint", "preprocessing", "code_data_availability"]
    item_labels = ["Gene", "Coeff.", "Formula", "Num.\ncutoff", "Cutoff\nderiv.", "Training", "Validation", "Endpoint", "Preproc.", "Code/\ndata"]
    rq = rq.sort_values(["Evidence_ladder_level", "Reporting_quality_score_0_10", "Year"], ascending=[False, False, True]).reset_index(drop=True)
    parts = [panel(5, 8, "a")]
    x0, y0, cw, ch = 35, 13, 4.7, 1.15
    for j, lab in enumerate(item_labels):
        for k, l in enumerate(lab.split("\n")):
            parts.append(text(x0 + j * cw + cw / 2, y0 - 5 + k * 2.3, l, size=1.25, fill=MUTED, anchor="middle"))
    for i, r in rq.iterrows():
        yy = y0 + i * ch
        level = int(r["Evidence_ladder_level"]) if pd.notna(r["Evidence_ladder_level"]) else 0
        parts.append(rect(8, yy, 5, ch - 0.1, [COL["light_grey"], "#dfeaf7", "#dfeaf7", "#dff1ef", "#f6e7d5", "#eadff4"][level]))
        for j, item in enumerate(items):
            fill = COL["green"] if r[item] == 1 else "#eeeeee"
            parts.append(rect(x0 + j * cw, yy, cw - 0.08, ch - 0.08, fill))
    parts.append(text(8, y0 - 5, "L", size=1.45, fill=MUTED))
    parts.append(text(8, y0 + len(rq) * ch + 4, "178 articles sorted by evidence level and reporting score", size=1.75, fill=MUTED))

    parts.append(panel(93, 8, "b"))
    x, y, w, h = 112, 15, 42, 55
    # histogram score by reproducibility
    bins = list(range(0, 11))
    for flag, color, xoff, lab in [(False, COL["grey"], -0.8, "Not reproducible"), (True, COL["teal"], 0.8, "Formula-reproducible")]:
        sub = rq[rq["Reproducible"].astype(str).str.lower().eq("yes") == flag]["Reporting_quality_score_0_10"]
        counts = sub.value_counts().reindex(bins, fill_value=0)
        maxc = max(1, counts.max())
        bw = w / len(bins) * 0.38
        for i, b in enumerate(bins):
            xx = x + (i + 0.5) * w / len(bins) + xoff - bw / 2
            hh = counts.loc[b] / maxc * h
            parts.append(rect(xx, y + h - hh, bw, hh, color, rx=0.25, stroke=None))
        parts.append(text(x + 4, y + (0 if flag else 5), lab, size=1.7, fill=color, weight="bold"))
    axis(parts, x, y, w, h, xticks=[0, 2, 4, 6, 8, 10], yticks=[], xlim=(0, 10), ylim=(0, 1), xlabel="Reporting score")

    parts.append(panel(93, 86, "c"))
    year = read("reporting_score_by_year.csv").dropna(subset=["Year"])
    x, y, w, h = 110, 90, 55, 34
    years = year["Year"].astype(int)
    axis(parts, x, y, w, h, xticks=[years.min(), 2020, years.max()], yticks=[0, 2, 4, 6, 8, 10], xlim=(years.min(), years.max()), ylim=(0, 10), xlabel="Year", ylabel="Mean score")
    pts = [(scale(r.Year, years.min(), years.max(), x, x + w), scale(r.Mean_reporting_score, 10, 0, y, y + h)) for r in year.itertuples()]
    parts.append(polyline(pts, COL["blue"], sw=0.5))
    for r, (xx, yy) in zip(year.itertuples(), pts):
        parts.append(circle(xx, yy, 0.65 + 0.12 * r.Formula_reproducible_count, COL["orange"], opacity=0.75))
    save("Figure2_reporting_failure_anatomy_v2.svg", W, 132, parts)


def figure3_v2():
    pair = read("geo_external_validation_results.csv")
    cohort = read("manuscript_geo_cohort_source_table.csv")
    order_df = prep_ordered_signatures()
    sigs = order_df["Signature_ID"].tolist()
    arts = order_df["Article_ID"].tolist()
    gses = ["GSE13213", "GSE30219", "GSE31210", "GSE3141", "GSE37745", "GSE42127", "GSE50081", "GSE68465", "GSE72094"]
    lookup = {(r.Signature_ID, r.GSE): r for r in pair.itertuples()}
    perf = read("all_57_GEO_performance_summary.csv").set_index("Signature_ID")
    endpoint = read("endpoint_mismatch_with_GEO_OS_validation.csv").set_index("Signature_ID")
    readiness = read("cutoff_reporting_audit_57.csv").set_index("Signature_ID")
    parts = [panel(4, 7, "a")]
    x0, y0, cw, ch = 28, 20, 8.0, 2.18
    for j, gse in enumerate(gses):
        parts.append(text(x0 + j * cw + cw / 2, y0 - 5, gse.replace("GSE", ""), size=1.35, anchor="middle", rotate=-45, fill=MUTED))
    # top platform strip
    platform_color = {"Affymetrix microarray": COL["blue"], "Illumina/other microarray": COL["teal"], "Other/unknown microarray": COL["grey"]}
    for j, gse in enumerate(gses):
        crow = cohort[cohort.GSE.eq(gse)].iloc[0]
        pclass = "Affymetrix microarray" if str(crow.GPL) in {"GPL570", "GPL96"} else "Illumina/other microarray"
        parts.append(rect(x0 + j * cw, y0 - 2.5, cw - 0.2, 1.8, platform_color[pclass]))
    for i, (sig, art) in enumerate(zip(sigs, arts)):
        yy = y0 + i * ch
        is_primary = art in PRIMARY
        parts.append(text(x0 - 2, yy + 1.55, art, size=1.35, anchor="end", weight="bold" if is_primary else "normal", fill=INK if is_primary else MUTED))
        for j, gse in enumerate(gses):
            r = lookup.get((sig, gse))
            xx = x0 + j * cw
            if r is None:
                fill = COL["missing"]
                parts.append(rect(xx, yy, cw - 0.18, ch - 0.12, fill))
                continue
            if r.status == "ok":
                fill = hr_color(r.Cox_continuous_HR_per_SD)
                parts.append(rect(xx, yy, cw - 0.18, ch - 0.12, fill, stroke="#ffffff", sw=0.08))
                if pd.notna(r.Cox_continuous_p) and r.Cox_continuous_p < 0.05:
                    parts.append(circle(xx + cw / 2, yy + ch / 2, 0.45, "#000000", opacity=0.72))
                if pd.notna(r.Cox_continuous_p) and r.Cox_continuous_p < 0.05 and pd.notna(r.Cox_high_vs_low_p) and r.Cox_high_vs_low_p < 0.05:
                    parts.append(circle(xx + cw / 2 + 1.1, yy + ch / 2, 0.38, "#ffffff", stroke="#000000", sw=0.18))
            else:
                fill = status_color(r.status)
                parts.append(rect(xx, yy, cw - 0.18, ch - 0.12, fill, stroke="#ffffff", sw=0.08))
                symbol = {"insufficient_gene_coverage": "×", "skipped_training_overlap": "T", "single_group_after_cutoff": "C", "score_contains_missing": "M"}.get(r.status, "")
                if symbol:
                    parts.append(text(xx + cw / 2, yy + 1.55, symbol, size=1.45, weight="bold", anchor="middle", fill=INK))
    # Side annotations
    annx = x0 + len(gses) * cw + 4
    headers = [("k", 8), ("C", 8), ("L", 7), ("E", 7), ("N", 7)]
    for idx, (h, width) in enumerate(headers):
        parts.append(text(annx + idx * 8, y0 - 5, h, size=1.4, weight="bold", anchor="middle", fill=MUTED))
    for i, sig in enumerate(sigs):
        yy = y0 + i * ch + 1.5
        if sig in perf.index:
            p = perf.loc[sig]
            vals = [p.evaluable_cohort_count, p.median_C_index, p.Evidence_ladder_level]
        else:
            vals = [0, np.nan, 1]
        ep_mis = endpoint.loc[sig, "Endpoint_mismatch_with_GEO_OS"] if sig in endpoint.index else False
        num_cut = bool(readiness.loc[sig, "reported_numeric_cutoff"]) if sig in readiness.index else False
        parts.append(text(annx, yy, str(int(vals[0])), size=1.25, anchor="middle"))
        parts.append(text(annx + 8, yy, "" if pd.isna(vals[1]) else f"{vals[1]:.2f}", size=1.25, anchor="middle"))
        parts.append(text(annx + 16, yy, f"L{int(vals[2])}", size=1.25, anchor="middle"))
        parts.append(text(annx + 24, yy, "Y" if ep_mis else "N", size=1.25, anchor="middle", fill=COL["orange"] if ep_mis else MUTED))
        parts.append(text(annx + 32, yy, "Y" if num_cut else "N", size=1.25, anchor="middle", fill=COL["blue"] if num_cut else MUTED))
    parts.append(text(annx - 1, y0 + len(sigs) * ch + 4, "k: evaluable GEOs; C: median C-index; L: ladder; E: endpoint mismatch; N: numeric cutoff", size=1.4, fill=MUTED))
    # Legend
    lx, ly = 18, 149
    parts.append(text(lx, ly, "Continuous HR/SD", size=1.7, weight="bold"))
    for i, (lab, col) in enumerate([("<1", hr_color(0.7)), ("1", hr_color(1.0)), (">1", hr_color(2.2))]):
        parts.append(rect(lx + i * 10, ly + 3, 8, 3, col))
        parts.append(text(lx + i * 10 + 4, ly + 8, lab, size=1.25, anchor="middle"))
    for i, (lab, col) in enumerate([("gene coverage", COL["coverage"]), ("training overlap", COL["overlap"]), ("cutoff collapse", COL["collapse"]), ("score missing", COL["score_missing"])]):
        parts.append(rect(58 + i * 25, ly + 3, 4, 3, col))
        parts.append(text(63 + i * 25, ly + 5.5, lab, size=1.25))
    parts.append(circle(158, ly + 4.5, 0.45, "#000000"))
    parts.append(text(160, ly + 5.2, "continuous p<0.05", size=1.25))
    save("Figure3_all57_external_transportability_atlas_v2.svg", W, 162, parts)


def figure4_v2():
    drift = read("cutoff_percentile_drift_long.csv")
    prop = read("reported_cutoff_high_risk_proportion_long.csv")
    audit = read("cutoff_reporting_audit_57.csv")
    parts = [panel(5, 8, "a")]
    # cutoff cascade
    steps = [
        ("Formula-reproducible", 57, COL["blue"]),
        ("Numeric cutoff", int(audit["reported_numeric_cutoff"].sum()), COL["blue"]),
        ("GEO applicable", int(audit["reported_cutoff_applicable_to_GEO"].sum()), COL["teal"]),
        ("No collapse", int(((audit["reported_cutoff_applicable_to_GEO"]) & (~audit["single_group_collapse_observed"])).sum()), COL["green"]),
        ("Strict primary", 1, COL["green"]),
    ]
    x0, y0, maxw, bh, gap = 52, 13, 85, 7.2, 4
    maxv = steps[0][1]
    for i, (lab, val, col) in enumerate(steps):
        yy = y0 + i * (bh + gap)
        parts.append(text(x0 - 2, yy + 4.8, lab, size=1.9, anchor="end"))
        parts.append(rect(x0, yy, scale(val, 0, maxv, 0, maxw), bh, col))
        parts.append(text(x0 + scale(val, 0, maxv, 0, maxw) + 1.5, yy + 4.8, str(val), size=1.9, weight="bold"))
    parts.append(panel(5, 73, "b"))
    x, y, w, h = 33, 77, 54, 38
    axis(parts, x, y, w, h, xticks=[0, 0.25, 0.5, 0.75, 1], xlim=(0, 1), ylim=(0, 1), xlabel="High-risk proportion")
    parts.append(line(scale(0.5, 0, 1, x, x + w), y, scale(0.5, 0, 1, x, x + w), y + h, stroke=COL["blue"], dash="1,1"))
    jitter = {}
    for r in prop.itertuples():
        hp = r.reported_cutoff_high_risk_proportion
        if pd.isna(hp):
            continue
        key = round(hp, 2)
        jitter[key] = jitter.get(key, 0) + 1
        xx = scale(hp, 0, 1, x, x + w)
        yy = y + 4 + (jitter[key] % 10) * 3.0
        color = COL["red"] if r.Article_ID == "A097" else COL["green"] if r.Status == "ok" else COL["orange"]
        parts.append(circle(xx, yy, 0.9, color, opacity=0.82))
    parts.append(text(31, 122, "Each dot is a reported-cutoff GEO test; blue line marks a balanced 50% split.", size=1.45, fill=MUTED))

    parts.append(panel(101, 73, "c"))
    a097 = drift[drift.Article_ID.eq("A097")]
    x, y, w, h = 125, 76, 40, 32
    axis(parts, x, y, w, h, xticks=[0, 0.5, 1], xlim=(0, 1), ylim=(0, 1), xlabel="High-risk proportion\nA097 cutoff 1.66")
    for i, r in enumerate(a097.itertuples()):
        yy = y + 3 + i * 4.5
        hp = r.High_n / r.N
        parts.append(text(x - 2, yy + 0.7, r.GSE, size=1.45, anchor="end"))
        parts.append(line(x, yy, scale(hp, 0, 1, x, x + w), yy, stroke="#cccccc"))
        parts.append(circle(scale(hp, 0, 1, x, x + w), yy, 1, COL["red"]))
    parts.append(text(102, 117, "A097: cutoff exists, but fixed-threshold GEO", size=1.55, fill=COL["red"], weight="bold"))
    parts.append(text(102, 121, "application collapses all tested cohorts.", size=1.55, fill=COL["red"], weight="bold"))

    parts.append(panel(5, 137, "d"))
    # A131 strict contrast
    pair = read("geo_external_validation_results.csv")
    a131 = pair[(pair.Article_ID.eq("A131")) & (pair.status.eq("ok"))].copy()
    x, y, w, h = 35, 139, 95, 26
    axis(parts, x, y, w, h, xticks=[0.5, 1, 2, 4], xlim=(math.log(0.5), math.log(4)), ylim=(0, 1), xlabel="A131 reported-cutoff HR (log scale)", xfmt=lambda v: f"{v:g}")
    for i, r in enumerate(a131.itertuples()):
        yy = y + 2 + i * 2.8
        hr = r.Cox_high_vs_low_HR
        lo = r.Cox_high_vs_low_CI95_low
        hi = r.Cox_high_vs_low_CI95_high
        parts.append(text(x - 2, yy + 0.6, r.GSE, size=1.35, anchor="end"))
        if pd.notna(lo) and pd.notna(hi):
            parts.append(line(scale(math.log(lo), math.log(0.5), math.log(4), x, x + w), yy, scale(math.log(hi), math.log(0.5), math.log(4), x, x + w), yy, stroke=MUTED, sw=0.35))
        parts.append(circle(scale(math.log(hr), math.log(0.5), math.log(4), x, x + w), yy, 0.9, COL["blue"]))
    parts.append(line(scale(math.log(1), math.log(0.5), math.log(4), x, x + w), y, scale(math.log(1), math.log(0.5), math.log(4), x, x + w), y + h, stroke="#999999", dash="1,1"))
    save("Figure4_cutoff_transportability_failure_anatomy_v2.svg", W, 174, parts)


def figure5_v2():
    pair = read("geo_external_validation_results.csv")
    inc = read("incremental_cindex_by_signature_cohort.csv")
    lrt = read("likelihood_ratio_test_clinical_vs_clinical_plus_risk.csv")
    auc = read("time_dependent_auc_primary_candidates.csv")
    brier = read("brier_score_primary_candidates.csv")
    parts = [panel(5, 8, "a")]
    main = pair[(pair.Article_ID.isin(PRIMARY)) & (pair.status.eq("ok"))].copy()
    x0, y0, cw, ch = 25, 16, 11.0, 5.0
    gses = sorted(main.GSE.unique())
    for j, gse in enumerate(gses):
        parts.append(text(x0 + j * cw + cw / 2, y0 - 4, gse.replace("GSE", ""), size=1.35, rotate=-45, anchor="middle", fill=MUTED))
    for i, art in enumerate(PRIMARY):
        parts.append(text(x0 - 2, y0 + i * ch + 3.4, art, size=1.55, anchor="end", weight="bold" if art == "A131" else "normal"))
        sub = main[main.Article_ID.eq(art)].set_index("GSE")
        for j, gse in enumerate(gses):
            xx, yy = x0 + j * cw, y0 + i * ch
            if gse not in sub.index:
                parts.append(rect(xx, yy, cw - 0.2, ch - 0.2, COL["missing"]))
                continue
            r = sub.loc[gse]
            parts.append(rect(xx, yy, cw - 0.2, ch - 0.2, hr_color(r.Cox_continuous_HR_per_SD), stroke="#ffffff", sw=0.1))
            if r.Cox_continuous_p < 0.05:
                parts.append(circle(xx + cw / 2, yy + ch / 2, 0.55, "#000000", opacity=0.75))

    parts.append(panel(5, 83, "b"))
    x, y, w, h = 38, 84, 52, 42
    axis(parts, x, y, w, h, xticks=[0, 0.05, 0.10, 0.15], xlim=(0, 0.15), ylim=(0, 1), xlabel="ΔC-index")
    for i, art in enumerate(PRIMARY):
        vals = inc[inc.Article_ID.eq(art)]["Delta_C_index"].dropna().values
        if len(vals) == 0:
            continue
        yy = y + 2 + i * 3.7
        parts.append(text(x - 2, yy + 0.7, art, size=1.55, anchor="end"))
        for v in vals:
            parts.append(circle(scale(v, 0, 0.15, x, x + w), yy + np.random.default_rng(abs(hash((art, round(v, 3)))) % 9999).uniform(-0.7, 0.7), 0.45, COL["teal"], opacity=0.45))
        med = np.median(vals)
        parts.append(circle(scale(med, 0, 0.15, x, x + w), yy, 1.0, COL["blue"] if art == "A131" else COL["teal"]))

    parts.append(panel(103, 83, "c"))
    x, y, w, h = 126, 84, 38, 42
    axis(parts, x, y, w, h, xticks=[0, 0.5, 1.0], xlim=(0, 1), ylim=(0, 1), xlabel="LRT p<0.05 rate")
    for i, art in enumerate(PRIMARY):
        rate = (lrt[lrt.Article_ID.eq(art)]["LR_test_p"] < 0.05).mean()
        yy = y + 2 + i * 3.7
        parts.append(text(x - 2, yy + 0.7, art, size=1.55, anchor="end"))
        parts.append(rect(x, yy - 0.8, scale(rate, 0, 1, 0, w), 1.8, COL["blue"] if art == "A131" else COL["teal"]))
        parts.append(text(x + scale(rate, 0, 1, 0, w) + 1, yy + 0.7, f"{rate:.2f}", size=1.35))

    parts.append(panel(5, 143, "d"))
    metrics = []
    for art in PRIMARY:
        metrics.append({
            "Article_ID": art,
            "AUC3": auc[(auc.Article_ID.eq(art)) & (auc.time_years.eq(3.0))]["dynamic_AUC_complete_case"].median(),
            "AUC5": auc[(auc.Article_ID.eq(art)) & (auc.time_years.eq(5.0))]["dynamic_AUC_complete_case"].median(),
            "Brier3": brier[(brier.Article_ID.eq(art)) & (brier.time_years.eq(3.0))]["IPCW_Brier_score_external_recalibrated"].median(),
            "Brier5": brier[(brier.Article_ID.eq(art)) & (brier.time_years.eq(5.0))]["IPCW_Brier_score_external_recalibrated"].median(),
        })
    met = pd.DataFrame(metrics)
    x0, y0, cw, ch = 40, 144, 18, 4.6
    cols = ["AUC3", "AUC5", "Brier3", "Brier5"]
    for j, lab in enumerate(["AUC 3y", "AUC 5y", "Brier 3y", "Brier 5y"]):
        parts.append(text(x0 + j * cw + cw / 2, y0 - 3, lab, size=1.55, fill=MUTED, anchor="middle"))
    for i, r in met.iterrows():
        parts.append(text(x0 - 2, y0 + i * ch + 3, r.Article_ID, size=1.45, anchor="end"))
        for j, c in enumerate(cols):
            v = r[c]
            if c.startswith("Brier"):
                norm = 1 - max(0, min(1, (v - 0.15) / 0.13))
            else:
                norm = max(0, min(1, (v - 0.55) / 0.18))
            parts.append(rect(x0 + j * cw, y0 + i * ch, cw - 0.2, ch - 0.2, rgb_interp((241, 241, 241), (47, 111, 179), norm), stroke="#ffffff", sw=0.1))
            parts.append(text(x0 + j * cw + cw / 2, y0 + i * ch + 3, f"{v:.2f}", size=1.35, anchor="middle", fill="#ffffff" if norm > 0.55 else INK))
    save("Figure5_primary_candidate_clinical_prediction_v2.svg", W, 195, parts)


def figure6_v2():
    heat = read("clinical_deployability_heatmap_data.csv")
    score = read("clinical_deployability_score_summary.csv")
    parts = [panel(5, 8, "a")]
    domains = list(heat.Deployability_domain.drop_duplicates())
    groups = {
        "Reproducibility": ["complete_formula", "gene_availability", "endpoint_compatibility_with_GEO_OS"],
        "Threshold": ["reported_numeric_cutoff", "strict_cutoff_testability", "no_single_group_collapse", "clinical_threshold_justification"],
        "External performance": ["continuous_GEO_performance", "clinical_adjusted_signal", "delta_C_index_positive", "time_AUC_available", "Brier_available", "external_recalibration_available", "TCGA_formula_direction_sanity", "leave_one_cohort_out_direction_retained"],
        "Implementation": ["original_code_or_calculator", "assay_standardization", "prospective_validation"],
    }
    ordered = [d for g in groups.values() for d in g]
    short = {
        "complete_formula": "Formula", "gene_availability": "Genes", "endpoint_compatibility_with_GEO_OS": "Endpoint",
        "reported_numeric_cutoff": "Cutoff", "strict_cutoff_testability": "Strict", "no_single_group_collapse": "No\ncollapse", "clinical_threshold_justification": "Threshold\nrationale",
        "continuous_GEO_performance": "GEO", "clinical_adjusted_signal": "Adj.", "delta_C_index_positive": "ΔC", "time_AUC_available": "AUC", "Brier_available": "Brier", "external_recalibration_available": "Recal.", "TCGA_formula_direction_sanity": "TCGA", "leave_one_cohort_out_direction_retained": "LOCO",
        "original_code_or_calculator": "Code", "assay_standardization": "Assay", "prospective_validation": "Prospective",
    }
    val = {(r.Article_ID, r.Deployability_domain): r.Heatmap_value for r in heat.itertuples()}
    x0, y0, cw, ch = 30, 26, 7.3, 5.7
    cursor = 0
    for gname, ds in groups.items():
        gx = x0 + cursor * cw
        gw = len(ds) * cw
        parts.append(rect(gx, y0 - 12, gw - 0.4, 4, {"Reproducibility": "#e8eef7", "Threshold": "#f6e7d5", "External performance": "#dff1ef", "Implementation": "#eeeeee"}[gname]))
        parts.append(text(gx + gw / 2, y0 - 9.2, gname, size=1.45, anchor="middle", weight="bold"))
        cursor += len(ds)
    for j, d in enumerate(ordered):
        for k, l in enumerate(short[d].split("\n")):
            parts.append(text(x0 + j * cw + cw / 2, y0 - 2.5 + k * 2.2, l, size=1.25, rotate=-35, anchor="middle", fill=MUTED))
    for i, art in enumerate(PRIMARY):
        yy = y0 + i * ch
        parts.append(text(x0 - 2, yy + 3.8, art, size=1.6, anchor="end", weight="bold" if art == "A131" else "normal"))
        for j, d in enumerate(ordered):
            v = val.get((art, d), 0)
            fill = COL["green"] if v >= 1 else "#d8b545" if v >= 0.5 else "#9ab7db" if v >= 0.25 else "#d9d9d9"
            parts.append(rect(x0 + j * cw, yy, cw - 0.25, ch - 0.25, fill, stroke="#ffffff", sw=0.08))
    parts.append(text(30, 88, "Green: yes; gold: partial; blue: manual check; grey: no.", size=1.55, fill=MUTED))

    parts.append(panel(5, 104, "b"))
    x, y, w, h = 38, 106, 52, 36
    axis(parts, x, y, w, h, xticks=[0, 6, 12, 18], xlim=(0, 18), ylim=(0, 1), xlabel="Readiness score")
    score["order"] = score.Article_ID.map({a: i for i, a in enumerate(PRIMARY)})
    for i, r in enumerate(score.sort_values("order").itertuples()):
        yy = y + 2 + i * 3.1
        parts.append(text(x - 2, yy + 0.7, r.Article_ID, size=1.45, anchor="end"))
        parts.append(rect(x, yy - 0.8, scale(r.Deployability_score_0_18, 0, 18, 0, w), 1.8, COL["blue"] if r.Article_ID == "A131" else COL["teal"]))
        parts.append(text(x + scale(r.Deployability_score_0_18, 0, 18, 0, w) + 1, yy + 0.7, f"{r.Deployability_score_0_18:.1f}", size=1.4))

    parts.append(panel(105, 104, "c"))
    x, y = 107, 106
    parts.append(rect(x, y, 58, 30, "#f5f5f5", stroke="#c8c8c8", rx=1.2))
    lines = [
        "Portable prognostic signal",
        "≠ cutoff transportability",
        "≠ clinical deployability",
    ]
    for i, l in enumerate(lines):
        parts.append(text(x + 29, y + 9 + i * 7, l, size=2.25, weight="bold" if i == 0 else "normal", anchor="middle", fill=COL["blue"] if i == 0 else COL["red"] if i == 1 else INK))
    parts.append(text(x + 29, y + 38, "Main message for precision-oncology readiness", size=1.55, anchor="middle", fill=MUTED))
    save("Figure6_clinical_deployability_gap_v2.svg", W, 150, parts)


def convert():
    zoom = 600 / 72
    mat = fitz.Matrix(zoom, zoom)
    for svg in SVG.glob("*.svg"):
        doc = fitz.open(str(svg))
        (PDF / f"{svg.stem}.pdf").write_bytes(doc.convert_to_pdf())
        pix = doc[0].get_pixmap(matrix=mat, alpha=False)
        png = PNG / f"{svg.stem}.png"
        pix.save(str(png))
        img = Image.open(png)
        img.save(TIFF / f"{svg.stem}.tiff", compression="tiff_lzw", dpi=(600, 600))
        doc.close()


def source_data_index():
    mapping = {
        "Figure1": ["evidence_ladder_summary_counts.csv", "manuscript_geo_cohort_source_table.csv", "dataset_signature_bipartite_network_edges.csv"],
        "Figure2": ["reporting_quality_score_all_articles.csv", "reporting_score_by_year.csv"],
        "Figure3": ["geo_external_validation_results.csv", "all_57_GEO_performance_summary.csv", "endpoint_mismatch_with_GEO_OS_validation.csv", "cutoff_reporting_audit_57.csv"],
        "Figure4": ["cutoff_reporting_audit_57.csv", "cutoff_percentile_drift_long.csv", "reported_cutoff_high_risk_proportion_long.csv", "geo_external_validation_results.csv"],
        "Figure5": ["geo_external_validation_results.csv", "incremental_cindex_by_signature_cohort.csv", "likelihood_ratio_test_clinical_vs_clinical_plus_risk.csv", "time_dependent_auc_primary_candidates.csv", "brier_score_primary_candidates.csv"],
        "Figure6": ["clinical_deployability_heatmap_data.csv", "clinical_deployability_score_summary.csv"],
    }
    rows = []
    for fig, files in mapping.items():
        d = SRC / fig
        d.mkdir(parents=True, exist_ok=True)
        for f in files:
            src = AUDIT / f if (AUDIT / f).exists() else WORK / f
            if src.exists():
                shutil.copy2(src, d / f)
            rows.append({"Figure": fig, "Source_file": f})
    pd.DataFrame(rows).to_csv(SRC / "figure_source_data_index.csv", index=False, encoding="utf-8")


def validate():
    files = sorted(SVG.glob("*.svg"))
    if len(files) != 6:
        raise RuntimeError(f"Expected 6 SVG files, found {len(files)}")
    rows = []
    for svg in files:
        stem = svg.stem
        t = svg.read_text(encoding="utf-8")
        if t.count("<text") < 20:
            raise RuntimeError(f"Too few editable text nodes: {svg}")
        row = {"Figure": stem, "svg_bytes": svg.stat().st_size, "text_nodes": t.count("<text")}
        for ext, directory in [("pdf", PDF), ("png", PNG), ("tiff", TIFF)]:
            p = directory / f"{stem}.{ext}"
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"Missing {p}")
            row[f"{ext}_bytes"] = p.stat().st_size
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "publication_figure_v2_manifest.csv", index=False, encoding="utf-8")


def contact_sheet():
    imgs = []
    for p in sorted(PNG.glob("*.png")):
        im = Image.open(p).convert("RGB")
        im.thumbnail((900, 710), Image.LANCZOS)
        canvas = Image.new("RGB", (940, 770), "white")
        canvas.paste(im, ((940 - im.width) // 2, 45))
        # use default bitmap font for label
        from PIL import ImageDraw
        d = ImageDraw.Draw(canvas)
        d.text((18, 15), p.stem, fill=(25, 25, 25))
        imgs.append(canvas)
    out = Image.new("RGB", (1880, 2310), "white")
    for i, im in enumerate(imgs):
        out.paste(im, ((i % 2) * 940, (i // 2) * 770))
    out.save(OUT / "publication_figures_v2_contact_sheet.png", dpi=(180, 180))


def refresh_package():
    target = PKG / "figures_nature_style_v2"
    if target.exists():
        shutil.rmtree(target)
    shutil.copytree(OUT, target)
    if ZIP.exists():
        ZIP.unlink()
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in PKG.rglob("*"):
            if p.is_file():
                zf.write(p, p.relative_to(PKG.parent))


def main():
    for d in [SVG, PDF, PNG, TIFF, SRC]:
        d.mkdir(parents=True, exist_ok=True)
    source_data_index()
    figure1_v2()
    figure2_v2()
    figure3_v2()
    figure4_v2()
    figure5_v2()
    figure6_v2()
    convert()
    validate()
    contact_sheet()
    refresh_package()
    print(f"High-density Nature-style v2 figures: {OUT}")
    print(f"Package refreshed: {ZIP}")


if __name__ == "__main__":
    main()
