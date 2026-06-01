import csv
import html
import json
import math
import shutil
import subprocess
import zipfile
from pathlib import Path

import fitz
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "outputs" / "luad_signature_extraction" / "precision_oncology_deployability_audit"
OUT = AUDIT / "publication_figures_nature_style"
SVG_DIR = OUT / "svg"
PNG_DIR = OUT / "png"
TIFF_DIR = OUT / "tiff"
PDF_DIR = OUT / "pdf"
HTML_DIR = OUT / "html"
SRC_DIR = OUT / "source_data"
PKG = ROOT / "outputs" / "luad_signature_extraction" / "precision_oncology_deployability_package"
ZIP = ROOT / "outputs" / "luad_signature_extraction" / "LUAD_precision_oncology_deployability_package.zip"

W = 183.0
INK = "#222222"
MUTED = "#6f6f6f"
GRID = "#d8d8d8"
LIGHT = "#f3f3f3"
PALETTE = {
    "blue": "#3b6fb6",
    "teal": "#1b9e91",
    "green": "#4c9f70",
    "orange": "#d98c29",
    "red": "#c84d4d",
    "purple": "#7c65b8",
    "gold": "#b69b30",
    "grey": "#8d8d8d",
    "light_blue": "#dce8f5",
    "light_teal": "#dff1ef",
    "light_orange": "#f6e7d5",
    "light_red": "#f3dada",
}

PRIMARY = ["A131", "A129", "A067", "A097", "A168", "A001", "A059", "A098", "A065", "A173"]
STRONG = {"A129", "A067", "A097", "A168", "A001", "A059"}
LOW = {"A065", "A173"}


def read(name):
    return pd.read_csv(AUDIT / name)


def esc(x):
    return html.escape(str(x), quote=True)


def fmt(x, digits=2):
    if pd.isna(x):
        return "NA"
    return f"{float(x):.{digits}f}"


def tag(name, attrs=None, content=None):
    attrs = attrs or {}
    a = " ".join(f'{k}="{esc(v)}"' for k, v in attrs.items() if v is not None)
    if content is None:
        return f"<{name} {a}/>"
    return f"<{name} {a}>{content}</{name}>"


def text(x, y, s, size=2.25, weight="normal", fill=INK, anchor="start", rotate=None):
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


def rect(x, y, w, h, fill, stroke=None, sw=0.25, rx=0):
    return tag(
        "rect",
        {
            "x": round(x, 3),
            "y": round(y, 3),
            "width": round(w, 3),
            "height": round(h, 3),
            "fill": fill,
            "stroke": stroke,
            "stroke-width": sw if stroke else None,
            "rx": rx,
        },
    )


def line(x1, y1, x2, y2, stroke=INK, sw=0.35, dash=None):
    return tag(
        "line",
        {
            "x1": round(x1, 3),
            "y1": round(y1, 3),
            "x2": round(x2, 3),
            "y2": round(y2, 3),
            "stroke": stroke,
            "stroke-width": sw,
            "stroke-dasharray": dash,
        },
    )


def circle(cx, cy, r, fill, stroke=None, sw=0.25, opacity=1):
    return tag(
        "circle",
        {
            "cx": round(cx, 3),
            "cy": round(cy, 3),
            "r": round(r, 3),
            "fill": fill,
            "stroke": stroke,
            "stroke-width": sw if stroke else None,
            "opacity": opacity,
        },
    )


def polyline(points, stroke, sw=0.45, fill="none"):
    pts = " ".join(f"{round(x,3)},{round(y,3)}" for x, y in points)
    return tag("polyline", {"points": pts, "fill": fill, "stroke": stroke, "stroke-width": sw, "stroke-linejoin": "round", "stroke-linecap": "round"})


def polygon(points, fill, stroke=None, sw=0.25, opacity=1):
    pts = " ".join(f"{round(x,3)},{round(y,3)}" for x, y in points)
    return tag("polygon", {"points": pts, "fill": fill, "stroke": stroke, "stroke-width": sw if stroke else None, "opacity": opacity})


def panel_label(x, y, letter):
    return text(x, y, letter, size=2.8, weight="bold")


def axis(parts, x, y, w, h, xticks=None, yticks=None, xlabel="", ylabel="", xlim=(0, 1), ylim=(0, 1), yfmt=None, xfmt=None):
    parts.append(line(x, y + h, x + w, y + h, sw=0.45))
    parts.append(line(x, y, x, y + h, sw=0.45))
    for val in yticks or []:
        yy = y + h - (val - ylim[0]) / (ylim[1] - ylim[0]) * h
        parts.append(line(x - 1.2, yy, x, yy, sw=0.3))
        parts.append(line(x, yy, x + w, yy, stroke="#eeeeee", sw=0.25))
        label = yfmt(val) if yfmt else str(val)
        parts.append(text(x - 2, yy + 0.7, label, size=1.85, fill=MUTED, anchor="end"))
    for val in xticks or []:
        xx = x + (val - xlim[0]) / (xlim[1] - xlim[0]) * w
        parts.append(line(xx, y + h, xx, y + h + 1.2, sw=0.3))
        label = xfmt(val) if xfmt else str(val)
        parts.append(text(xx, y + h + 3.8, label, size=1.85, fill=MUTED, anchor="middle"))
    if xlabel:
        parts.append(text(x + w / 2, y + h + 8.2, xlabel, size=2.1, anchor="middle"))
    if ylabel:
        parts.append(text(x - 10.0, y + h / 2, ylabel, size=2.1, anchor="middle", rotate=-90))


def save_svg(name, width, height, parts):
    body = "\n".join(parts)
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}mm" height="{height}mm" viewBox="0 0 {width} {height}">
<metadata>{{"generator":"build_nature_style_publication_figures.py","style":"Nature-style editable SVG; sans-serif text; RGB"}}</metadata>
<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>
{body}
</svg>
'''
    path = SVG_DIR / name
    path.write_text(svg, encoding="utf-8")
    return path


def scale(v, lo, hi, out_lo, out_hi):
    if hi == lo:
        return (out_lo + out_hi) / 2
    return out_lo + (float(v) - lo) / (hi - lo) * (out_hi - out_lo)


def wrap_label(s, width=18):
    words = str(s).split()
    lines = []
    cur = ""
    for word in words:
        if len(cur) + len(word) + 1 <= width:
            cur = (cur + " " + word).strip()
        else:
            if cur:
                lines.append(cur)
            cur = word
    if cur:
        lines.append(cur)
    return lines


def figure1():
    counts = read("evidence_ladder_summary_counts.csv")
    sig = counts[counts.Entity.eq("signature")].copy()
    sig["short"] = sig["Evidence_ladder_label"].str.replace("_", " ", regex=False).str.replace("L0 ", "L0\n", regex=False)
    labels = ["L0\nnot reproducible", "L1\nformula", "L2\nscoreable", "L3\nscore portable", "L4\ncutoff", "L5\ndeployability\nevaluation"]
    vals = sig["Count"].tolist()
    parts = [panel_label(6, 9, "a")]
    x0, y0, bw, bh, gap = 14, 14, 24, 18, 4
    colors = [PALETTE["light_red"], PALETTE["light_blue"], PALETTE["light_blue"], PALETTE["light_teal"], PALETTE["light_orange"], "#eadff4"]
    for i, (lab, val) in enumerate(zip(labels, vals)):
        x = x0 + i * (bw + gap)
        parts.append(rect(x, y0, bw, bh, colors[i], stroke="#bdbdbd", rx=1.5))
        for k, line_ in enumerate(lab.split("\n")):
            parts.append(text(x + bw / 2, y0 + 4.2 + 3.5 * k, line_, size=1.9, weight="bold" if k == 0 else "normal", anchor="middle"))
        parts.append(text(x + bw / 2, y0 + bh - 2.2, str(int(val)), size=3.1, weight="bold", fill=PALETTE["blue"], anchor="middle"))
        if i < len(vals) - 1:
            parts.append(line(x + bw + 0.8, y0 + bh / 2, x + bw + gap - 1.0, y0 + bh / 2, stroke=MUTED, sw=0.45))
            parts.append(polygon([(x + bw + gap - 1, y0 + bh / 2 - 1.2), (x + bw + gap + 0.4, y0 + bh / 2), (x + bw + gap - 1, y0 + bh / 2 + 1.2)], MUTED))
    parts.append(text(6, 42, "Evidence levels are claims about deployability, not binary validity labels.", size=2.15, weight="bold"))

    parts.append(panel_label(6, 56, "b"))
    trans = read("evidence_ladder_transition_matrix.csv")
    stages = ["Full-text\narticles", "Formula-\nreproducible", "Externally\nscoreable", "Score\nportable", "Cutoff\ncandidate", "Deployability\nevaluation"]
    nums = [178, 57, 38, 28, 1, 1]
    x0, y0, ww, hh = 17, 61, 145, 34
    maxn = max(nums)
    for i, (st, n) in enumerate(zip(stages, nums)):
        x = x0 + i * (ww / (len(nums) - 1))
        r = 2.2 + 8 * math.sqrt(n / maxn)
        parts.append(circle(x, y0 + hh / 2, r, PALETTE["blue"] if i < 3 else PALETTE["teal"], opacity=0.82))
        for k, l in enumerate(st.split("\n")):
            parts.append(text(x, y0 + hh + 6 + k * 3.0, l, size=1.7, fill=INK, anchor="middle"))
        parts.append(text(x, y0 + hh / 2 + 0.7, str(n), size=2.0 if n > 9 else 2.3, weight="bold", fill="#ffffff", anchor="middle"))
        if i < len(nums) - 1:
            parts.append(line(x + r + 1, y0 + hh / 2, x0 + (i + 1) * (ww / (len(nums) - 1)) - r - 1, y0 + hh / 2, stroke="#9a9a9a", sw=0.45))

    parts.append(panel_label(6, 108, "c"))
    checklist = read("clinical_deployability_score_summary.csv")
    checklist["order"] = checklist["Article_ID"].map({a: i for i, a in enumerate(PRIMARY)})
    checklist = checklist.sort_values("order")
    x, y, w, h = 34, 108, 128, 30
    axis(parts, x, y, w, h, xticks=[8, 10, 12, 14, 16], xlim=(8, 16), ylim=(0, 1), xlabel="Deployability readiness score (0–18)")
    for i, r in enumerate(checklist.itertuples()):
        yy = y + 3 + i * 2.55
        xx = scale(r.Deployability_score_0_18, 8, 16, x, x + w)
        color = PALETTE["blue"] if r.Article_ID == "A131" else PALETTE["teal"] if r.Article_ID in STRONG else PALETTE["orange"] if r.Article_ID == "A098" else PALETTE["grey"]
        parts.append(text(x - 2, yy + 0.7, r.Article_ID, size=1.75, anchor="end"))
        parts.append(line(x, yy, xx, yy, stroke="#c8c8c8", sw=0.35))
        parts.append(circle(xx, yy, 1.35, color))
    save_svg("Figure1_evidence_ladder_workflow.svg", W, 148, parts)


def figure2():
    rq = read("reporting_quality_score_all_articles.csv")
    ryear = read("reporting_score_by_year.csv").dropna(subset=["Year"])
    parts = [panel_label(6, 8, "a")]
    items = ["gene_list", "coefficients", "formula", "numeric_cutoff", "cutoff_derivation", "training_cohort", "validation_cohort", "endpoint", "preprocessing", "code_data_availability"]
    labels = ["Gene\nlist", "Coeffi-\ncients", "Formula", "Numeric\ncutoff", "Cutoff\nderivation", "Training\ncohort", "Validation\ncohort", "Endpoint", "Pre-\nprocessing", "Code/data\navailability"]
    rep = rq.groupby(rq["Reproducible"].astype(str).str.lower().eq("yes"))[items].mean()
    x0, y0, cellw, cellh = 46, 13, 9.8, 9.0
    groups = [(False, "Not formula-reproducible"), (True, "Formula-reproducible")]
    for j, (flag, glab) in enumerate(groups):
        parts.append(text(14, y0 + j * cellh + 5.2, glab, size=2.0, anchor="start"))
        for i, item in enumerate(items):
            val = rep.loc[flag, item] if flag in rep.index else 0
            fill = interp_color(val, low=(243, 218, 218), high=(76, 159, 112))
            parts.append(rect(x0 + i * cellw, y0 + j * cellh, cellw - 0.5, cellh - 0.5, fill))
            parts.append(text(x0 + i * cellw + cellw / 2, y0 + j * cellh + 5.3, f"{val:.0%}", size=1.7, anchor="middle", fill="#ffffff" if val > 0.65 else INK, weight="bold"))
    for i, lab in enumerate(labels):
        for k, l in enumerate(lab.split("\n")):
            parts.append(text(x0 + i * cellw + cellw / 2, y0 - 5.5 + 2.8 * k, l, size=1.55, anchor="middle", fill=MUTED))

    parts.append(panel_label(6, 47, "b"))
    # Score distributions by ladder level
    levels = rq.groupby("Evidence_ladder_label")["Reporting_quality_score_0_10"].agg(["median", "mean", "count"]).reset_index()
    order = ["L0_non_reproducible_or_insufficiently_reported", "L1_formula_reproducible", "L2_externally_scoreable", "L3_score_portable", "L5_clinical_deployability_candidate"]
    levels = levels[levels["Evidence_ladder_label"].isin(order)].set_index("Evidence_ladder_label").loc[order].reset_index()
    x, y, w, h = 23, 54, 67, 36
    axis(parts, x, y, w, h, yticks=[0, 2, 4, 6, 8, 10], ylim=(0, 10), ylabel="Reporting score")
    bw = w / len(levels) * 0.58
    for i, r in enumerate(levels.itertuples()):
        xx = x + (i + 0.5) * w / len(levels) - bw / 2
        barh = r.mean / 10 * h
        parts.append(rect(xx, y + h - barh, bw, barh, PALETTE["blue"] if i < 3 else PALETTE["teal"], rx=0.5))
        parts.append(text(xx + bw / 2, y + h + 4, f"L{i}", size=1.8, anchor="middle"))
        parts.append(text(xx + bw / 2, y + h - barh - 1.2, str(int(r.count)), size=1.6, anchor="middle", fill=MUTED))

    parts.append(panel_label(101, 47, "c"))
    x, y, w, h = 112, 54, 58, 36
    years = ryear["Year"].astype(int).tolist()
    mean = ryear["Mean_reporting_score"].tolist()
    formula = ryear["Formula_reproducible_count"].tolist()
    axis(parts, x, y, w, h, yticks=[0, 2, 4, 6, 8, 10], ylim=(0, 10), ylabel="", xticks=[min(years), 2020, max(years)], xlim=(min(years), max(years)), xlabel="Publication year")
    pts = [(scale(yr, min(years), max(years), x, x + w), y + h - sc / 10 * h) for yr, sc in zip(years, mean)]
    parts.append(polyline(pts, PALETTE["blue"], sw=0.55))
    for (xx, yy), n in zip(pts, formula):
        parts.append(circle(xx, yy, 0.9 + 0.25 * n, PALETTE["orange"], opacity=0.72))
    parts.append(text(112, 95, "Dot size: formula-reproducible articles", size=1.8, fill=MUTED))
    save_svg("Figure2_reporting_reproducibility_landscape.svg", W, 102, parts)


def interp_color(v, low=(255, 255, 255), high=(0, 0, 0)):
    v = max(0, min(1, float(v)))
    rgb = tuple(round(low[i] + (high[i] - low[i]) * v) for i in range(3))
    return "#%02x%02x%02x" % rgb


def figure3():
    perf = read("all_57_GEO_performance_summary.csv")
    parts = [panel_label(6, 8, "a")]
    x, y, w, h = 24, 14, 76, 55
    axis(parts, x, y, w, h, xticks=[0.50, 0.55, 0.60, 0.65, 0.70], yticks=[0, 0.25, 0.5, 0.75, 1], xlim=(0.48, 0.72), ylim=(0, 1), xlabel="Median C-index", ylabel="Continuous p<0.05 rate", xfmt=lambda v: f"{v:.2f}", yfmt=lambda v: f"{v:.2f}")
    for r in perf.itertuples():
        xx = scale(r.median_C_index, 0.48, 0.72, x, x + w)
        yy = y + h - scale(r.continuous_significant_rate, 0, 1, 0, h)
        primary = r.Article_ID in PRIMARY
        color = PALETTE["blue"] if r.Article_ID == "A131" else PALETTE["teal"] if r.Article_ID in STRONG else PALETTE["orange"] if r.Article_ID == "A098" else PALETTE["grey"]
        parts.append(circle(xx, yy, 1.45 if primary else 0.9, color, opacity=0.95 if primary else 0.45))
        if primary:
            parts.append(text(xx + 1.5, yy - 1.0, r.Article_ID, size=1.65, weight="bold", fill=INK))
    parts.append(text(28, 76, "Primary candidates labelled; unlabelled points are other formula-reproducible signatures.", size=1.8, fill=MUTED))

    parts.append(panel_label(111, 8, "b"))
    score = read("all_57_signature_scoreability.csv")
    statuses = [
        ("ok_tests", "Evaluable", PALETTE["green"]),
        ("insufficient_gene_coverage_tests", "Insufficient gene coverage", PALETTE["red"]),
        ("skipped_training_overlap_tests", "Training overlap excluded", PALETTE["orange"]),
        ("single_group_after_cutoff_tests", "Single-group collapse", PALETTE["purple"]),
        ("score_contains_missing_tests", "Missing score component", PALETTE["grey"]),
    ]
    totals = [int(score[c].sum()) for c, _, _ in statuses]
    x0, y0, maxw, bh, gap = 116, 18, 52, 6.2, 3.4
    maxv = max(totals)
    for i, ((col, lab, color), val) in enumerate(zip(statuses, totals)):
        yy = y0 + i * (bh + gap)
        parts.append(text(x0 - 2, yy + 4.4, lab, size=1.75, anchor="end"))
        parts.append(rect(x0, yy, maxw * val / maxv, bh, color))
        parts.append(text(x0 + maxw * val / maxv + 1.5, yy + 4.4, str(val), size=1.8, weight="bold"))

    parts.append(panel_label(6, 91, "c"))
    top = perf.sort_values(["Evidence_ladder_level", "median_C_index"], ascending=[False, False]).head(15)
    x, y, w, h = 37, 93, 124, 42
    axis(parts, x, y, w, h, xticks=[0.50, 0.60, 0.70], xlim=(0.48, 0.72), ylim=(0, 1), xlabel="Median C-index")
    for i, r in enumerate(top.itertuples()):
        yy = y + 2 + i * 2.55
        xx = scale(r.median_C_index, 0.48, 0.72, x, x + w)
        parts.append(text(x - 2, yy + 0.7, r.Article_ID, size=1.65, anchor="end"))
        parts.append(line(x, yy, xx, yy, stroke="#cfcfcf", sw=0.32))
        color = PALETTE["blue"] if r.Article_ID == "A131" else PALETTE["teal"] if r.Article_ID in STRONG else PALETTE["grey"]
        parts.append(circle(xx, yy, 1.1, color))
    save_svg("Figure3_all57_score_portability_landscape.svg", W, 145, parts)


def figure4():
    cutoff = read("cutoff_transportability_summary.csv")
    drift = read("cutoff_percentile_drift_long.csv")
    parts = [panel_label(6, 8, "a")]
    labels, vals, colors = [], [], []
    for r in cutoff.itertuples():
        if not r.reported_numeric_cutoff:
            lab = "No numeric cutoff"
            col = PALETTE["grey"]
        elif r.reported_cutoff_applicable_to_GEO and not r.single_group_collapse_observed:
            lab = "Applicable; no collapse"
            col = PALETTE["green"]
        elif r.single_group_collapse_observed:
            lab = "Collapse observed"
            col = PALETTE["red"]
        else:
            lab = "Numeric but not GEO-applicable"
            col = PALETTE["orange"]
        labels.append(lab)
        vals.append(int(r.Signature_count))
        colors.append(col)
    # combine repeated labels
    comb = {}
    for lab, val, col in zip(labels, vals, colors):
        comb.setdefault(lab, [0, col])[0] += val
    x0, y0, maxw, bh, gap = 65, 15, 86, 8, 5
    maxv = max(v for v, _ in comb.values())
    for i, (lab, (val, col)) in enumerate(comb.items()):
        yy = y0 + i * (bh + gap)
        parts.append(text(x0 - 2, yy + 5.3, lab, size=2.0, anchor="end"))
        parts.append(rect(x0, yy, maxw * val / maxv, bh, col))
        parts.append(text(x0 + maxw * val / maxv + 1.5, yy + 5.3, str(val), size=2.0, weight="bold"))

    parts.append(panel_label(6, 66, "b"))
    a097 = drift[drift["Article_ID"].eq("A097")].copy()
    x, y, w, h = 39, 69, 59, 40
    axis(parts, x, y, w, h, xticks=[0, 0.5, 1.0], yticks=[], xlim=(0, 1), ylim=(0, 1), xlabel="High-risk proportion under cutoff 1.66")
    for i, r in enumerate(a097.itertuples()):
        yy = y + 3 + i * 4.6
        hp = r.High_n / r.N if r.N else np.nan
        xx = scale(hp, 0, 1, x, x + w)
        parts.append(text(x - 2, yy + 0.7, r.GSE, size=1.65, anchor="end"))
        parts.append(line(x, yy, xx, yy, stroke="#cccccc", sw=0.35))
        parts.append(circle(xx, yy, 1.2, PALETTE["red"]))
    parts.append(text(107, 79, "A097: verified TCGA/training cutoff,", size=2.1, weight="bold"))
    parts.append(text(107, 84, "but all tested GEO cohorts collapse", size=2.1, weight="bold", fill=PALETTE["red"]))
    parts.append(text(107, 89, "to all-high or all-low groups.", size=2.1, weight="bold", fill=PALETTE["red"]))

    parts.append(panel_label(6, 123, "c"))
    prop = read("reported_cutoff_high_risk_proportion_long.csv")
    prop = prop[prop["reported_cutoff_high_risk_proportion"].notna()].copy()
    x, y, w, h = 29, 124, 126, 22
    axis(parts, x, y, w, h, xticks=[0, 0.25, 0.5, 0.75, 1.0], xlim=(0, 1), ylim=(0, 1), xlabel="Reported-cutoff high-risk proportion across GEO tests")
    parts.append(line(scale(0.5, 0, 1, x, x + w), y, scale(0.5, 0, 1, x, x + w), y + h, stroke=PALETTE["blue"], sw=0.4, dash="1.2,1.2"))
    for i, r in enumerate(prop.itertuples()):
        xx = scale(r.reported_cutoff_high_risk_proportion, 0, 1, x, x + w)
        yy = y + 3 + (i % 6) * 2.8 + (i // 6) * 0.1
        color = PALETTE["red"] if r.Article_ID == "A097" else PALETTE["orange"] if r.Status != "ok" else PALETTE["green"]
        parts.append(circle(xx, yy, 0.9, color, opacity=0.8))
    save_svg("Figure4_cutoff_transportability_failure_anatomy.svg", W, 155, parts)


def figure5():
    auc = read("time_dependent_auc_primary_candidates.csv")
    brier = read("brier_score_primary_candidates.csv")
    inc = read("incremental_cindex_by_signature_cohort.csv")
    lrt = read("likelihood_ratio_test_clinical_vs_clinical_plus_risk.csv")
    rows = []
    for art in PRIMARY:
        rows.append({
            "Article_ID": art,
            "AUC_1y": auc[(auc.Article_ID.eq(art)) & (auc.time_years.eq(1.0))]["dynamic_AUC_complete_case"].median(),
            "AUC_3y": auc[(auc.Article_ID.eq(art)) & (auc.time_years.eq(3.0))]["dynamic_AUC_complete_case"].median(),
            "AUC_5y": auc[(auc.Article_ID.eq(art)) & (auc.time_years.eq(5.0))]["dynamic_AUC_complete_case"].median(),
            "Brier_3y": brier[(brier.Article_ID.eq(art)) & (brier.time_years.eq(3.0))]["IPCW_Brier_score_external_recalibrated"].median(),
            "Brier_5y": brier[(brier.Article_ID.eq(art)) & (brier.time_years.eq(5.0))]["IPCW_Brier_score_external_recalibrated"].median(),
            "Delta_C": inc[inc.Article_ID.eq(art)]["Delta_C_index"].median(),
            "LRT_rate": (lrt[lrt.Article_ID.eq(art)]["LR_test_p"] < 0.05).mean(),
        })
    df = pd.DataFrame(rows)
    df.to_csv(SRC_DIR / "Figure5_primary_prediction_performance_summary.csv", index=False, encoding="utf-8")
    parts = [panel_label(6, 8, "a")]
    metrics = ["AUC_1y", "AUC_3y", "AUC_5y", "Brier_3y", "Brier_5y", "Delta_C", "LRT_rate"]
    labels = ["AUC\n1 y", "AUC\n3 y", "AUC\n5 y", "Brier\n3 y", "Brier\n5 y", "ΔC", "LRT\np<0.05"]
    x0, y0, cw, ch = 38, 15, 14, 5.8
    for i, art in enumerate(PRIMARY):
        parts.append(text(x0 - 2, y0 + i * ch + 4, art, size=1.8, anchor="end", weight="bold" if art == "A131" else "normal"))
        for j, m in enumerate(metrics):
            val = df.loc[df.Article_ID.eq(art), m].iloc[0]
            if "Brier" in m:
                norm = 1 - max(0, min(1, (val - 0.12) / 0.18))
            elif m == "Delta_C":
                norm = max(0, min(1, val / 0.12))
            else:
                norm = max(0, min(1, (val - 0.50) / 0.25))
            fill = interp_color(norm, low=(242, 242, 242), high=(59, 111, 182))
            parts.append(rect(x0 + j * cw, y0 + i * ch, cw - 0.4, ch - 0.4, fill, stroke="#ffffff", sw=0.2))
            parts.append(text(x0 + j * cw + cw / 2, y0 + i * ch + 3.8, fmt(val, 2), size=1.55, anchor="middle", fill="#ffffff" if norm > 0.55 else INK))
    for j, lab in enumerate(labels):
        for k, l in enumerate(lab.split("\n")):
            parts.append(text(x0 + j * cw + cw / 2, y0 - 5 + k * 2.6, l, size=1.65, anchor="middle", fill=MUTED))
    parts.append(text(38, 81, "Brier values are external recalibration assessments; lower is better.", size=1.85, fill=MUTED))

    parts.append(panel_label(6, 96, "b"))
    x, y, w, h = 35, 98, 54, 36
    axis(parts, x, y, w, h, yticks=[0.50, 0.60, 0.70, 0.80], ylim=(0.5, 0.82), xticks=[1, 3, 5], xlim=(0.6, 5.4), xlabel="Years", ylabel="Median AUC", yfmt=lambda v: f"{v:.2f}")
    for art in ["A131", "A129", "A067", "A097", "A168", "A001", "A059"]:
        sub = df[df.Article_ID.eq(art)].iloc[0]
        pts = [(scale(1, 0.6, 5.4, x, x + w), scale(sub.AUC_1y, 0.82, 0.5, y, y + h)), (scale(3, 0.6, 5.4, x, x + w), scale(sub.AUC_3y, 0.82, 0.5, y, y + h)), (scale(5, 0.6, 5.4, x, x + w), scale(sub.AUC_5y, 0.82, 0.5, y, y + h))]
        color = PALETTE["blue"] if art == "A131" else PALETTE["teal"]
        parts.append(polyline(pts, color, sw=0.42))
        for xx, yy in pts:
            parts.append(circle(xx, yy, 0.8, color))

    parts.append(panel_label(103, 96, "c"))
    d = df[["Article_ID", "Delta_C"]].copy()
    x, y, w, h = 124, 98, 40, 36
    axis(parts, x, y, w, h, xticks=[0, 0.05, 0.10, 0.15], xlim=(0, 0.15), ylim=(0, 1), xlabel="Median ΔC-index")
    for i, r in enumerate(d.itertuples()):
        yy = y + 2 + i * 3.1
        xx = scale(r.Delta_C, 0, 0.15, x, x + w)
        parts.append(text(x - 2, yy + 0.7, r.Article_ID, size=1.65, anchor="end"))
        parts.append(line(x, yy, xx, yy, stroke="#cccccc", sw=0.32))
        parts.append(circle(xx, yy, 1.0, PALETTE["blue"] if r.Article_ID == "A131" else PALETTE["teal"]))
    save_svg("Figure5_primary_prediction_performance.svg", W, 145, parts)


def figure6():
    heat = read("clinical_deployability_heatmap_data.csv")
    domains = list(heat["Deployability_domain"].drop_duplicates())
    # Short labels for publication.
    domain_short = {
        "complete_formula": "Formula",
        "gene_availability": "Genes",
        "endpoint_compatibility_with_GEO_OS": "Endpoint",
        "reported_numeric_cutoff": "Cutoff\nreported",
        "strict_cutoff_testability": "Strict\ncutoff",
        "no_single_group_collapse": "No\ncollapse",
        "continuous_GEO_performance": "GEO\nscore",
        "clinical_adjusted_signal": "Adjusted\nsignal",
        "delta_C_index_positive": "ΔC",
        "time_AUC_available": "AUC",
        "Brier_available": "Brier",
        "external_recalibration_available": "Recal.",
        "TCGA_formula_direction_sanity": "TCGA\nsanity",
        "leave_one_cohort_out_direction_retained": "LOCO",
        "original_code_or_calculator": "Code/\ncalc.",
        "assay_standardization": "Assay\nstd.",
        "prospective_validation": "Prospective",
        "clinical_threshold_justification": "Threshold\nrationale",
    }
    val = {(r.Article_ID, r.Deployability_domain): r.Heatmap_value for r in heat.itertuples()}
    parts = [panel_label(6, 8, "a")]
    x0, y0, cw, ch = 31, 24, 7.5, 6.2
    for j, dom in enumerate(domains):
        label = domain_short.get(dom, dom)
        for k, l in enumerate(label.split("\n")):
            parts.append(text(x0 + j * cw + cw / 2, y0 - 10 + k * 2.6, l, size=1.45, anchor="middle", fill=MUTED, rotate=-35))
    for i, art in enumerate(PRIMARY):
        y = y0 + i * ch
        parts.append(text(x0 - 2, y + 4.1, art, size=1.85, anchor="end", weight="bold" if art == "A131" else "normal"))
        for j, dom in enumerate(domains):
            v = val.get((art, dom), 0)
            fill = "#4c9f70" if v >= 0.99 else "#d8b545" if v >= 0.49 else "#9bbce3" if v >= 0.24 else "#d9d9d9"
            parts.append(rect(x0 + j * cw, y, cw - 0.35, ch - 0.35, fill, stroke="#ffffff", sw=0.18))
    parts.append(text(31, 92, "Green: yes; gold: partial; blue: manual check; grey: no.", size=1.85, fill=MUTED))

    parts.append(panel_label(6, 100, "b"))
    checklist = read("clinical_deployability_score_summary.csv")
    checklist["order"] = checklist["Article_ID"].map({a: i for i, a in enumerate(PRIMARY)})
    checklist = checklist.sort_values("order")
    x, y, w, h = 36, 101, 118, 32
    axis(parts, x, y, w, h, xticks=[0, 6, 12, 18], xlim=(0, 18), ylim=(0, 1), xlabel="Deployability readiness score")
    for i, r in enumerate(checklist.itertuples()):
        yy = y + 2 + i * 2.8
        xx = scale(r.Deployability_score_0_18, 0, 18, x, x + w)
        parts.append(text(x - 2, yy + 0.7, r.Article_ID, size=1.65, anchor="end"))
        parts.append(rect(x, yy - 1, xx - x, 2.0, PALETTE["blue"] if r.Article_ID == "A131" else PALETTE["teal"]))
        parts.append(text(xx + 1.3, yy + 0.7, fmt(r.Deployability_score_0_18, 1), size=1.6, fill=INK))
    save_svg("Figure6_clinical_deployability_gap_heatmap.svg", W, 142, parts)


def copy_source_data():
    mapping = {
        "Figure1_evidence_ladder_workflow": ["evidence_ladder_summary_counts.csv", "evidence_ladder_transition_matrix.csv", "clinical_deployability_score_summary.csv"],
        "Figure2_reporting_reproducibility_landscape": ["reporting_quality_score_all_articles.csv", "reporting_score_by_year.csv"],
        "Figure3_all57_score_portability_landscape": ["all_57_GEO_performance_summary.csv", "all_57_signature_scoreability.csv"],
        "Figure4_cutoff_transportability_failure_anatomy": ["cutoff_transportability_summary.csv", "cutoff_percentile_drift_long.csv", "reported_cutoff_high_risk_proportion_long.csv"],
        "Figure5_primary_prediction_performance": ["time_dependent_auc_primary_candidates.csv", "brier_score_primary_candidates.csv", "incremental_cindex_by_signature_cohort.csv", "likelihood_ratio_test_clinical_vs_clinical_plus_risk.csv"],
        "Figure6_clinical_deployability_gap_heatmap": ["clinical_deployability_heatmap_data.csv", "clinical_deployability_score_summary.csv"],
    }
    rows = []
    for fig, files in mapping.items():
        figdir = SRC_DIR / fig
        figdir.mkdir(parents=True, exist_ok=True)
        for file in files:
            shutil.copy2(AUDIT / file, figdir / file)
            rows.append({"Figure": fig, "Source_file": file})
    pd.DataFrame(rows).to_csv(SRC_DIR / "figure_source_data_index.csv", index=False, encoding="utf-8")


def convert_with_pymupdf():
    zoom = 600 / 72
    matrix = fitz.Matrix(zoom, zoom)
    for svg in SVG_DIR.glob("*.svg"):
        doc = fitz.open(str(svg))
        (PDF_DIR / f"{svg.stem}.pdf").write_bytes(doc.convert_to_pdf())
        pix = doc[0].get_pixmap(matrix=matrix, alpha=False)
        png_path = PNG_DIR / f"{svg.stem}.png"
        pix.save(str(png_path))
        img = Image.open(png_path)
        img.save(TIFF_DIR / f"{svg.stem}.tiff", compression="tiff_lzw", dpi=(600, 600))
        doc.close()


def refresh_package():
    fig_pkg = PKG / "figures_nature_style"
    if fig_pkg.exists():
        shutil.rmtree(fig_pkg)
    shutil.copytree(OUT, fig_pkg)
    if ZIP.exists():
        ZIP.unlink()
    with zipfile.ZipFile(ZIP, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in PKG.rglob("*"):
            if path.is_file():
                zf.write(path, path.relative_to(PKG.parent))


def validate():
    expected = [
        "Figure1_evidence_ladder_workflow",
        "Figure2_reporting_reproducibility_landscape",
        "Figure3_all57_score_portability_landscape",
        "Figure4_cutoff_transportability_failure_anatomy",
        "Figure5_primary_prediction_performance",
        "Figure6_clinical_deployability_gap_heatmap",
    ]
    rows = []
    for base in expected:
        row = {"Figure": base}
        for ext, directory in [("svg", SVG_DIR), ("png", PNG_DIR), ("tiff", TIFF_DIR), ("pdf", PDF_DIR)]:
            p = directory / f"{base}.{ext}"
            if not p.exists() or p.stat().st_size == 0:
                raise RuntimeError(f"Missing {p}")
            row[f"{ext}_bytes"] = p.stat().st_size
            if ext == "svg":
                text_ = p.read_text(encoding="utf-8")
                if "<text" not in text_:
                    raise RuntimeError(f"SVG text not editable: {p}")
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "publication_figure_export_manifest.csv", index=False, encoding="utf-8")


def main():
    for d in [SVG_DIR, PNG_DIR, TIFF_DIR, PDF_DIR, HTML_DIR, SRC_DIR]:
        d.mkdir(parents=True, exist_ok=True)
    copy_source_data()
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    figure6()
    convert_with_pymupdf()
    validate()
    refresh_package()
    print(f"Nature-style publication figures: {OUT}")
    print(f"Package refreshed: {ZIP}")


if __name__ == "__main__":
    main()
