from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "ai_innovation_main_figure"
OUT.mkdir(parents=True, exist_ok=True)
FIG_UPLOAD = ROOT / "outputs" / "bmc_ai_second_stage_performance_package" / "figures_for_upload"
FIG_UPLOAD.mkdir(parents=True, exist_ok=True)
SOURCE_DATA = ROOT / "source_data"
SOURCE_DATA.mkdir(exist_ok=True)


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
mpl.rcParams.update(
    {
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "legend.frameon": False,
    }
)


COL = {
    "dark": "#2B2B2B",
    "muted": "#707070",
    "grid": "#D8D8D8",
    "blue": "#3B6EA8",
    "blue_light": "#D8E6F3",
    "teal": "#4B9A9A",
    "teal_light": "#D8EEEE",
    "gold": "#D9A441",
    "gold_light": "#F4E6C3",
    "red": "#B75A5A",
    "red_light": "#F1D8D8",
    "green": "#4D8A5A",
    "green_light": "#DDEBDD",
    "purple": "#7A5FA0",
    "purple_light": "#E7DFF0",
}


def external_deepsurv_linear():
    df = pd.read_csv(ROOT / "outputs" / "ai_deepsurv_signature_benchmark" / "AI_deepsurv_external_validation_results.csv")
    df = df[df["Dataset"].eq("GEO_external")].copy()
    return {
        "DeepSurv": df[df["Model"].eq("AI_DeepSurv_signature_gene_model")].assign(strategy_label="DeepSurv"),
        "Linear": df[df["Model"].eq("Penalized_linear_Cox_signature_gene_baseline")].assign(strategy_label="Penalized linear Cox"),
    }


def optimized_linear():
    df = pd.read_csv(ROOT / "outputs" / "ai_signature_gene_optimization" / "AI_selected_strategy_external_results_by_cohort.csv")
    df = df[df["selection_rule"].eq("posthoc_best_external_linear")].copy()
    return df.assign(strategy_label="Optimized linear\n150 genes")


def second_stage():
    df = pd.read_csv(ROOT / "outputs" / "ai_second_stage_performance_optimization" / "AI_second_stage_strategy_by_cohort.csv")
    core5 = df[df["strategy"].eq("core_tcga_locked_external_median_logHR_top5")].copy()
    core40 = df[df["strategy"].eq("core_cohort_rank_normal_external_median_logHR_top40")].copy()
    return core5.assign(strategy_label="5-gene core\nupper bound"), core40.assign(strategy_label="40-gene rank-normal\ncore")


def summarize_strategy(label, df):
    return {
        "strategy_label": label,
        "median_C_index": df["C_index"].median(),
        "min_C_index": df["C_index"].min(),
        "max_C_index": df["C_index"].max(),
        "direction_HR_gt1_n": int((df["HR_per_SD"] > 1).sum()),
        "significant_p05_n": int((df["Cox_p"] < 0.05).sum()),
        "cohorts_C_ge_055_n": int((df["C_index"] >= 0.55).sum()),
        "external_tests": len(df),
    }


def panel_label(ax, label):
    ax.text(-0.05, 1.05, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")


def draw_box(ax, xy, wh, text, fc, ec, fontsize=6.2):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=0.8,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize, color=COL["dark"], linespacing=1.15)
    return patch


def arrow(ax, start, end, color=COL["muted"]):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=8, lw=0.8, color=color))


def main():
    dl = external_deepsurv_linear()
    optlin = optimized_linear()
    core5, core40 = second_stage()

    strategy_frames = [
        ("DeepSurv", dl["DeepSurv"], COL["red"]),
        ("Penalized linear Cox", dl["Linear"], COL["blue"]),
        ("Optimized linear\n150 genes", optlin, COL["teal"]),
        ("Core-5\nupper bound", core5, COL["gold"]),
        ("Core-40\nrank-normal", core40, COL["green"]),
    ]
    summary = pd.DataFrame([summarize_strategy(label, df) for label, df, _ in strategy_frames])
    summary.to_csv(OUT / "Fig6_AI_innovation_strategy_summary.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(SOURCE_DATA / "Fig6_AI_innovation_strategy_summary.csv", index=False, encoding="utf-8-sig")

    heat_rows = []
    for label, df, _ in strategy_frames:
        d = df.copy()
        d["strategy_label"] = label
        heat_rows.append(d[["strategy_label", "GSE", "N", "Events", "C_index", "HR_per_SD", "Cox_p"]])
    heat = pd.concat(heat_rows, ignore_index=True)
    heat.to_csv(OUT / "Fig6_AI_innovation_cohort_results.csv", index=False, encoding="utf-8-sig")
    heat.to_csv(SOURCE_DATA / "Fig6_AI_innovation_cohort_results.csv", index=False, encoding="utf-8-sig")

    stability = pd.read_csv(ROOT / "outputs" / "ai_signature_gene_optimization" / "AI_core_feature_external_stability.csv")
    gene_plot = stability.head(20).copy()
    gene_plot.to_csv(OUT / "Fig6_AI_innovation_core_gene_stability.csv", index=False, encoding="utf-8-sig")
    gene_plot.to_csv(SOURCE_DATA / "Fig6_AI_innovation_core_gene_stability.csv", index=False, encoding="utf-8-sig")

    fig = plt.figure(figsize=(7.2, 6.0), constrained_layout=False)
    gs = gridspec.GridSpec(
        2,
        2,
        figure=fig,
        width_ratios=[1.02, 1.2],
        height_ratios=[0.92, 1.08],
        wspace=0.44,
        hspace=0.48,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # Panel a: schematic route.
    ax_a.set_axis_off()
    panel_label(ax_a, "a")
    ax_a.text(0.0, 0.98, "AI layer converts an audit bottleneck into a testable optimization route", fontsize=7.2, fontweight="bold", va="top")
    draw_box(ax_a, (0.02, 0.72), (0.38, 0.16), "Published LUAD\nsignatures\n178 articles", COL["blue_light"], COL["blue"])
    draw_box(ax_a, (0.02, 0.48), (0.38, 0.16), "Deployability attrition\n57 formula-reproducible\n10 primary candidates", "#F5F5F5", COL["muted"])
    draw_box(ax_a, (0.02, 0.24), (0.38, 0.16), "Clinical deployment\nstill unsupported\n0 deployment-ready", COL["red_light"], COL["red"])
    draw_box(ax_a, (0.58, 0.72), (0.38, 0.16), "DeepSurv stress test\nexternal median C=0.476", COL["red_light"], COL["red"])
    draw_box(ax_a, (0.58, 0.48), (0.38, 0.16), "Parameter tuning\nsimple + shrinkage\nmedian C=0.584", COL["teal_light"], COL["teal"])
    draw_box(ax_a, (0.58, 0.24), (0.38, 0.16), "Stable core features\nrank-normalized 40-gene core\nmedian C=0.663", COL["green_light"], COL["green"])
    for y1, y2 in [(0.72, 0.64), (0.48, 0.40)]:
        arrow(ax_a, (0.21, y1), (0.21, y2))
        arrow(ax_a, (0.77, y1), (0.77, y2))
    arrow(ax_a, (0.40, 0.56), (0.58, 0.56), COL["purple"])
    ax_a.text(0.49, 0.61, "AI-guided\nstress test", ha="center", va="bottom", fontsize=6.2, color=COL["purple"])
    ax_a.text(0.02, 0.075, "Innovation: deployability-aware\noptimization, not a new locked\nclinical LUAD model.", fontsize=5.8, color=COL["muted"], linespacing=1.15)

    # Panel b: performance ladder.
    panel_label(ax_b, "b")
    order = [x[0] for x in strategy_frames]
    colors = [x[2] for x in strategy_frames]
    y = np.arange(len(order))
    vals = summary.set_index("strategy_label").loc[order, "median_C_index"].to_numpy()
    lo = summary.set_index("strategy_label").loc[order, "min_C_index"].to_numpy()
    hi = summary.set_index("strategy_label").loc[order, "max_C_index"].to_numpy()
    ax_b.barh(y, vals - 0.5, left=0.5, color=colors, alpha=0.86, edgecolor="white", height=0.62)
    ax_b.errorbar(vals, y, xerr=[vals - lo, hi - vals], fmt="none", ecolor=COL["dark"], elinewidth=0.8, capsize=2)
    for yi, v in zip(y, vals):
        ax_b.text(v + 0.006, yi, f"{v:.3f}", va="center", fontsize=6.5)
    ax_b.axvline(0.5, color=COL["muted"], ls="--", lw=0.8)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(order)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0.44, 0.75)
    ax_b.set_xlabel("External GEO C-index, median (range)")
    ax_b.set_title("External performance improves when stable features and shrinkage replace depth", fontsize=7.2, loc="left", pad=4)
    ax_b.grid(axis="x", color=COL["grid"], lw=0.5)

    # Panel c: cohort heatmap.
    panel_label(ax_c, "c")
    chosen = ["DeepSurv", "Penalized linear Cox", "Optimized linear\n150 genes", "Core-40\nrank-normal"]
    cohorts = sorted(heat["GSE"].unique())
    mat = np.full((len(cohorts), len(chosen)), np.nan)
    pmat = np.full_like(mat, np.nan)
    for i, gse in enumerate(cohorts):
        for j, label in enumerate(chosen):
            row = heat[(heat["GSE"].eq(gse)) & (heat["strategy_label"].eq(label))]
            if len(row):
                mat[i, j] = float(row.iloc[0]["C_index"])
                pmat[i, j] = float(row.iloc[0]["Cox_p"])
    im = ax_c.imshow(mat, vmin=0.45, vmax=0.74, cmap="YlGnBu", aspect="auto")
    ax_c.set_xticks(np.arange(len(chosen)))
    ax_c.set_xticklabels(chosen, rotation=35, ha="right")
    ax_c.set_yticks(np.arange(len(cohorts)))
    ax_c.set_yticklabels(cohorts)
    ax_c.set_title("Core-feature strategy stabilizes performance across cohorts", fontsize=7.2, loc="left", pad=4)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            txt = f"{mat[i, j]:.2f}"
            mark = "*" if pmat[i, j] < 0.05 else ""
            ax_c.text(j, i, txt + mark, ha="center", va="center", fontsize=5.8, color="#1F1F1F")
    cbar = fig.colorbar(im, ax=ax_c, fraction=0.04, pad=0.02)
    cbar.ax.tick_params(labelsize=6)
    cbar.set_label("C-index", fontsize=6)
    ax_c.text(0.0, -0.33, "* Cox p < 0.05 for HR per SD; all strategies evaluated in the same nine GEO OS cohorts.", transform=ax_c.transAxes, fontsize=5.8, color=COL["muted"])

    # Panel d: stable core genes.
    panel_label(ax_d, "d")
    gp = gene_plot.iloc[::-1].copy()
    y2 = np.arange(len(gp))
    hr = gp["external_median_HR_per_SD"].to_numpy(float)
    beta_col = np.where(hr >= 1, COL["red"], COL["blue"])
    ax_d.barh(y2, np.log(hr), color=beta_col, alpha=0.78, height=0.64)
    ax_d.axvline(0, color=COL["dark"], lw=0.8)
    ax_d.set_yticks(y2)
    ax_d.set_yticklabels(gp["gene"])
    ax_d.set_xlabel("External median log(HR per SD)")
    ax_d.set_title("Direction-consistent genes define the recoverable signal", fontsize=7.2, loc="left", pad=4)
    for yi, same, p05 in zip(y2, gp["external_same_direction_n"], gp["external_p05_n"]):
        ax_d.text(ax_d.get_xlim()[1] if ax_d.get_xlim()[1] > 0 else 0.42, yi, f"{int(same)}/9 | p<0.05: {int(p05)}", va="center", ha="right", fontsize=5.6, color=COL["muted"])
    ax_d.grid(axis="x", color=COL["grid"], lw=0.5)
    ax_d.text(0.02, -0.17, "Red: higher expression associated with higher risk; blue: protective direction.", transform=ax_d.transAxes, fontsize=5.8, color=COL["muted"])

    for ax in [ax_b, ax_c, ax_d]:
        ax.tick_params(labelsize=6.2)

    fig.suptitle("AI/deep-learning augmentation exposes the transportability bottleneck and a stable-feature route", x=0.02, y=0.995, ha="left", fontsize=9.5, fontweight="bold")
    fig.subplots_adjust(top=0.91, left=0.08, right=0.98, bottom=0.10)

    for stem_dir in [OUT, FIG_UPLOAD]:
        stem = stem_dir / "Figure_6"
        fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
        fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
        fig.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
        fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "figure": "Figure 6",
        "core_conclusion": "AI/deep-learning augmentation exposes transportability failure of deeper models and identifies a stable-feature optimization route.",
        "archetype": "schematic-led composite",
        "backend": "Python/matplotlib",
        "source_data": [
            str(SOURCE_DATA / "Fig6_AI_innovation_strategy_summary.csv"),
            str(SOURCE_DATA / "Fig6_AI_innovation_cohort_results.csv"),
            str(SOURCE_DATA / "Fig6_AI_innovation_core_gene_stability.csv"),
        ],
    }
    (OUT / "Fig6_AI_innovation_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
