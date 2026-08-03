"""
engulfment_plots.py
-------------------
Generate SuperPlots from a pre-computed engulfment_results.csv.

Provides:
  * the primary genotype comparison (WT vs Opn4 KO), as in the main pipeline, and
  * an inline breakdown by sex, requested during peer review.

Group coding: in the CSV, groups are stored as A / B (kept separate during
blinded analysis). A = WT, B = Opn4 KO.

Channel coding (kept separate during blinded analysis):
  red   = VGluT2      (engulfed presynaptic puncta)
  green = microglia        (the engulfing volume; reference for engulfment %)
  blue  = CTB         (engulfed anterograde tracer)

Engulfment is the fraction of microglial (green) volume also occupied by VGluT2 (red)
or CTB (blue) signal.

Author: Colenso M. Speer (cspeer@umd.edu)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# --------------------------------------------------------------------------
# Genotype coding. A / B are the blinded labels stored in the CSV.
# --------------------------------------------------------------------------
GROUP_LABELS = {"A": "WT", "B": "Opn4 KO"}

# Sex assignment per animal (the `replicate` column). Edit if reassigned.
SEX_MAP = {
    "Sample_1A": "M", "Sample_2A": "M", "Sample_3A": "M",
    "Sample_4A": "M", "Sample_5A": "M", "Sample_6A": "F", "Sample_7A": "F",
    "Sample_1B": "M", "Sample_2B": "M", "Sample_3B": "M", "Sample_4B": "M",
    "Sample_5B": "M", "Sample_6B": "M", "Sample_7B": "F",
}

# Metrics, with display titles and y-axis labels (channels correctly named).
METRIC_INFO = {
    "red_VGluT2_engulfment_pct":      ("VGluT2 Engulfment",       "% of microglial volume"),
    "blue_CTB_engulfment_pct":        ("CTB Engulfment",          "% of microglial volume"),
    "red_VGluT2_engulfed_objects":    ("VGluT2 Engulfed Objects", "number of objects"),
    "blue_CTB_engulfed_objects":      ("CTB Engulfed Objects",    "number of objects"),
    "red_VGluT2_engulfed_volume_um3": ("VGluT2 Engulfed Volume",  "volume (\u00b5m\u00b3)"),
    "blue_CTB_engulfed_volume_um3":   ("CTB Engulfed Volume",     "volume (\u00b5m\u00b3)"),
    "red_VGluT2_density_per_1000um3": ("VGluT2 Object Density",    "objects / 1000 \u00b5m\u00b3"),
}

# Genotype colours for the primary comparison.
GENOTYPE_COLOR = {"A": (0.35, 0.35, 0.35), "B": (0.75, 0.12, 0.12)}

# Colourblind-friendly palette for individual animals.
ANIMAL_PALETTE = [
    (0.230, 0.299, 0.754), (0.706, 0.016, 0.150), (0.336, 0.706, 0.184),
    (0.800, 0.475, 0.655), (1.000, 0.501, 0.000), (0.400, 0.200, 0.600),
    (0.200, 0.630, 0.792), (0.900, 0.745, 0.000),
]


def load_data(csv_path, sex_map=None):
    """Load CSV; attach genotype and sex labels. Fail loudly if any animal is
    unmapped so nothing drops silently."""
    if sex_map is None:
        sex_map = SEX_MAP
    df = pd.read_csv(csv_path)
    df["genotype"] = df["group"].map(GROUP_LABELS)
    df["sex"] = df["replicate"].map(sex_map)
    unmapped = sorted(df.loc[df["sex"].isna(), "replicate"].unique())
    if unmapped:
        raise ValueError(f"Replicates missing from SEX_MAP: {unmapped}")
    df["subgroup"] = df["group"].astype(str) + "-" + df["sex"].astype(str)
    return df


def _beeswarm_x(values, center, width):
    values = np.asarray(values, float)
    n = len(values)
    x = np.full(n, center, float)
    if n <= 1:
        return x
    order = np.argsort(values)
    yr = (values.max() - values.min()) or 1.0
    placed_y, placed_x = [], []
    for idx in order:
        y = values[idx]
        close = [px for py, px in zip(placed_y, placed_x) if abs(py - y) < 0.035 * yr]
        if not close:
            xi = center
        else:
            step = width / 8.0
            xi = center
            for k in range(1, 40):
                for cand in (center + k * step, center - k * step):
                    if all(abs(cand - c) > step * 0.9 for c in close):
                        xi = cand
                        break
                else:
                    continue
                break
        placed_y.append(y)
        placed_x.append(xi)
        x[idx] = xi
    return np.clip(x, center - width / 2, center + width / 2)


CONNECTING_LINE_ALPHA = 0.25
CONNECTING_LINE_WIDTH = 0.5


def _draw_column(ax, sub, x_center, metric, animal_palette=True, base_color=None):
    """Draw one column: field swarm coloured by animal, thin lines connecting each
    field to its animal mean (as in the original SuperPlot), per-animal means, and
    a group mean +/- SD square. Returns animal means."""
    animals = sorted(sub["replicate"].unique())
    means = []
    swarm_w = 0.55
    mean_x_off = swarm_w * 0.55   # x offset of the animal-mean marker

    # Pre-compute positions so lines can be drawn first (behind the points).
    per_animal = []
    for j, animal in enumerate(animals):
        col = ANIMAL_PALETTE[j % len(ANIMAL_PALETTE)] if animal_palette else base_color
        vals = sub.loc[sub["replicate"] == animal, metric].values
        xs = _beeswarm_x(vals, x_center, swarm_w)
        m = float(np.mean(vals))
        means.append(m)
        per_animal.append((col, vals, xs, m))

    # 1) connecting lines (drawn first, low alpha, animal colour)
    for col, vals, xs, m in per_animal:
        mean_x = x_center + mean_x_off
        for xi, yi in zip(xs, vals):
            ax.plot([xi, mean_x], [yi, m], "-",
                    color=(col[0], col[1], col[2], CONNECTING_LINE_ALPHA),
                    linewidth=CONNECTING_LINE_WIDTH, zorder=1)

    # 2) field points
    for col, vals, xs, m in per_animal:
        ax.scatter(xs, vals, s=11, color=col, alpha=0.35, edgecolors="none", zorder=2)

    # 3) animal means
    for col, vals, xs, m in per_animal:
        ax.scatter(x_center + mean_x_off, m, s=60, color=col,
                   edgecolors="black", linewidths=0.8, zorder=4)

    means = np.array(means)
    gm = float(means.mean())
    gsd = float(means.std(ddof=1)) if len(means) > 1 else 0.0
    ax.errorbar(x_center + 0.72, gm, yerr=gsd, fmt="s", color="black",
                markerfacecolor="none", markersize=7, capsize=4,
                linewidth=1.2, zorder=5)
    return means


def _legend(fig, y=-0.06):
    handles = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markersize=6, alpha=0.5, label="imaging field"),
        Line2D([0], [0], color="gray", alpha=0.5, linewidth=0.8,
               label="field \u2192 animal mean"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markeredgecolor="k", markersize=9, label="animal mean"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="none",
               markeredgecolor="k", markersize=8, label="group mean \u00b1 SD"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, y))


def _finish_axis(ax, positions, xticklabels, y_min, y_max):
    yr = (y_max - y_min) or 1.0
    ax.set_xticks(positions)
    ax.set_xticklabels(xticklabels, fontsize=9)
    ax.set_xlim(positions.min() - 0.6, positions.max() + 1.1)
    ax.set_ylim(y_min - 0.08 * yr, y_max + 0.10 * yr)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")


# --------------------------------------------------------------------------
# Primary comparison: WT vs Opn4 KO (all animals pooled by genotype)
# --------------------------------------------------------------------------
def superplot_genotype(df, metric, outfile=None, dpi=300):
    title, ylab = METRIC_INFO[metric]
    fig, ax = plt.subplots(figsize=(6.0, 4.8))
    positions = np.array([0.0, 1.6])
    all_vals = df[metric].values
    y_min, y_max = np.nanmin(all_vals), np.nanmax(all_vals)
    labels = []
    for pos, g in zip(positions, ["A", "B"]):
        sub = df[df["group"] == g]
        _draw_column(ax, sub, pos, metric, animal_palette=True)
        labels.append(f"{GROUP_LABELS[g]}\n({sub['replicate'].nunique()} animals\n{len(sub)} fields)")
    _finish_axis(ax, positions, labels, y_min, y_max)
    ax.set_ylabel(ylab, fontsize=10)
    ax.set_title(title, fontsize=12, fontweight="bold")
    _legend(fig)
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
        fig.savefig(outfile.replace(".png", ".pdf"), bbox_inches="tight")
        fig.savefig(outfile.replace(".png", ".eps"), format="eps", bbox_inches="tight")
    return fig


# --------------------------------------------------------------------------
# Sex breakdown (inline in the review section)
# --------------------------------------------------------------------------
def superplot_split(df, metric, layout="separate", outfile=None, dpi=300):
    """layout='separate' -> Male | Female panels (WT vs KO within each).
       layout='four'     -> WT-M, WT-F, KO-M, KO-F in one axis."""
    title, ylab = METRIC_INFO[metric]
    sex_name = {"M": "Male", "F": "Female"}
    all_vals = df[metric].values
    y_min, y_max = np.nanmin(all_vals), np.nanmax(all_vals)

    if layout == "separate":
        fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.6), sharey=True)
        for ax, sx in zip(axes, ["M", "F"]):
            positions = np.array([0.0, 1.6])
            labels = []
            for pos, g in zip(positions, ["A", "B"]):
                sub = df[(df["group"] == g) & (df["sex"] == sx)]
                if len(sub):
                    _draw_column(ax, sub, pos, metric, animal_palette=True)
                labels.append(f"{GROUP_LABELS[g]}\n({sub['replicate'].nunique()} animals\n{len(sub)} fields)")
            _finish_axis(ax, positions, labels, y_min, y_max)
            ax.set_title(sex_name[sx], fontsize=11, fontweight="bold")
        axes[0].set_ylabel(ylab, fontsize=10)
        fig.suptitle(f"{title} \u2014 by genotype and sex", fontsize=12,
                     fontweight="bold", y=1.02)

    elif layout == "four":
        fig, ax = plt.subplots(figsize=(9.5, 4.8))
        cols = [("A", "M"), ("A", "F"), ("B", "M"), ("B", "F")]
        positions = np.arange(len(cols), dtype=float) * 1.0
        labels = []
        for pos, (g, sx) in zip(positions, cols):
            sub = df[(df["group"] == g) & (df["sex"] == sx)]
            if len(sub):
                _draw_column(ax, sub, pos, metric, animal_palette=True)
            labels.append(f"{GROUP_LABELS[g]}\n{sex_name[sx]}\n({sub['replicate'].nunique()} an. / {len(sub)} f.)")
        _finish_axis(ax, positions, labels, y_min, y_max)
        ax.set_ylabel(ylab, fontsize=10)
        ax.set_title(f"{title} \u2014 by genotype and sex", fontsize=12, fontweight="bold")
    else:
        raise ValueError("layout must be 'separate' or 'four'")

    _legend(fig)
    fig.tight_layout()
    if outfile:
        fig.savefig(outfile, dpi=dpi, bbox_inches="tight")
        fig.savefig(outfile.replace(".png", ".pdf"), bbox_inches="tight")
        fig.savefig(outfile.replace(".png", ".eps"), format="eps", bbox_inches="tight")
    return fig


if __name__ == "__main__":
    df = load_data("data/engulfment_results.csv")
    for metric in METRIC_INFO:
        superplot_genotype(df, metric, outfile=f"figures/{metric}_genotype.png")
        superplot_split(df, metric, "separate", outfile=f"figures/{metric}_split_separate.png")
        superplot_split(df, metric, "four", outfile=f"figures/{metric}_split_four.png")
    plt.close("all")
    print("figures written")
