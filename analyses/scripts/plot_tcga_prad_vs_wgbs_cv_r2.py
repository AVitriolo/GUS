"""Compare cv_r2 distributions across four conditions:
  - WGBS complete (Prostate_GAS main pipeline, full CpG set)
  - WGBS CpGs in probes (Prostate_GAS, HM450-probe-restricted, same TSS window)
  - TCGA PRAD (full 533-sample cohort)
  - TCGA PRAD subsampled (342 samples, matching Prostate_GAS's total N)

Each group's distribution is independent (TCGA PRAD uses gene symbols, WGBS
uses transcript IDs - there is no shared identifier space to pair rows on), so
significance between the two within-pipeline comparisons (WGBS complete vs.
CpGs-in-probes; TCGA full vs. subsampled) is assessed with a two-sided
Mann-Whitney U test. The two WGBS boxes are drawn closer to each other than to
the two TCGA boxes, visually grouping each pipeline's pair.

Usage (2 Mb k-capped main run; TCGA uses the same +-1 Mb TSS window):
  python analyses/scripts/plot_tcga_prad_vs_wgbs_cv_r2.py \
    --performance_dir_tcga_prad=/home/ieo7293/TCGA_PRAD/results_prad_1mb/performance \
    --performance_dir_tcga_prad_subsampled=/home/ieo7293/TCGA_PRAD/results_prad_1mb_n342/performance \
    --performance_dir_wgbs_complete=results/performance_2mb_full_kcapped_cv3r5 \
    --wgbs_complete_prefix=hg38_v38_pc_50000_25_10_3 \
    --performance_dir_wgbs_probes=results/performance_2mb_hm450_kcapped_cv3r5 \
    --wgbs_probes_prefix=hg38_v38_pc_50000_25_10_3 \
    --output_path=results/plots/wgbs_2mb_kcapped_vs_tcga_prad_cv_r2_boxplot.pdf
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import mannwhitneyu

# Each pipeline's pair shares one hue from the validated categorical palette
# (slot 3 aqua for WGBS, slot 1 blue for TCGA); the second member of each pair
# is a 40%-white tint of that same hue, so the two pairs read as two groups
# rather than four unrelated series.
COLORS = {
    "WGBS complete": "#1baf7a",
    "WGBS CpGs in probes": "#76cfaf",
    "TCGA PRAD": "#2a78d6",
    "TCGA PRAD subsampled": "#7faee6",
}
ORDER = list(COLORS.keys())
# Tick labels carry the sample size; the per-box "n" annotation counts models.
LABELS = {
    "WGBS complete": "WGBS\nall CpGs\n(342 samples)",
    "WGBS CpGs in probes": "WGBS\nHM450 CpGs only\n(342 samples)",
    "TCGA PRAD": "TCGA PRAD\nHM450\n(533 samples)",
    "TCGA PRAD subsampled": "TCGA PRAD\nsubsampled\n(342 samples)",
}
# Tight gap within each pipeline's pair (0.8), wider gap between pipelines (1.4).
POSITIONS = [0.0, 0.8, 2.2, 3.0]
PAIRS = [("WGBS complete", "WGBS CpGs in probes"), ("TCGA PRAD", "TCGA PRAD subsampled")]


def read_cv_r2(path):
    with open(path) as fh:
        for line in fh:
            key, _, value = line.strip().partition("\t")
            if key == "cv_r2":
                try:
                    return float(value)
                except ValueError:
                    return None
    return None


def collect(pattern, group):
    rows = []
    for f in sorted(glob.glob(pattern)):
        r2 = read_cv_r2(f)
        if r2 is not None:
            rows.append((group, r2))
    return rows


def stars(p):
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def p_label(p):
    return f"p < 1e-300" if p == 0.0 else f"p = {p:.2e}"


def draw_bracket(ax, x1, x2, y, height, label):
    ax.plot([x1, x1, x2, x2], [y, y + height, y + height, y], color="black", linewidth=1)
    ax.text((x1 + x2) / 2, y + height, label, ha="center", va="bottom", fontsize=8)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--performance_dir_tcga_prad", required=True)
    p.add_argument("--performance_dir_tcga_prad_subsampled", required=True)
    p.add_argument("--performance_dir_wgbs_complete", required=True)
    p.add_argument("--wgbs_complete_prefix", required=True)
    p.add_argument("--performance_dir_wgbs_probes", required=True)
    p.add_argument("--wgbs_probes_prefix", required=True)
    p.add_argument("--output_path", required=True)
    p.add_argument("--no_titles", action="store_true", help="omit the title (report version: it goes in the caption)")
    args = p.parse_args()

    rows = []
    rows += collect(os.path.join(args.performance_dir_tcga_prad, "*_SVR.tsv"), "TCGA PRAD")
    rows += collect(os.path.join(args.performance_dir_tcga_prad_subsampled, "*_SVR.tsv"), "TCGA PRAD subsampled")
    rows += collect(os.path.join(args.performance_dir_wgbs_complete, f"{args.wgbs_complete_prefix}_*"), "WGBS complete")
    rows += collect(os.path.join(args.performance_dir_wgbs_probes, f"{args.wgbs_probes_prefix}_*"), "WGBS CpGs in probes")

    df = pd.DataFrame(rows, columns=["group", "cv_r2"])
    if df.empty:
        print("No cv_r2 values found in any of the four directories - skipping plot.")
        return

    data = {g: df.loc[df["group"] == g, "cv_r2"].to_numpy() for g in ORDER}
    counts = {g: len(data[g]) for g in ORDER}
    medians = {g: pd.Series(data[g]).median() for g in ORDER}

    # The tail of poorly-fit genes/TxIDs runs to very large negative cv_r2 and
    # would otherwise compress the informative region near 0 into a sliver.
    lo, hi = df["cv_r2"].quantile([0.01, 0.99])
    span = hi - lo
    ylim_bottom = lo - 0.1 * span
    n_below = int((df["cv_r2"] < ylim_bottom).sum())

    os.makedirs(os.path.dirname(args.output_path), exist_ok=True)

    # Drawn at roughly report size with the same 8/9 pt text as the coherence
    # figure (plot_wgbs_vs_tcga_feature_coherence.py), so both print alike at 100%.
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8})
    fig, ax = plt.subplots(figsize=(5.5, 3.9))
    bxp = ax.boxplot(
        [data[g] for g in ORDER], positions=POSITIONS, widths=0.6,
        patch_artist=True, showfliers=True,
        # markerfacecolor must be set: with the default ("none") and no edge, fliers are invisible.
        flierprops={"markersize": 2, "alpha": 0.2, "markeredgewidth": 0, "marker": "o",
                    "markerfacecolor": "#52514e"},
        boxprops={"linewidth": 1.2}, medianprops={"color": "black", "linewidth": 1.2},
        whiskerprops={"linewidth": 1.2}, capprops={"linewidth": 1.2},
    )
    for patch, g in zip(bxp["boxes"], ORDER):
        patch.set_facecolor(COLORS[g])

    # Headroom is measured from the highest mark actually drawn (whisker cap or
    # flier), not a fixed fraction of the data span, so the per-box annotation
    # row never collides with whichever box happens to reach highest.
    max_whisker_top = max(cap.get_ydata()[0] for cap in bxp["caps"])
    max_mark_top = max(max_whisker_top, max(data[g].max() for g in ORDER))
    ann_y = max_mark_top + 0.06 * span
    bracket_y = ann_y + 0.14 * span
    bracket_h = 0.04 * span
    ylim = (ylim_bottom, bracket_y + bracket_h + 0.06 * span)

    ax.axhline(0, color="grey", linewidth=1, linestyle="--", zorder=0)
    ax.set_ylim(ylim)
    ax.set_xlim(POSITIONS[0] - 0.6, POSITIONS[-1] + 0.6)

    ax.set_xticks(POSITIONS)
    ax.set_xticklabels([LABELS[g] for g in ORDER])
    ax.set_xlabel("")
    ax.set_ylabel("CV R²")

    for x, g in zip(POSITIONS, ORDER):
        # kept to ~13 characters per line so neighbouring boxes' labels don't touch
        ax.text(x, ann_y, f"{counts[g]:,} models\nmedian {medians[g]:.3f}",
                 ha="center", va="bottom", fontsize=7, color="dimgrey")
    for g1, g2 in PAIRS:
        x1 = POSITIONS[ORDER.index(g1)]
        x2 = POSITIONS[ORDER.index(g2)]
        _, pval = mannwhitneyu(data[g1], data[g2], alternative="two-sided")
        draw_bracket(ax, x1, x2, bracket_y, bracket_h, f"{stars(pval)}  ({p_label(pval)})")
        print(f"Mann-Whitney U, {g1} vs {g2}: {p_label(pval)}, {stars(pval)}")

    title = "CV R² distribution: WGBS vs. TCGA PRAD, full vs. CpG/sample-restricted"
    subtitle = f"y-axis clipped to 1st-99th percentile ({n_below} points below {ylim[0]:.2f} not shown)"
    if not args.no_titles:
        ax.set_title(f"{title}\n{subtitle}")

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.output_path, dpi=300)
    print(f"Wrote {args.output_path}")
    print(pd.Series(counts).reindex(ORDER).to_string())


if __name__ == "__main__":
    main()
