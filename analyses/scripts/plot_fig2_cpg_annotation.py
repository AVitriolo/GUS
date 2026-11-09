"""Report Fig. 2 panels: where the selected (functional) CpGs are and which motifs they carry.

Redraws the plots of plot_cpg_annotation.R and plot_motif_enrichment.R from
their TSVs, keeping their colours (dark / light grey for selected / candidate
CpGs, the blue ramp of the repeat pie, the light-yellow -> purple q-value ramp
of the motif bars), in the style of the Fig. 1 panels: DejaVu Sans, 8 pt
labels and ticks, 7 pt annotations, left and bottom axis lines only, no grid,
300 dpi. Sizes are on the same scale as Fig. 1B-E (10.1-in rows).

Outputs (in --out_dir):
  fig2_genomic_annotation.{png,pdf}  % of selected vs candidate CpGs per cCRE / gene-structure category
  fig2_repeat_classes.{png,pdf}      same per RepeatMasker class, plus the genomic region of the selected CpGs in repeats
  fig2_motif_enrichment.{png,pdf}    top --n non-redundant motifs by enrichment, coloured by -log10 q

Usage:
  python analyses/scripts/plot_fig2_cpg_annotation.py \
    --annotation_dir=results/annotation_w20 \
    --motif_tsv=results/gimme_w20/enrichment_2mb_full_kcapped_cv3r5_w20/top50_enriched_motifs_nonredundant_bar.tsv \
    --out_dir=results/plots --n=35
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch

# colours of plot_cpg_annotation.R / plot_motif_enrichment.R
SEL, CAND = "#4a4a48", "#c9c8c3"
INK, INK2 = "#0b0b0b", "#52514e"
PIE_RAMP = ["#0d366b", "#1c5cab", "#3987e5", "#6da7ec", "#9ec5f4", "#cde2fb"]   # promoter (dark) -> intergenic (light)
MOTIF_RAMP = LinearSegmentedColormap.from_list(
    "yellow_purple", ["#f3e1a9", "#e3c9a5", "#d1b0a0", "#b88f9a", "#a06e93", "#723287"])
Q_CAP = 200                                                                    # colour scale cap (-log10 q)

GENOMIC = ["Promoter (PLS)", "Proximal enhancer (pELS)", "Distal enhancer (dELS)", "Exon", "Intron", "Intergenic"]
REPEAT_CLASSES = ["No repeat", "LINE", "SINE", "LTR", "DNA transposon",
                  "Simple, low-complexity,\nsatellite", "Other"]


FG = "Selected"          # --fg_label: name of the foreground set in legends and notes
PREFIX = "fig2_"         # --prefix: output file prefix


def style(font="DejaVu Sans"):
    plt.rcParams.update({"font.family": font, "pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 8, "axes.titlesize": 9,
                         "axes.labelsize": 8, "xtick.labelsize": 8, "ytick.labelsize": 8,
                         "legend.fontsize": 8})


def despine(ax):
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "n.s."


def paired_bars(ax, labels, sel, cand, note=None, marks=None, note_top=False,
                legend=None, xlabel="% of CpGs"):
    """Selected (dark) above candidate (light) per category, ratio written at the end of each pair
    (followed by a significance mark when given)."""
    legend = legend or (f"{FG} CpGs", "Candidate CpGs")
    y = np.arange(len(labels))[::-1]
    h = 0.36
    ax.barh(y + h / 2, sel, height=h, color=SEL, edgecolor="white", linewidth=1)
    ax.barh(y - h / 2, cand, height=h, color=CAND, edgecolor="white", linewidth=1)
    xmax = max(max(sel), max(cand))
    marks = marks or [""] * len(labels)
    for yi, s, c, mk in zip(y, sel, cand, marks):
        ax.text(max(s, c) + 0.015 * xmax, yi, f"{s / c:.2f}× {mk}".strip(), ha="left", va="center", fontsize=7, color=INK2)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, xmax * (1.24 if any(marks) else 1.16))
    ax.set_xlabel(xlabel)
    ax.legend(handles=[Patch(color=SEL, label=legend[0]), Patch(color=CAND, label=legend[1])],
              frameon=False, ncol=2, loc="lower left", bbox_to_anchor=(0, 1.0), borderaxespad=0.2,
              handlelength=1.2, columnspacing=1.2)
    if note:   # lower right is empty when the last categories are short bars; else upper right
        ax.text(0.98, 0.97 if note_top else 0.03, note, transform=ax.transAxes, ha="right",
                va="top" if note_top else "bottom", fontsize=7, color=INK2)
    despine(ax)


def perm_panel(tf, out_dir, name, unit):
    """Observed share of selected points per category vs the mean of the random sets of candidate CpGs
    (same number, drawn per chromosome; permtest_cpg_localization.R): bars, folds and marks from one test."""
    t = pd.read_csv(tf, sep="\t").set_index("category").reindex(GENOMIC)
    n, nt = int(t.n_points.iloc[0]), int(t.ntimes.iloc[0])
    fig, ax = plt.subplots(figsize=(4.6, 2.9), layout="constrained")
    paired_bars(ax, GENOMIC, 100 * t.observed.to_numpy() / n, 100 * t.permuted_mean.to_numpy() / n,
                note=f"*** p < 0.001 ({nt:,} permutations of\n{n:,} random candidate CpGs)",
                marks=[stars(p) for p in t.p_empirical], note_top=True,
                legend=(f"{FG} {unit}", "Random candidate CpGs"), xlabel=f"% of {unit}")
    save(fig, out_dir, name)


def genomic_panel(ann_dir, out_dir):
    """Region-level panel (main) and single-CpG panel from the permutation tests when they exist;
    otherwise the descriptive % of selected vs candidate CpGs."""
    d = os.path.join(ann_dir, "localization_test")
    tests = [("regioner_permtest.tsv", f"{PREFIX}genomic_annotation", "regions"),
             ("regioner_permtest_cpgs.tsv", f"{PREFIX}genomic_annotation_cpgs", "CpGs")]
    if not any(os.path.exists(os.path.join(d, f)) for f, _, _ in tests):
        g = pd.read_csv(os.path.join(ann_dir, "genomic_annotation_counts.tsv"), sep="\t")
        pct = g.pivot(index="category", columns="set", values="pct").reindex(GENOMIC)
        fig, ax = plt.subplots(figsize=(4.6, 2.9), layout="constrained")
        paired_bars(ax, GENOMIC, pct["selected"].to_numpy(), pct["candidates"].to_numpy())
        save(fig, out_dir, f"{PREFIX}genomic_annotation")
        return
    for f, name, unit in tests:
        if os.path.exists(os.path.join(d, f)):
            perm_panel(os.path.join(d, f), out_dir, name, unit)


def repeat_class(c):
    c = c.rstrip("?")
    if c in ("LINE", "SINE", "LTR", "No repeat"):
        return c
    if c == "DNA":
        return "DNA transposon"
    if c in ("Simple_repeat", "Low_complexity", "Satellite"):
        return REPEAT_CLASSES[5]
    return "Other"


def repeat_panel(ann_dir, out_dir):
    r = pd.read_csv(os.path.join(ann_dir, "repeat_class_counts.tsv"), sep="\t")
    r = r[r.set.isin(["selected", "candidates"])].copy()
    r["cls"] = r.repClass.map(repeat_class)
    pct = r.groupby(["cls", "set"]).pct.sum().unstack().reindex(REPEAT_CLASSES)
    in_rep = r[r.cls != "No repeat"].groupby("set").pct.sum()
    note = f"In a repeat: {in_rep['selected']:.1f}% of {FG.lower()}\nvs {in_rep['candidates']:.1f}% of candidate CpGs"

    # genomic region of the selected CpGs that lie in a repeat
    cb = pd.read_csv(os.path.join(ann_dir, "selected_cpgs_category_by_repeat_class.tsv"), sep="\t")
    cb = cb.set_index("category").reindex(GENOMIC)
    n_rep = cb.drop(columns="No repeat").sum(axis=1)
    share = 100 * n_rep / n_rep.sum()

    # bars on the left; pie and its legend stacked on the right, the legend in its
    # own (empty) axes so the layout engine reserves room for it
    fig = plt.figure(figsize=(5.6, 3.0), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[3.2, 2.4], height_ratios=[1.25, 1])
    ax = fig.add_subplot(gs[:, 0])
    paired_bars(ax, REPEAT_CLASSES, pct["selected"].to_numpy(), pct["candidates"].to_numpy(), note=note)
    axp = fig.add_subplot(gs[0, 1])
    # clockwise from 12 o'clock, promoter first (as coord_polar(direction = -1) in the R plot)
    axp.pie(n_rep.to_numpy(), colors=PIE_RAMP, startangle=90, counterclock=False,
            wedgeprops={"edgecolor": "white", "linewidth": 0.6})
    axp.set_aspect("equal")
    axl = fig.add_subplot(gs[1, 1])
    axl.axis("off")
    leg = axl.legend(handles=[Patch(color=c, label=f"{cat}  {s:.1f}%") for c, cat, s in zip(PIE_RAMP, GENOMIC, share)],
                     title=f"{FG} CpGs in repeats\n(n = {int(n_rep.sum()):,})", title_fontsize=8,
                     frameon=False, loc="upper center", fontsize=7, handlelength=1.0, handleheight=0.9,
                     labelspacing=0.3)
    leg._legend_box.align = "left"
    save(fig, out_dir, f"{PREFIX}repeat_classes")


def motif_panel(motif_tsv, out_dir, n, xlabel):
    d = pd.read_csv(motif_tsv, sep="\t").sort_values("rank").head(n)
    mlq = -np.log10(np.maximum(d["q-value (BH)"].to_numpy(), 1e-300))
    lo = np.floor(mlq.min() / 10) * 10
    norm = Normalize(vmin=lo, vmax=Q_CAP, clip=True)
    enr = d["Enr. at 1% FPR"].to_numpy()

    # ~0.14 in per motif row keeps the 8 pt labels apart
    fig, ax = plt.subplots(figsize=(5.0, 0.14 * n + 0.75), layout="constrained")
    y = np.arange(n)[::-1]
    ax.barh(y, enr - 1, left=1, height=0.78, color=MOTIF_RAMP(norm(np.minimum(mlq, Q_CAP))),
            edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(d["name"].str.replace(r"\s+\[", " [", regex=True))
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlim(1, enr.max() * 1.02)
    ax.set_xlabel(xlabel)
    despine(ax)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=MOTIF_RAMP)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.02, aspect=18, shrink=0.45, anchor=(0, 0.9))
    ticks = [t for t in np.arange(lo, Q_CAP + 1, 50) if lo <= t <= Q_CAP]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"≥{Q_CAP}" if t == Q_CAP else f"{t:.0f}" for t in ticks])
    cbar.set_label("−log10 q (BH)")
    cbar.outline.set_visible(False)
    save(fig, out_dir, f"{PREFIX}motif_enrichment")


def save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(os.path.join(out_dir, f"{name}.{ext}"), dpi=300)
    plt.close(fig)
    print(f"Wrote {os.path.join(out_dir, name)}.{{png,pdf}}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--annotation_dir", required=True)
    p.add_argument("--motif_tsv", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--n", type=int, default=35, help="number of top motifs to show")
    p.add_argument("--fg_label", default="Selected", help="foreground set name, e.g. Stable")
    p.add_argument("--prefix", default="fig2_")
    p.add_argument("--font", default="DejaVu Sans", help="e.g. 'Nimbus Sans' (Helvetica metrics)")
    p.add_argument("--motif_xlabel", default="Enrichment at 1% FPR\n(selected vs matched unselected CpG windows)")
    args = p.parse_args()
    global FG, PREFIX
    FG, PREFIX = args.fg_label, args.prefix
    style(args.font)
    genomic_panel(args.annotation_dir, args.out_dir)
    repeat_panel(args.annotation_dir, args.out_dir)
    motif_panel(args.motif_tsv, args.out_dir, min(args.n, len(pd.read_csv(args.motif_tsv, sep="\t"))),
                args.motif_xlabel)


if __name__ == "__main__":
    main()
