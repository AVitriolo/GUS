"""Feature-selection coherence between the WGBS models and TCGA PRAD.

For every WGBS transcript model, the selected CpGs that sit on an Illumina HM450
probe are mapped to probe IDs and compared with the probes selected by the TCGA
PRAD model of the same gene (transcript -> gene via the GENCODE conversion
table). Coherence = % of our probe-mapped selected CpGs also selected by TCGA.

Chance expectation, per transcript: if TCGA had picked its features at random
among its own candidate probes, a probe-mapped CpG of ours would be selected
with probability (TCGA selected / TCGA candidates) when it is a TCGA candidate,
and never otherwise, so
    expected = (% of our probe-mapped CpGs that are TCGA candidates)
               * (TCGA selected / TCGA candidates).
Observed vs. expected is compared per transcript with a two-sided Wilcoxon
signed-rank test.

Transcripts are "comparable" when they have >= 1 probe-mapped selected CpG and
their gene has a TCGA model; the others are counted as 0% in the all-TxIDs panel.

Usage (2 Mb k-capped main run vs. TCGA PRAD, same +-1 Mb TSS window):
  python analyses/scripts/plot_wgbs_vs_tcga_feature_coherence.py \
    --sel_feats_dir_wgbs=results/sel_feats_2mb_full_kcapped_cv3r5 \
    --performance_dir_wgbs=results/performance_2mb_full_kcapped_cv3r5 \
    --sel_feats_dir_tcga=/home/ieo7293/TCGA_PRAD/results_prad_1mb/sel_feats \
    --data_dir_tcga=/home/ieo7293/TCGA_PRAD/data_prad_1mb \
    --manifest=resources/illumina_probes/HM450.hg38.manifest.tsv.gz \
    --conversion_table=resources/conversion_table_gencodeV38.tsv \
    --output_path=results/plots/wgbs_2mb_kcapped_vs_tcga_prad_feature_coherence.pdf \
    --output_table=results/plots/wgbs_2mb_kcapped_vs_tcga_prad_feature_coherence.tsv
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
from scipy.stats import wilcoxon

# Validated categorical slot 3 (aqua) = WGBS, as in plot_tcga_prad_vs_wgbs_cv_r2.py;
# the chance expectation is a reference, not a series, so it stays in muted ink.
WGBS = "#1baf7a"
CHANCE = "#898781"
INK = "#0b0b0b"
NOT_COMPARABLE = ["#c3c2b7", "#898781"]
CV_BINS = [-np.inf, 0, 0.1, 0.2, 0.3, np.inf]
CV_BIN_LABELS = ["≤ 0", "0–0.1", "0.1–0.2", "0.2–0.3", "> 0.3"]


def read_lines(path):
    with open(path) as fh:
        return [line.strip() for line in fh if line.strip()]


def read_cv_r2(path):
    with open(path) as fh:
        for line in fh:
            key, _, value = line.strip().partition("\t")
            if key == "cv_r2":
                try:
                    return float(value)
                except ValueError:
                    return np.nan
    return np.nan


def pct_p_label(p):
    return "p < 1e-300" if p == 0.0 else f"p = {p:.1e}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sel_feats_dir_wgbs", required=True)
    p.add_argument("--performance_dir_wgbs", required=True)
    p.add_argument("--sel_feats_dir_tcga", required=True)
    p.add_argument("--data_dir_tcga", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--conversion_table", required=True)
    p.add_argument("--output_path", required=True)
    p.add_argument("--output_table", required=True)
    p.add_argument("--no_titles", action="store_true", help="omit panel titles (report version: they go in the caption)")
    p.add_argument("--match", choices=["target", "footprint"], default="target",
                   help="target: a CpG is a probe only at the probe's target CpG; footprint: any CpG inside the probe's "
                        "50-nt footprint maps to that probe (see probe_footprint.py)")
    p.add_argument("--cv_above_mean", action="store_true",
                   help="keep only WGBS models with CV R² above the mean CV R² of all models")
    args = p.parse_args()

    # Manifest CpG_beg is 0-based; the pipeline's CpGIDs are 1-based (+1, see
    # cluster_probe_coherence.R: cg00000029 = chr16 53434200).
    man = pd.read_csv(args.manifest, sep="\t", usecols=["CpG_chrm", "CpG_beg", "Probe_ID"]).dropna()
    pos2probe = {f"{c}_{int(b) + 1}": pid for c, b, pid in zip(man.CpG_chrm, man.CpG_beg, man.Probe_ID)}
    print(f"Loaded {len(pos2probe)} HM450 probe positions")
    if args.match == "footprint":
        from probe_footprint import load_footprints, cpgs_to_probes
        footprints = load_footprints(args.manifest)

    conv = pd.read_csv(args.conversion_table, sep=" ")
    tx2gene = dict(zip(conv.transcript_id.str.replace(r"\.\d+$", "", regex=True), conv.gene_name))

    tcga_sel, tcga_cand = {}, {}
    for f in glob.glob(os.path.join(args.sel_feats_dir_tcga, "*_SVR.tsv")):
        tcga_sel[os.path.basename(f)[:-len("_SVR.tsv")]] = set(read_lines(f))
    for f in glob.glob(os.path.join(args.data_dir_tcga, "*_train.tsv")):
        with open(f) as fh:
            header = fh.readline().rstrip("\n").split("\t")
        tcga_cand[os.path.basename(f)[:-len("_train.tsv")]] = {x for x in header if x.startswith("cg")}
    print(f"TCGA PRAD: {len(tcga_sel)} genes with selected features, {len(tcga_cand)} with candidates")

    rows = []
    for f in sorted(glob.glob(os.path.join(args.sel_feats_dir_wgbs, "*"))):
        tx = f.rsplit("_", 1)[1]
        feats = read_lines(f)
        # features are <TxID>_<chr>_<pos>
        if args.match == "footprint":
            probes = {pr for hit in cpgs_to_probes((x.split("_", 1)[1] for x in feats), footprints).values() for pr in hit}
        else:
            probes = {pos2probe[x.split("_", 1)[1]] for x in feats if x.split("_", 1)[1] in pos2probe}
        gene = tx2gene.get(tx)
        cv_r2 = read_cv_r2(os.path.join(args.performance_dir_wgbs, os.path.basename(f)))
        observed = expected = np.nan
        if not probes:
            status = "no_probe"
        elif gene not in tcga_sel or not tcga_cand.get(gene):
            status = "gene_not_in_tcga"
        else:
            status = "comparable"
            observed = 100 * len(probes & tcga_sel[gene]) / len(probes)
            expected = (100 * len(probes & tcga_cand[gene]) / len(probes)
                        * len(tcga_sel[gene]) / len(tcga_cand[gene]))
        rows.append((tx, gene, len(feats), len(probes), status, cv_r2, observed, expected))

    df = pd.DataFrame(rows, columns=["TxID", "gene", "n_selected", "n_probe_mapped", "status",
                                     "cv_r2", "coherence_pct", "expected_pct"])
    if args.cv_above_mean:
        mean_cv = df.cv_r2.mean()
        df = df[df.cv_r2 > mean_cv].reset_index(drop=True)
        print(f"Kept {len(df)} models with CV R² > mean ({mean_cv:.4f})")
    os.makedirs(os.path.dirname(args.output_table), exist_ok=True)
    df.to_csv(args.output_table, sep="\t", index=False)

    comp = df[df.status == "comparable"].copy()
    all_obs = df.coherence_pct.fillna(0)
    all_exp = df.expected_pct.fillna(0)
    _, pval = wilcoxon(comp.coherence_pct, comp.expected_pct, alternative="two-sided")
    comp["cv_bin"] = pd.cut(comp.cv_r2, CV_BINS, labels=CV_BIN_LABELS)
    by_bin = comp.groupby("cv_bin", observed=False).agg(
        n=("TxID", "size"), observed=("coherence_pct", "mean"), expected=("expected_pct", "mean"))
    by_bin = by_bin[by_bin.n > 0]          # with --cv_above_mean the low bins are empty

    # Same 8-9 pt text as plot_tcga_prad_vs_wgbs_cv_r2.py (5.5 in wide) and
    # plot_cv_vs_test_r2.py (4.6 in wide). The figure is as wide as those two
    # side by side, so scaling the row and this figure to the page width uses
    # one factor and the text stays the same size; the columns line up under
    # them. The comparable-TxID histogram is the main panel, the two bar
    # charts sit smaller on its right.
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8})
    fig = plt.figure(figsize=(5.5 + 4.6, 4.4), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[5.5, 4.6])
    ax_hist = fig.add_subplot(gs[:, 0])
    ax_comp = fig.add_subplot(gs[0, 1])
    ax_perf = fig.add_subplot(gs[1, 1])

    # coherence histogram, comparable TxIDs; observed mean vs. chance mean
    bins = np.linspace(0, 100, 21)
    ax_hist.hist(comp.coherence_pct, bins=bins, color=WGBS, edgecolor="white", linewidth=1)
    # the mean lines stop at 85% of the height so they don't run through the
    # legend; the 0% bin stops at ~70%, so the top-left corner is free for it
    ax_hist.axvline(comp.coherence_pct.mean(), ymax=0.85, color=INK, linestyle="--", linewidth=1.2,
                    label=f"observed mean {comp.coherence_pct.mean():.1f}%")
    ax_hist.axvline(comp.expected_pct.mean(), ymax=0.85, color=CHANCE, linestyle=":", linewidth=1.6,
                    label=f"chance mean {comp.expected_pct.mean():.1f}%")
    ax_hist.legend(frameon=False, loc="upper left")
    ax_hist.set_xlim(-2, 102)
    ax_hist.set_xlabel("% of probe-mapped selected CpGs\nalso selected by TCGA (same gene)")
    ax_hist.set_ylabel("Transcript models")
    if not args.no_titles:
        ax_hist.set_title(f"Comparable transcripts (n = {len(comp):,})\n"
                          f"observed vs. chance: Wilcoxon {pct_p_label(pval)}")

    # which TxIDs are comparable - horizontal so the category names stay legible
    counts = df.status.value_counts().reindex(["comparable", "no_probe", "gene_not_in_tcga"], fill_value=0)
    labels = ["Comparable", "No selected CpG\non a 450K probe", "Gene not in\nTCGA models"]
    y = np.arange(len(counts))[::-1]
    ax_comp.barh(y, counts.values, color=[WGBS] + NOT_COMPARABLE, height=0.65, edgecolor="white", linewidth=1)
    for yi, v in zip(y, counts.values):
        ax_comp.text(v + 0.02 * counts.max(), yi, f"{v:,}", ha="left", va="center")
    ax_comp.set_yticks(y)
    ax_comp.set_yticklabels(labels)
    ax_comp.set_xlim(0, counts.max() * 1.35)
    ax_comp.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v / 1000:.0f}k" if v else "0"))
    ax_comp.set_xlabel("Transcript models")
    if not args.no_titles:
        ax_comp.set_title(f"Comparable to TCGA? (n = {len(df):,})")

    # coherence by model performance (comparable TxIDs)
    x = np.arange(len(by_bin))
    w = 0.38
    ax_perf.bar(x - w / 2 - 0.01, by_bin.observed, width=w, color=WGBS, label="observed")
    ax_perf.bar(x + w / 2 + 0.01, by_bin.expected, width=w, color=CHANCE, label="chance")
    for xi, (_, r) in zip(x, by_bin.iterrows()):
        ax_perf.text(xi - w / 2 - 0.01, r.observed, f"{r.observed:.0f}", ha="center", va="bottom", fontsize=7)
        ax_perf.text(xi + w / 2 + 0.01, r.expected, f"{r.expected:.0f}", ha="center", va="bottom", fontsize=7,
                     color="#52514e")
    ax_perf.set_xticks(x)
    ax_perf.set_xticklabels([f"{lab}\n{n:,}" for lab, n in zip(by_bin.index, by_bin.n)], fontsize=7.5)
    ax_perf.set_xlabel("WGBS model CV R² (models per bin)")
    ax_perf.set_ylabel("Mean coherence (%)")
    ax_perf.set_ylim(0, 100)
    ax_perf.legend(frameon=False, loc="upper left", ncol=2)
    if not args.no_titles:
        ax_perf.set_title("By model performance")

    for ax in (ax_hist, ax_comp, ax_perf):
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
    fig.savefig(args.output_path, dpi=300)
    print(f"Wrote {args.output_path}")

    print(df.status.value_counts().to_string())
    print(f"All TxIDs:        observed mean {all_obs.mean():.1f}%, chance mean {all_exp.mean():.1f}%")
    print(f"Comparable TxIDs: observed mean {comp.coherence_pct.mean():.1f}% "
          f"(median {comp.coherence_pct.median():.1f}), chance mean {comp.expected_pct.mean():.1f}% "
          f"(median {comp.expected_pct.median():.1f}); Wilcoxon {pct_p_label(pval)}")
    print(by_bin.round(1).to_string())


if __name__ == "__main__":
    main()
