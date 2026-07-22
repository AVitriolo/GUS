#!/usr/bin/env python
# ============================================================
# plot_meth_combat_qc.py
#
# Regenerate the before/after QC PCA plots for the methylation
# ComBat correction WITHOUT re-running the correction.
#
# Reads the uncorrected and corrected H5s, picks the top-variable
# CpGs per chromosome (never realises the full matrix), fits PCA on
# TRAIN and projects TEST, and draws the 2x2 before/after panels.
#
# Usage:
#   python plot_meth_combat_qc.py \
#     --in_h5_before  resources/integrated_wgbs_hdf5/integrated_wgbs.h5 \
#     --in_h5_after   resources/integrated_wgbs_hdf5/integrated_wgbs_combat.h5 \
#     --in_rowranges  resources/integrated_wgbs_hdf5/meth_rowRanges.tsv \
#     --in_coldata    resources/integrated_wgbs_hdf5/meth_colData.tsv \
#     --out_plot_ds   results/plots/batch_correct/meth_combat_by_dataset.pdf \
#     --out_plot_cond results/plots/batch_correct/meth_combat_by_condition.pdf
# ============================================================
import argparse
import os

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument("--in_h5_before", required=True)
p.add_argument("--in_h5_after", required=True)
p.add_argument("--in_rowranges", required=True)
p.add_argument("--in_coldata", required=True)
p.add_argument("--out_plot_ds", required=True)
p.add_argument("--out_plot_cond", required=True)
p.add_argument("--n_top_pca", type=int, default=5000,
               help="global number of top-variable CpGs used for the PCA")
p.add_argument("--per_chrom_top", type=int, default=5000,
               help="candidates kept per chromosome before the global trim")
args = p.parse_args()

coldata = pd.read_csv(args.in_coldata, sep="\t", index_col=0)
samples = list(coldata.index)
batch = coldata["dataset"].to_numpy()
condition = coldata["condition"].to_numpy()
is_train = coldata["split"].to_numpy() == "train"

chrom_per_cpg = pd.read_csv(args.in_rowranges, sep="\t")["seqnames"].to_numpy()

print(f"Samples: {len(samples)} (train={is_train.sum()}, test={(~is_train).sum()})")

with h5py.File(args.in_h5_before, "r") as f:
    shp_b = f["M"].shape
with h5py.File(args.in_h5_after, "r") as f:
    shp_a = f["M"].shape
print(f"before H5 (h5py view): {shp_b}")
print(f"after  H5 (h5py view): {shp_a}")
assert shp_b == shp_a, "before/after H5 shapes differ"
assert shp_b[0] == len(samples), \
    f"H5 axis0={shp_b[0]} but colData has {len(samples)} samples"
n_cpg = shp_b[1]
assert len(chrom_per_cpg) == n_cpg, "rowRanges length != H5 CpG axis"

# ---- collect top-variable CpGs per chromosome from BOTH files ----
keep_b, keep_a, keep_v = [], [], []

with h5py.File(args.in_h5_before, "r") as fb, h5py.File(args.in_h5_after, "r") as fa:
    for chrom in pd.unique(chrom_per_cpg):
        idx = np.where(chrom_per_cpg == chrom)[0]
        lo, hi = idx.min(), idx.max() + 1
        assert np.array_equal(idx, np.arange(lo, hi)), f"{chrom}: not contiguous"

        Mb = fb["M"][:, lo:hi].T          # CpGs x samples
        v = Mb.var(axis=1)
        k = min(args.per_chrom_top, len(v))
        sel = np.argpartition(v, -k)[-k:]

        keep_b.append(Mb[sel, :])
        keep_a.append(fa["M"][:, lo:hi].T[sel, :])
        keep_v.append(v[sel])
        print(f"[{chrom}] kept {k} of {len(v)} CpGs", flush=True)
        del Mb

Mb = np.vstack(keep_b)
Ma = np.vstack(keep_a)
vv = np.concatenate(keep_v)
del keep_b, keep_a, keep_v

k = min(args.n_top_pca, Mb.shape[0])
gsel = np.argpartition(vv, -k)[-k:]
Mb, Ma = Mb[gsel, :], Ma[gsel, :]
print(f"PCA on {Mb.shape[0]} top-variable CpGs")


def frame(mat, stage):
    """PCA fit on TRAIN, TEST projected onto the same components."""
    Xtr = mat[:, is_train].T
    ctr = Xtr.mean(axis=0)
    Xtr_c = Xtr - ctr
    _, _, Vt = np.linalg.svd(Xtr_c, full_matrices=False)
    comp = Vt.T[:, :2]
    s_tr = Xtr_c @ comp
    s_te = (mat[:, ~is_train].T - ctr) @ comp

    rows = []
    for i, gi in enumerate(np.where(is_train)[0]):
        rows.append((s_tr[i, 0], s_tr[i, 1], condition[gi], batch[gi], "train", stage))
    for i, gi in enumerate(np.where(~is_train)[0]):
        rows.append((s_te[i, 0], s_te[i, 1], condition[gi], batch[gi], "test", stage))
    return pd.DataFrame(rows, columns=["PC1", "PC2", "condition", "dataset", "set", "stage"])


df = pd.concat([frame(Mb, "before"), frame(Ma, "after")], ignore_index=True)


def draw(color_col, palette, title, path):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for r, st in enumerate(["before", "after"]):
        for c, se in enumerate(["train", "test"]):
            ax = axes[c, r]
            sub = df[(df["stage"] == st) & (df["set"] == se)]
            for lvl, col in palette.items():
                m = sub[color_col] == lvl
                ax.scatter(sub[m]["PC1"], sub[m]["PC2"], s=12, alpha=.7,
                           color=col, label=lvl)
            ax.set_title(f"{se} - {st}", fontsize=10)
            ax.set_xlabel("PC1")
            ax.set_ylabel("PC2")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    print(f"  {path}")


draw("dataset", {"CPGEA": "#4575b4", "MCRPC": "#d73027"},
     "Methylation ComBat - by DATASET (want: separation shrinks after)",
     args.out_plot_ds)
draw("condition", {"Normal": "#2166ac", "Tumour": "#d6604d", "metastasis": "#1a9641"},
     "Methylation ComBat - by CONDITION",
     args.out_plot_cond)

print("\nDone.")
