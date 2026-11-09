#!/usr/bin/env python
# ============================================================
# combat_correction_meth.py — LEAKAGE-SAFE methylation batch correction
#
# ComBat (batch = dataset CPGEA/MCRPC, mod = None), PER CHROMOSOME:
#   - parameters FIT on TRAIN samples only
#   - the SAME train-derived parameters applied to TEST
#
# Per-chromosome because the full M matrix (~18.7M CpGs x 351) cannot be
# held dense. Batch effects are genome-wide, so per-chromosome empirical-
# Bayes priors are near-identical to a genome-wide fit (each chromosome
# still has >= tens of thousands of CpGs for stable estimation).
#
# Uses combat_fit_transform.py (validated against inmoose.pycombat_norm).
# Reads/writes only the M assay (corrected M does not back-transform to
# valid beta; the downstream pipeline uses M).
#
# NOTE: mod=None removes the dataset axis and, since dataset is confounded
# with metastasis, the metastasis biology with it. N<->T is the clean contrast.
# ============================================================
import argparse
import os

import h5py
import numpy as np
import pandas as pd

from combat_fit_transform import combat_fit, combat_apply

p = argparse.ArgumentParser()
p.add_argument("--in_h5", required=True)
p.add_argument("--in_rowranges", required=True,
               help="TSV with a 'seqnames' column, one row per CpG (chrom order)")
p.add_argument("--in_coldata", required=True,
               help="TSV: index=sample, columns include condition, dataset, split")
p.add_argument("--out_h5", required=True)
p.add_argument("--out_params", required=True)
p.add_argument("--out_plot_ds", default=None)
p.add_argument("--out_plot_cond", default=None)
p.add_argument("--n_top_pca", type=int, default=5000,
               help="top-variable CpGs used for the QC PCA")
args = p.parse_args()

os.makedirs(os.path.dirname(args.out_h5) or ".", exist_ok=True)

coldata = pd.read_csv(args.in_coldata, sep="\t", index_col=0)
samples = list(coldata.index)
batch = coldata["dataset"].to_numpy()
split = coldata["split"].to_numpy()
is_train = split == "train"

chrom_per_cpg = pd.read_csv(args.in_rowranges, sep="\t")["seqnames"].to_numpy()

print(f"Samples: {len(samples)} (train={is_train.sum()}, test={(~is_train).sum()})")
print(f"Train batches: {dict(pd.Series(batch[is_train]).value_counts())}")
print(f"Test  batches: {dict(pd.Series(batch[~is_train]).value_counts())}")

missing = set(np.unique(batch[~is_train])) - set(np.unique(batch[is_train]))
if missing:
    raise ValueError(f"Batch(es) {missing} in TEST but not TRAIN.")

with h5py.File(args.in_h5, "r") as fin:
    h5_shape = fin["M"].shape
# The integrated H5 is stored SAMPLES x CpGs (R column-major -> h5py sees it
# this way). ComBat needs CpGs x samples, so we slice the CpG (column) axis
# per chromosome and transpose each chunk.
assert h5_shape[0] == len(samples), \
    f"H5 axis0={h5_shape[0]} but colData has {len(samples)} samples"
n_cpg = h5_shape[1]
assert len(chrom_per_cpg) == n_cpg, \
    f"rowRanges length {len(chrom_per_cpg)} != H5 CpG axis {n_cpg}"

if os.path.exists(args.out_h5):
    os.remove(args.out_h5)

fit_store = {}
batches_ref = None

want_pca = bool(args.out_plot_ds or args.out_plot_cond)
per_chrom_top = args.n_top_pca
pca_before, pca_after, pca_var = [], [], []

with h5py.File(args.in_h5, "r") as fin, h5py.File(args.out_h5, "w") as fout:
    # Match the ORIGINAL layout exactly. The integration writes via R
    # (rhdf5) with dims = c(CpGs, samples); R is column-major, so h5py
    # sees the transpose: (samples, CpGs). We must reproduce that same
    # h5py layout so the R downstream reads the corrected H5 as CpGs x
    # samples, identically to the uncorrected one.
    dout = fout.create_dataset("M", shape=(len(samples), n_cpg), dtype="float64",
                               chunks=(len(samples), min(10000, n_cpg)))

    for chrom in pd.unique(chrom_per_cpg):
        idx = np.where(chrom_per_cpg == chrom)[0]
        lo, hi = idx.min(), idx.max() + 1
        assert np.array_equal(idx, np.arange(lo, hi)), \
            f"{chrom}: CpGs not contiguous in the H5."

        print(f"[{chrom}] {len(idx)} CpGs - reading ...", flush=True)
        # h5py view is (samples, CpGs); slice the CpG columns, transpose
        # to CpGs x samples for ComBat.
        M = fin["M"][:, lo:hi].T

        fit = combat_fit(M[:, is_train], batch[is_train],
                         covar_mod=None, par_prior=True)
        M_corr = combat_apply(M, batch, fit, covar_mod=None)
        # write back in the original (samples, CpGs) h5py layout
        dout[:, lo:hi] = M_corr.T

        batches_ref = fit["batches"]
        fit_store[f"{chrom}_gamma"] = fit["gamma_star"]
        fit_store[f"{chrom}_delta"] = fit["delta_star"]

        gb = abs(M[:, batch == "CPGEA"].mean() - M[:, batch == "MCRPC"].mean())
        ga = abs(M_corr[:, batch == "CPGEA"].mean() - M_corr[:, batch == "MCRPC"].mean())
        print(f"[{chrom}] batch gap {gb:.4f} -> {ga:.4f}", flush=True)

        if want_pca:
            v = M.var(axis=1)
            k = min(per_chrom_top, len(v))
            sel = np.argpartition(v, -k)[-k:]
            pca_before.append(M[sel, :])
            pca_after.append(M_corr[sel, :])
            pca_var.append(v[sel])

        del M, M_corr

np.savez(args.out_params, batches=batches_ref, **fit_store)
print(f"\nSaved:\n  {args.out_h5}\n  {args.out_params}")

# ---- QC PCA: fit on TRAIN, project TEST (before vs after) ----
if want_pca:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    Mb = np.vstack(pca_before)
    Ma = np.vstack(pca_after)
    vv = np.concatenate(pca_var)

    k = min(args.n_top_pca, Mb.shape[0])
    gsel = np.argpartition(vv, -k)[-k:]
    Mb, Ma = Mb[gsel, :], Ma[gsel, :]

    condition = coldata["condition"].to_numpy()
    sample_arr = np.array(samples)

    def frame(mat, stage):
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
                ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        axes[0, 0].legend(fontsize=8)
        fig.suptitle(title); fig.tight_layout()
        fig.savefig(path); plt.close(fig)

    if args.out_plot_ds:
        os.makedirs(os.path.dirname(args.out_plot_ds) or ".", exist_ok=True)
        draw("dataset", {"CPGEA": "#4575b4", "MCRPC": "#d73027"},
             "Methylation ComBat - by DATASET (want: separation shrinks after)",
             args.out_plot_ds)
        print(f"  {args.out_plot_ds}")
    if args.out_plot_cond:
        os.makedirs(os.path.dirname(args.out_plot_cond) or ".", exist_ok=True)
        draw("condition", {"Normal": "#2166ac", "Tumour": "#d6604d", "metastasis": "#1a9641"},
             "Methylation ComBat - by CONDITION",
             args.out_plot_cond)
        print(f"  {args.out_plot_cond}")

print("\nREMINDER: mod=None removes the dataset axis and the confounded")
print("metastasis biology with it. N<->T (both CPGEA) is the clean contrast.")