#!/usr/bin/env python
# ============================================================
# combat_correction_rna.py — LEAKAGE-SAFE RNA batch correction
#
# ComBat (batch = dataset CPGEA/MCRPC, mod = None), whole matrix.
#   - parameters FIT on TRAIN samples only
#   - the SAME train-derived parameters applied to TEST
#
# Uses combat_fit_transform.py (validated against inmoose.pycombat_norm
# to < 1e-7 relative). Standard ComBat (sva / pycombat_norm) estimates and
# applies in one call and so cannot be used for train/test discipline.
#
# NOTE ON CONFOUNDING: dataset is aliased with metastasis (all MCRPC = M).
# With mod=None this removes the dataset axis AND the metastasis biology.
# N<->T (both CPGEA) is the clean contrast.
# ============================================================
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from combat_fit_transform import combat_fit, combat_apply

p = argparse.ArgumentParser()
p.add_argument("--counts_train", required=True,
               help="log2 normalized TRAIN counts (genes x samples)")
p.add_argument("--counts_test", required=True,
               help="log2 normalized TEST counts (genes x samples)")
p.add_argument("--out_train", required=True)
p.add_argument("--out_test", required=True)
p.add_argument("--out_params", default=None)
p.add_argument("--out_plot_ds", default=None)
p.add_argument("--out_plot_cond", default=None)
args = p.parse_args()

for f in (args.out_train, args.out_test):
    os.makedirs(os.path.dirname(f) or ".", exist_ok=True)

train = pd.read_csv(args.counts_train, sep="\t", index_col=0)
test = pd.read_csv(args.counts_test, sep="\t", index_col=0)
assert list(train.index) == list(test.index), \
    "Train and test gene lists differ - they must match (list is train-derived)."

print(f"Train: {train.shape[0]} genes x {train.shape[1]} samples")
print(f"Test : {test.shape[0]} genes x {test.shape[1]} samples")


def batch_of(s):
    return "MCRPC" if s.startswith("DTB") else "CPGEA"


def cond_of(s):
    if s.startswith("N"):
        return "Normal"
    if s.startswith("T"):
        return "Tumour"
    if s.startswith("DTB"):
        return "metastasis"
    return "unknown"


b_train = np.array([batch_of(s) for s in train.columns])
b_test = np.array([batch_of(s) for s in test.columns])
print(f"Train batches: {dict(pd.Series(b_train).value_counts())}")
print(f"Test  batches: {dict(pd.Series(b_test).value_counts())}")

missing = set(np.unique(b_test)) - set(np.unique(b_train))
if missing:
    raise ValueError(f"Batch(es) {missing} in TEST but not TRAIN - "
                     "cannot estimate their parameters from training data.")

# ---- FIT on TRAIN, APPLY to both ----
print("\nFitting ComBat on TRAIN (batch=dataset, mod=None) ...")
fit = combat_fit(train.values, b_train, covar_mod=None, par_prior=True)

print("Applying TRAIN-fitted parameters to train and test ...")
train_corr = combat_apply(train.values, b_train, fit, covar_mod=None)
test_corr = combat_apply(test.values, b_test, fit, covar_mod=None)

pd.DataFrame(train_corr, index=train.index, columns=train.columns) \
    .to_csv(args.out_train, sep="\t")
pd.DataFrame(test_corr, index=test.index, columns=test.columns) \
    .to_csv(args.out_test, sep="\t")

if args.out_params:
    np.savez(args.out_params,
             **{k: v for k, v in fit.items() if isinstance(v, np.ndarray)})

# ---- QC PCA: fit on TRAIN, project TEST (before vs after) ----
def pca_frame(mat_train, mat_test, stage, n_top=2000):
    rv = mat_train.var(axis=1)
    vg = rv.sort_values(ascending=False).head(min(n_top, mat_train.shape[0])).index
    Xtr = mat_train.loc[vg].T.values
    ctr = Xtr.mean(axis=0)
    Xtr_c = Xtr - ctr
    U, S, Vt = np.linalg.svd(Xtr_c, full_matrices=False)
    scores_tr = Xtr_c @ Vt.T[:, :2]
    Xte_c = mat_test.loc[vg].T.values - ctr
    scores_te = Xte_c @ Vt.T[:, :2]

    rows = []
    for i, s in enumerate(mat_train.columns):
        rows.append((scores_tr[i, 0], scores_tr[i, 1], cond_of(s), batch_of(s), "train", stage))
    for i, s in enumerate(mat_test.columns):
        rows.append((scores_te[i, 0], scores_te[i, 1], cond_of(s), batch_of(s), "test", stage))
    return pd.DataFrame(rows, columns=["PC1", "PC2", "condition", "dataset", "set", "stage"])


if args.out_plot_ds or args.out_plot_cond:
    tr_df = pd.DataFrame(train.values, index=train.index, columns=train.columns)
    te_df = pd.DataFrame(test.values, index=test.index, columns=test.columns)
    trc_df = pd.DataFrame(train_corr, index=train.index, columns=train.columns)
    tec_df = pd.DataFrame(test_corr, index=test.index, columns=test.columns)

    df = pd.concat([pca_frame(tr_df, te_df, "before"),
                    pca_frame(trc_df, tec_df, "after")], ignore_index=True)

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
        fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(path)
        plt.close(fig)

    if args.out_plot_ds:
        os.makedirs(os.path.dirname(args.out_plot_ds) or ".", exist_ok=True)
        draw("dataset", {"CPGEA": "#4575b4", "MCRPC": "#d73027"},
             "RNA ComBat - by DATASET (want: separation shrinks after)",
             args.out_plot_ds)
    if args.out_plot_cond:
        os.makedirs(os.path.dirname(args.out_plot_cond) or ".", exist_ok=True)
        draw("condition", {"Normal": "#2166ac", "Tumour": "#d6604d", "metastasis": "#1a9641"},
             "RNA ComBat - by CONDITION",
             args.out_plot_cond)

print(f"\nSaved:\n  {args.out_train}\n  {args.out_test}")
print("\nREMINDER: mod=None removes the dataset axis and, because dataset is")
print("confounded with metastasis, the metastasis biology with it.")