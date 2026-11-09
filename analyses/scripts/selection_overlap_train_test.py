"""Overlap between CpGs selected on TRAIN (the saved sel_feats) and CpGs selected on TEST with the
same procedure, vs model performance.

The models select features with SelectKBest(f_regression) on the train samples, k tuned in CV.
Here the same selection (f_regression = top-k by |Pearson r| with expression, same k, same
candidate CpGs) is repeated on the TEST samples, and for each gene we report
    overlap = |selected_train ∩ selected_test| / k
against chance = k / n_candidates. Train selection is also recomputed here as a sanity check
that it reproduces the saved sel_feats.

Outputs (in out_dir): overlap_genes.tsv and selection_overlap_vs_performance.png/.pdf
"""
import argparse
import glob
import os

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--sel_dir", required=True)
p.add_argument("--perf_dir", required=True)
p.add_argument("--cand_dir", required=True)
p.add_argument("--h5", required=True)
p.add_argument("--positions", required=True, help="h5_cpg_positions.tsv (chr, pos in h5 column order)")
p.add_argument("--coldata", required=True, help="meth_colData.tsv (row order = h5 rows)")
p.add_argument("--counts_train", required=True)
p.add_argument("--counts_test", required=True)
p.add_argument("--prefix", default="hg38_v38_pc_50000_25_10")
p.add_argument("--out_dir", required=True)
p.add_argument("--shard", default="0/1", help="i/n: process every n-th gene starting at i; n>1 writes a shard table only")
p.add_argument("--plot_only", action="store_true", help="merge overlap_genes_shard*.tsv in out_dir and plot")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
shard_i, n_shards = map(int, a.shard.split("/"))


def abs_corr(X, y):
    Xc = X - X.mean(0)
    yc = y - y.mean()
    den = np.sqrt((Xc ** 2).sum(0) * (yc ** 2).sum())
    return np.abs(np.divide(Xc.T @ yc, den, out=np.zeros(X.shape[1]), where=den > 0))


def top_k(score, k):
    return set(np.argsort(-score, kind="stable")[:k])


def read_perf(TxID):
    f = os.path.join(a.perf_dir, f"{a.prefix}_3_{TxID}")
    s = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
    return s["cv_r2"], s["test_r2"]


def compute_genes():
    pos = pd.read_csv(a.positions, sep="\t")
    chr_cols = {c: (g.index.to_numpy(), g["pos"].to_numpy()) for c, g in pos.groupby("chr")}
    del pos

    samples = pd.read_csv(a.coldata, sep="\t", index_col=0).index.to_numpy()
    row_of = {s: i for i, s in enumerate(samples)}
    counts = {}
    for name, path in [("train", a.counts_train), ("test", a.counts_test)]:
        c = pd.read_csv(path, sep="\t", index_col=0)
        c = c[[s for s in c.columns if s in row_of]]
        counts[name] = (c, np.array([row_of[s] for s in c.columns]))

    rows = []
    with h5py.File(a.h5, "r") as h5:
        d = h5["M"]
        sel_files = sorted(glob.glob(os.path.join(a.sel_dir, "*")))[shard_i::n_shards]
        for i, f in enumerate(sel_files):
            TxID = os.path.basename(f).split("_")[-1]
            sel = {l.strip().replace(f"{TxID}_", "", 1) for l in open(f) if l.strip()}
            cand_file = os.path.join(a.cand_dir, f"{a.prefix}_{TxID}.bed")
            if not sel or not os.path.exists(cand_file) or TxID not in counts["test"][0].index:
                continue
            cand = pd.read_csv(cand_file, sep="\t", header=None, dtype=str)[1].drop_duplicates()
            chrom = cand.str.split("_").str[0]
            if chrom.nunique() != 1:
                continue
            cols_c, pos_c = chr_cols[chrom.iloc[0]]
            cpos = cand.str.split("_").str[1].astype(int).to_numpy()
            j = np.searchsorted(pos_c, cpos)
            ok = (j < len(pos_c)) & (pos_c[np.minimum(j, len(pos_c) - 1)] == cpos)
            cand, cols = cand[ok].to_numpy(), cols_c[j[ok]]
            order = np.argsort(cols)
            cand, cols = cand[order], cols[order]

            block = d[:, cols[0]:cols[-1] + 1][:, cols - cols[0]]      # samples x candidates
            k = len(sel)
            sel_idx = {n for n, c in enumerate(cand) if c in sel}
            res = {"TxID": TxID, "k": k, "n_cand": len(cand)}
            for name, (c, r) in counts.items():
                y = c.loc[TxID].to_numpy(dtype=np.float64)
                res[f"top_{name}"] = top_k(abs_corr(block[r].astype(np.float64), y), k)
            res["train_reproduced"] = len(res["top_train"] & sel_idx) / k
            res["overlap"] = len(res.pop("top_test") & sel_idx) / k
            res.pop("top_train")
            res["chance"] = k / len(cand)
            res["cv_r2"], res["test_r2"] = read_perf(TxID)
            rows.append(res)
            if i % 100 == 0:
                print(f"shard {shard_i}: {i}/{len(sel_files)}", flush=True)

    genes = pd.DataFrame(rows)
    genes[["overlap", "chance", "train_reproduced"]] *= 100
    return genes


if a.plot_only:
    genes = pd.concat([pd.read_csv(f, sep="\t")
                       for f in sorted(glob.glob(os.path.join(a.out_dir, "overlap_genes_shard*.tsv")))])
else:
    genes = compute_genes()
    if n_shards > 1:
        genes.to_csv(os.path.join(a.out_dir, f"overlap_genes_shard{shard_i}.tsv"), sep="\t", index=False)
        raise SystemExit
genes.to_csv(os.path.join(a.out_dir, "overlap_genes.tsv"), sep="\t", index=False)

# ---- plot ----
edges = [-np.inf, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, np.inf]
labels = ["<0", "0–0.1", "0.1–0.2", "0.2–0.3", "0.3–0.4", "0.4–0.5", "0.5–0.6", "0.6–0.7", "0.7–0.8", "0.8–1"]
genes["bin"] = pd.cut(genes["cv_r2"], edges, labels=labels, right=False)
g = genes.groupby("bin", observed=False)
summ = pd.DataFrame({"n_genes": g.size(),
                     "overlap_mean": g["overlap"].mean(), "overlap_se": g["overlap"].sem(),
                     "overlap_median": g["overlap"].median(),
                     "chance_mean": g["chance"].mean(), "median_k": g["k"].median()})
summ.to_csv(os.path.join(a.out_dir, "overlap_bins.tsv"), sep="\t")

fig, ax = plt.subplots(figsize=(9, 4.6))
x = np.arange(len(labels))
m, h = summ["overlap_mean"], 1.96 * summ["overlap_se"]
ax.errorbar(x, m, yerr=[m - (m - h).clip(lower=0), (m + h).clip(upper=100) - m], color="#2a78d6", lw=2,
            marker="o", ms=6, capsize=3, label="selected in train AND in test (same k)")
ax.plot(x, summ["chance_mean"], color="#eb6834", lw=2, marker="o", ms=6, label="expected by chance (k / n candidates)")
for xi, (m, v) in enumerate(zip(summ["overlap_mean"], summ["chance_mean"])):
    ax.annotate(f"{m:.0f}%", (xi, m), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8)
ax.set_xticks(x, [f"{l}\nn={n}\nk={int(kk)}" for l, n, kk in zip(labels, summ["n_genes"], summ["median_k"])], fontsize=7.5)
ax.set_xlabel("model CV R²   (n = genes, k = median no. selected CpGs)")
ax.set_ylabel("% of train-selected CpGs\nalso selected on test set")
ax.set_title("Overlap of selected CpGs: train vs test (SelectKBest f_regression)", fontsize=10)
ax.set_ylim(-3, 100)
ax.grid(axis="y", color="#e5e5e5", lw=0.8)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(a.out_dir, f"selection_overlap_vs_performance.{ext}"), dpi=200)
print(summ.round(2))
print("train selection reproduced (%):", genes["train_reproduced"].describe().round(2).to_dict())
