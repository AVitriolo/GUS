"""Test-set replication of model-selected CpGs vs model performance.

For every (gene, CpG) pair in the index built by test_replication_index.R (selected CpGs + an
equal-size random baseline of non-selected candidate CpGs from the same 2Mb window), computes the
partial correlation between methylation and expression adjusted for condition (Normal/Tumour/
metastasis), or with --no_adjust the plain Pearson correlation, separately in TRAIN and TEST samples. A CpG "replicates" if the test association is
significant (p < alpha) with the same sign as in train. The test samples took no part in the
model's feature selection, so this is an independent check of each selected CpG.

Outputs (in out_dir): pairs.tsv (per CpG-gene pair), genes.tsv (per gene: cv/test R2, % replicated),
and replication_vs_performance.png/.pdf (% replicated across CV R2 bins, selected vs random).
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
from scipy import stats
from statsmodels.stats.multitest import multipletests

p = argparse.ArgumentParser()
p.add_argument("--index", required=True)
p.add_argument("--h5", required=True)
p.add_argument("--coldata", required=True, help="meth_colData.tsv (row order = h5 rows)")
p.add_argument("--counts_train", required=True)
p.add_argument("--counts_test", required=True)
p.add_argument("--perf_dir", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--alpha", type=float, default=0.05)
p.add_argument("--no_adjust", action="store_true",
               help="plain Pearson correlation (no condition adjustment), as the unadjusted models")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

idx = pd.read_csv(a.index, sep="\t")
coldata = pd.read_csv(a.coldata, sep="\t", index_col=0)

# ---- read the needed h5 columns, one chunk at a time ----
cols = np.unique(idx["h5_col"].to_numpy())
with h5py.File(a.h5, "r") as f:
    d = f["M"]
    assert d.shape[0] == len(coldata)
    chunk = d.chunks[1]
    meth = np.empty((d.shape[0], len(cols)), dtype=np.float32)
    blocks = cols // chunk
    for b in np.unique(blocks):
        sel = np.where(blocks == b)[0]
        lo = b * chunk
        meth[:, sel] = d[:, lo:min(lo + chunk, d.shape[1])][:, cols[sel] - lo]
col_pos = pd.Series(np.arange(len(cols)), index=cols)
print(f"read {len(cols)} CpG columns", flush=True)

samples = coldata.index.to_numpy()
cond = coldata["condition"].to_numpy()
if a.no_adjust:
    cond = np.full(len(cond), "all")                     # one group: demeaning only = plain correlation


def load_counts(path):
    c = pd.read_csv(path, sep="\t", index_col=0)
    return c[[s for s in c.columns if s in set(samples)]]


sets = {}
for name, path in [("train", a.counts_train), ("test", a.counts_test)]:
    c = load_counts(path)
    rows = np.array([np.where(samples == s)[0][0] for s in c.columns])
    sets[name] = (c, rows, cond[rows])


def demean_by_group(x, g):
    """Residuals of x on condition dummies (= subtract per-condition means), column-wise."""
    x = x.astype(np.float64).copy()
    for lvl in np.unique(g):
        m = g == lvl
        x[m] -= x[m].mean(axis=0)
    return x


def partial_r(y, X, g):
    ry = demean_by_group(y[:, None], g)[:, 0]
    rX = demean_by_group(X, g)
    num = rX.T @ ry
    den = np.sqrt((rX ** 2).sum(0) * (ry ** 2).sum())
    r = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    df = len(y) - 1 - len(np.unique(g))
    t = r * np.sqrt(df / np.clip(1 - r ** 2, 1e-12, None))
    return r, 2 * stats.t.sf(np.abs(t), df)


out = []
for TxID, sub in idx.groupby("TxID", sort=False):
    res = sub[["TxID", "CpGID", "set"]].copy()
    for name, (c, rows, g) in sets.items():
        if TxID not in c.index:
            break
        y = c.loc[TxID].to_numpy(dtype=np.float64)
        X = meth[np.ix_(rows, col_pos[sub["h5_col"].to_numpy()].to_numpy())]
        res[f"r_{name}"], res[f"p_{name}"] = partial_r(y, X, g)
    else:
        out.append(res)
pairs = pd.concat(out, ignore_index=True)
pairs["replicated"] = (pairs["p_test"] < a.alpha) & (np.sign(pairs["r_test"]) == np.sign(pairs["r_train"]))
for s in ["selected", "random"]:
    m = pairs["set"] == s
    pairs.loc[m, "q_test"] = multipletests(pairs.loc[m, "p_test"], method="fdr_bh")[1]
pairs["replicated_fdr"] = (pairs["q_test"] < a.alpha) & (np.sign(pairs["r_test"]) == np.sign(pairs["r_train"]))
pairs.to_csv(os.path.join(a.out_dir, "pairs.tsv"), sep="\t", index=False)

# ---- per gene summary + performance ----
perf = {}
for f in glob.glob(os.path.join(a.perf_dir, "*")):
    s = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
    perf[os.path.basename(f).split("_")[-1]] = (s.get("cv_r2", np.nan), s.get("test_r2", np.nan))
genes = (pairs.groupby(["TxID", "set"])
         .agg(n=("replicated", "size"), pct_rep=("replicated", "mean"), pct_rep_fdr=("replicated_fdr", "mean"))
         .unstack("set"))
genes.columns = [f"{m}_{s}" for m, s in genes.columns]
genes["cv_r2"] = [perf.get(t, (np.nan, np.nan))[0] for t in genes.index]
genes["test_r2"] = [perf.get(t, (np.nan, np.nan))[1] for t in genes.index]
genes[[c for c in genes.columns if c.startswith("pct")]] *= 100
genes.to_csv(os.path.join(a.out_dir, "genes.tsv"), sep="\t")

# ---- plot: % of selected CpGs replicated in test, across CV R2 bins ----
edges = [-np.inf, 0, 0.1, 0.2, 0.3, 0.4, 0.5, np.inf]
labels = ["<0", "0–0.1", "0.1–0.2", "0.2–0.3", "0.3–0.4", "0.4–0.5", "≥0.5"]
pairs = pairs.merge(genes[["cv_r2"]], left_on="TxID", right_index=True)
pairs["bin"] = pd.cut(pairs["cv_r2"], edges, labels=labels, right=False)


def summarise(col):
    g = pairs.groupby(["bin", "set"], observed=False)[col]
    k, n = g.sum().unstack(), g.size().unstack()
    pct = 100 * k / n
    se = 100 * np.sqrt((k / n) * (1 - k / n) / n)
    return pct, se, n


ngenes = genes.assign(bin=pd.cut(genes["cv_r2"], edges, labels=labels, right=False)).groupby("bin", observed=False).size()
pairs.groupby(["bin", "set"], observed=False)["replicated"].agg(["sum", "size", "mean"]).to_csv(
    os.path.join(a.out_dir, "bins.tsv"), sep="\t")

colors = {"selected": "#2a78d6", "random": "#eb6834"}
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, col, title in [(axes[0], "replicated", f"test p < {a.alpha}, same sign as train"),
                       (axes[1], "replicated_fdr", f"test BH-FDR < {a.alpha}, same sign as train")]:
    pct, se, n = summarise(col)
    x = np.arange(len(labels))
    for s, lab in [("selected", "selected CpGs"), ("random", "random candidate CpGs (same window)")]:
        ax.errorbar(x, pct[s], yerr=1.96 * se[s], color=colors[s], lw=2, marker="o", ms=6,
                    capsize=3, label=lab)
    ax.set_xticks(x, [f"{l}\n(n={ngenes[l]})" for l in labels], fontsize=8)
    ax.set_xlabel("model CV R² (train only); n = genes")
    ax.set_title(title, fontsize=10)
    ax.grid(axis="y", color="#e5e5e5", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].set_ylabel("% CpGs replicated in test set")
axes[0].legend(frameon=False, fontsize=8, loc="upper left")
fig.suptitle("Test-set replication of CpG–expression association ("
             + ("no adjustment" if a.no_adjust else "adjusted for condition") + ")", fontsize=11)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(a.out_dir, f"replication_vs_performance.{ext}"), dpi=200)
print(pairs.groupby(["bin", "set"], observed=False)["replicated"].mean().unstack().mul(100).round(1))
print(ngenes)
