"""Feature-selector benchmark, step 1: SelectKBest and elastic net, full train + half-samples.

For each benchmark gene (candidate CpGs = the gene's 2Mb set, train samples only, no covariate
adjustment):
  SKB   SelectKBest(f_regression) with the gene's tuned k (= saved selected CpGs), as in the SVR models
  ENET  ElasticNet tuned exactly as model_selection_enet_repeated.py (StandardScaler, standardised y,
        l1_ratio x alpha grid, RepeatedKFold 3x5 with the pipeline seed); selected = non-zero coefficients
Each selector is run on the full train set ("full") and, with the full-data k / (alpha, l1_ratio)
fixed, on n_pairs complementary pairs of condition-stratified half-samples.

Also writes the shared inputs for benchmark_susie.R / benchmark_evaluate.py / benchmark_summary.R:
h5 columns of the candidates, train/test expression, the half-samples, and co-methylated regions
(consecutive CpGs <= max_gap bp apart with adjacent Pearson r >= corr in train; as stability_regions_corr.py).
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.linear_model import ElasticNet, enet_path
from sklearn.metrics import r2_score
from sklearn.model_selection import RepeatedKFold
from sklearn.preprocessing import StandardScaler

p = argparse.ArgumentParser()
p.add_argument("--genes", required=True, help="benchmark gene list (TxIDs)")
p.add_argument("--sel_dir", required=True)
p.add_argument("--cand_dir", required=True)
p.add_argument("--h5", required=True)
p.add_argument("--positions", required=True)
p.add_argument("--coldata", required=True)
p.add_argument("--counts_train", required=True)
p.add_argument("--counts_test", required=True)
p.add_argument("--prefix", default="hg38_v38_pc_50000_25_10")
p.add_argument("--n_pairs", type=int, default=50)
p.add_argument("--max_gap", type=int, default=300)
p.add_argument("--corr", type=float, default=0.5)
p.add_argument("--skip_enet", action="store_true", help="SelectKBest only (elastic net tuning is ~30 min/gene)")
p.add_argument("--permute", action="store_true",
               help="null run: shuffle train expression within condition (use with --n_pairs 0); keeps the "
                    "between-condition expression differences, i.e. tests what is selected via tissue state alone")
p.add_argument("--permute_all", action="store_true",
               help="null run: shuffle train expression across all samples (no association left)")
p.add_argument("--shard", default="0/1")
p.add_argument("--out_dir", required=True)
a = p.parse_args()
shard_i, n_shards = map(int, a.shard.split("/"))

SEED = 48047510                                   # model_selection_enet_repeated.py
L1_RATIOS = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
ALPHAS = np.logspace(1, -2.5, 15)

pos = pd.read_csv(a.positions, sep="\t")
chr_cols = {c: (g.index.to_numpy(), g["pos"].to_numpy()) for c, g in pos.groupby("chr")}
del pos
coldata = pd.read_csv(a.coldata, sep="\t", index_col=0)
row_of = {s: i for i, s in enumerate(coldata.index)}
counts = {}
for name, path in [("train", a.counts_train), ("test", a.counts_test)]:
    c = pd.read_csv(path, sep="\t", index_col=0)
    c = c[[s for s in c.columns if s in row_of]]
    counts[name] = (c, np.array([row_of[s] for s in c.columns]))
cond = coldata["condition"].to_numpy()[counts["train"][1]]


def enet_tune(X, y):
    """Pipeline's hyperparameter search; returns (l1_ratio, alpha, cv_r2)."""
    scores = []
    for tr, va in RepeatedKFold(n_splits=3, n_repeats=5, random_state=SEED).split(X):
        xs = StandardScaler().fit(X[tr])
        Xtr, Xva = xs.transform(X[tr]), xs.transform(X[va])
        mu, sd = y[tr].mean(), y[tr].std()
        out = np.empty((len(L1_RATIOS), len(ALPHAS)))
        for i, l1 in enumerate(L1_RATIOS):
            _, coefs, _ = enet_path(Xtr, (y[tr] - mu) / sd, l1_ratio=l1, alphas=ALPHAS, max_iter=5000)
            pred = Xva @ coefs * sd + mu
            out[i] = [r2_score(y[va], pred[:, j]) for j in range(len(ALPHAS))]
        scores.append(out)
    scores = np.mean(scores, axis=0)
    i, j = np.unravel_index(np.argmax(scores), scores.shape)
    return L1_RATIOS[i], ALPHAS[j], scores[i, j]


def enet_select(X, y, l1, alpha):
    Xs = StandardScaler().fit_transform(X)
    ys = (y - y.mean()) / y.std()
    return np.flatnonzero(ElasticNet(alpha=alpha, l1_ratio=l1, max_iter=5000).fit(Xs, ys).coef_)


def skb_select(X, y, k):
    return SelectKBest(f_regression, k=k).fit(X, y).get_support(indices=True)


def halves(rng):
    for _ in range(a.n_pairs):
        h = np.zeros(len(cond), dtype=bool)
        for lvl in np.unique(cond):
            idx = rng.permutation(np.flatnonzero(cond == lvl))
            h[idx[: len(idx) // 2]] = True
        yield np.flatnonzero(h)
        yield np.flatnonzero(~h)


def join(v):
    return ",".join(map(str, v))


genes = [g.strip() for g in open(a.genes) if g.strip()][shard_i::n_shards]
fh = {n: open(os.path.join(a.out_dir, f"bench_{n}_shard{shard_i}.tsv"), "w")
      for n in ["selections", "genes", "inputs", "halves"]}
with h5py.File(a.h5, "r") as h5:
    d = h5["M"]
    for gi, TxID in enumerate(genes):
        sel = [l.strip() for l in open(os.path.join(a.sel_dir, f"{a.prefix}_3_{TxID}")) if l.strip()]
        k = len(sel)
        cand = pd.read_csv(os.path.join(a.cand_dir, f"{a.prefix}_{TxID}.bed"), sep="\t", header=None,
                           dtype=str)[1].drop_duplicates()
        cols_c, pos_c = chr_cols[cand.str.split("_").str[0].iloc[0]]
        cpos = cand.str.split("_").str[1].astype(int).to_numpy()
        j = np.searchsorted(pos_c, cpos)
        ok = (j < len(pos_c)) & (pos_c[np.minimum(j, len(pos_c) - 1)] == cpos)
        cpos, cols = cpos[ok], cols_c[j[ok]]
        order = np.argsort(cpos)
        cpos, cols = cpos[order], cols[order]

        block = d[:, cols.min():cols.max() + 1][:, cols - cols.min()]
        Xtr = block[counts["train"][1]].astype(np.float64)
        ytr = counts["train"][0].loc[TxID].to_numpy(dtype=np.float64)
        yte = counts["test"][0].loc[TxID].to_numpy(dtype=np.float64)
        if a.permute_all:
            ytr = ytr[np.random.default_rng([SEED, int(TxID[4:]), 2]).permutation(len(ytr))]
        if a.permute:
            prng = np.random.default_rng([SEED, int(TxID[4:]), 1])
            for lvl in np.unique(cond):
                m = np.flatnonzero(cond == lvl)
                ytr[m] = ytr[prng.permutation(m)]

        # co-methylated regions on train methylation
        Xc = Xtr - Xtr.mean(0)
        sd = np.sqrt((Xc ** 2).sum(0))
        with np.errstate(invalid="ignore", divide="ignore"):
            r_adj = (Xc[:, :-1] * Xc[:, 1:]).sum(0) / (sd[:-1] * sd[1:])
        region = np.cumsum(np.concatenate([[True], (np.diff(cpos) > a.max_gap) | ~(r_adj >= a.corr)]))

        fh["inputs"].write(f"{TxID}\tcols\t{join(cols)}\n{TxID}\tpos\t{join(cpos)}\n{TxID}\tregion\t{join(region)}\n"
                           f"{TxID}\ty_train\t{join(ytr)}\n{TxID}\ty_test\t{join(yte)}\n")

        l1, alpha, cv_r2 = (np.nan, np.nan, np.nan) if a.skip_enet else enet_tune(Xtr, ytr)
        fh["genes"].write(f"{TxID}\t{len(cols)}\t{k}\t{l1}\t{alpha}\t{cv_r2}\n")
        fh["selections"].write(f"{TxID}\tSKB\tfull\t{join(skb_select(Xtr, ytr, k))}\n")
        if not a.skip_enet:
            fh["selections"].write(f"{TxID}\tENET\tfull\t{join(enet_select(Xtr, ytr, l1, alpha))}\n")

        rng = np.random.default_rng([SEED, int(TxID[4:])])
        for b, h in enumerate(halves(rng)):
            fh["halves"].write(f"{TxID}\t{b}\t{join(h)}\n")
            fh["selections"].write(f"{TxID}\tSKB\t{b}\t{join(skb_select(Xtr[h], ytr[h], k))}\n")
            if not a.skip_enet:
                fh["selections"].write(f"{TxID}\tENET\t{b}\t{join(enet_select(Xtr[h], ytr[h], l1, alpha))}\n")
        for f in fh.values():
            f.flush()
        print(f"shard {shard_i}: {gi + 1}/{len(genes)} {TxID} k={k}"
              + ("" if a.skip_enet else f" enet l1={l1} alpha={alpha:.3g} cv_r2={cv_r2:.3f}"), flush=True)
for f in fh.values():
    f.close()
# train / test sample h5 rows, shared by all genes (benchmark_susie.R reads X in this order)
if shard_i == 0:
    for name in ["train", "test"]:
        np.savetxt(os.path.join(a.out_dir, f"bench_{name}_rows.txt"), counts[name][1], fmt="%d")
