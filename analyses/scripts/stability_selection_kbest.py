"""Stability selection with the models' own selector, step 1: subsample selections.

For each gene, repeats the models' feature selection - sklearn SelectKBest(f_regression) with the
gene's tuned k (= number of saved selected CpGs) - on B complementary pairs of half-samples of the
TRAIN samples (Shah & Samworth 2013), each half stratified by condition. Candidate CpGs = the
gene's 2Mb candidate set, as in the models.

Writes one line per subsample: TxID, subsample, comma-separated 0-based candidate indices, plus a
candidates table (TxID, n_cand, k, index, CpGID, saved_selected) for CpGs selected at least once. stability_selection_summary.R then computes
the Nogueira stability index (stabm) and the stabs error bounds.
"""
import argparse
import glob
import os

import h5py
import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_regression

p = argparse.ArgumentParser()
p.add_argument("--sel_dir", required=True)
p.add_argument("--cand_dir", required=True)
p.add_argument("--h5", required=True)
p.add_argument("--positions", required=True, help="h5_cpg_positions.tsv (chr, pos in h5 column order)")
p.add_argument("--coldata", required=True, help="meth_colData.tsv (row order = h5 rows)")
p.add_argument("--counts_train", required=True)
p.add_argument("--prefix", default="hg38_v38_pc_50000_25_10")
p.add_argument("--n_pairs", type=int, default=100)
p.add_argument("--seed", type=int, default=1)
p.add_argument("--shard", default="0/1", help="i/n: process every n-th gene starting at i")
p.add_argument("--out_dir", required=True)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
shard_i, n_shards = map(int, a.shard.split("/"))

pos = pd.read_csv(a.positions, sep="\t")
chr_cols = {c: (g.index.to_numpy(), g["pos"].to_numpy()) for c, g in pos.groupby("chr")}
del pos

coldata = pd.read_csv(a.coldata, sep="\t", index_col=0)
row_of = {s: i for i, s in enumerate(coldata.index)}
counts = pd.read_csv(a.counts_train, sep="\t", index_col=0)
counts = counts[[s for s in counts.columns if s in row_of]]
rows = np.array([row_of[s] for s in counts.columns])
cond = coldata["condition"].to_numpy()[rows]


def complementary_pairs(rng):
    """n_pairs x 2 half-samples, each half stratified by condition."""
    for _ in range(a.n_pairs):
        half = np.zeros(len(rows), dtype=bool)
        for lvl in np.unique(cond):
            idx = rng.permutation(np.where(cond == lvl)[0])
            half[idx[: len(idx) // 2]] = True
        yield half
        yield ~half


sel_files = sorted(glob.glob(os.path.join(a.sel_dir, "*")))[shard_i::n_shards]
f_sel = open(os.path.join(a.out_dir, f"subsample_selections_shard{shard_i}.tsv"), "w")
f_cand = open(os.path.join(a.out_dir, f"candidates_shard{shard_i}.tsv"), "w")
with h5py.File(a.h5, "r") as h5:
    d = h5["M"]
    for i, f in enumerate(sel_files):
        TxID = os.path.basename(f).split("_")[-1]
        sel = {l.strip().replace(f"{TxID}_", "", 1) for l in open(f) if l.strip()}
        cand_file = os.path.join(a.cand_dir, f"{a.prefix}_{TxID}.bed")
        if not sel or not os.path.exists(cand_file) or TxID not in counts.index:
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

        X = d[:, cols[0]:cols[-1] + 1][:, cols - cols[0]][rows]     # train samples x candidates
        y = counts.loc[TxID].to_numpy(dtype=np.float64)
        k = len(sel)
        rng = np.random.default_rng([a.seed, int(TxID[4:])])
        ever = set()
        for b, half in enumerate(complementary_pairs(rng)):
            s = SelectKBest(f_regression, k=k).fit(X[half], y[half]).get_support(indices=True)
            ever.update(s)
            f_sel.write(f"{TxID}\t{b}\t{','.join(map(str, s))}\n")
        # only CpGs selected at least once or in the saved selection; n_cand carries p
        for n, c in enumerate(cand):
            if n in ever or c in sel:
                f_cand.write(f"{TxID}\t{len(cand)}\t{k}\t{n}\t{c}\t{int(c in sel)}\n")
        if i % 50 == 0:
            print(f"shard {shard_i}: {i}/{len(sel_files)}", flush=True)
f_sel.close()
f_cand.close()
