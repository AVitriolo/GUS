"""Stability selection on the models' own FEATURE matrices (cluster_cpgs: 1 -> features are CpG clusters).

Same procedure as analyses/scripts/stability_selection_kbest.py (Shah & Samworth 2013 complementary pairs):
for each gene, SelectKBest(f_regression) with the gene's tuned k (= number of saved selected features) is refit
on 2 x n_pairs half-samples of the TRAIN samples, each half stratified by condition. The difference is the
candidate set: the columns of the model's training matrix (data/xgb_clustered<TAG>/..._train.txt, the exact
input of model_selection; StandardScaler before SelectKBest does not change F-scores), so a candidate is a
cluster of co-methylated CpGs and no extra region step is needed.

Output, in the format of stability_selection_kbest.py (read unchanged by stability_selection_summary.R):
  subsample_selections_shard<i>.tsv  TxID, subsample, comma-separated 0-based candidate indices
  candidates_shard<i>.tsv            TxID, n_cand, k, index, featureID, saved_selected  (featureID = cluster ID)
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd
from sklearn.feature_selection import SelectKBest, f_regression

p = argparse.ArgumentParser()
p.add_argument("--sel_dir", required=True, help="results/sel_feats<TAG> (one file per model, lines <TxID>_<feature>)")
p.add_argument("--xgb_dir", required=True, help="data/xgb_clustered<TAG> (<prefix>_<TxID>_xgb_input_clustered_train.txt)")
p.add_argument("--prefix", required=True, help="e.g. hg38_v38_pc_1000000_25_10")
p.add_argument("--coldata", required=True, help="meth_colData.tsv (sample -> condition)")
p.add_argument("--n_pairs", type=int, default=100)
p.add_argument("--seed", type=int, default=1)
p.add_argument("--shard", default="0/1", help="i/n: process every n-th gene starting at i")
p.add_argument("--out_dir", required=True)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
shard_i, n_shards = map(int, a.shard.split("/"))

condition = pd.read_csv(a.coldata, sep="\t", index_col=0)["condition"]

sel_files = sorted(glob.glob(os.path.join(a.sel_dir, "*")))[shard_i::n_shards]
f_sel = open(os.path.join(a.out_dir, f"subsample_selections_shard{shard_i}.tsv"), "w")
f_cand = open(os.path.join(a.out_dir, f"candidates_shard{shard_i}.tsv"), "w")
n_done = 0
for i, f in enumerate(sel_files):
    TxID = os.path.basename(f).split("_")[-1]
    sel = {l.strip() for l in open(f) if l.strip()}
    train_file = os.path.join(a.xgb_dir, f"{a.prefix}_{TxID}_xgb_input_clustered_train.txt")
    if not sel or not os.path.exists(train_file):
        continue
    M = pd.read_csv(train_file, sep="\t", index_col=0)
    feats = [c for c in M.columns if c != TxID]
    X = M[feats].to_numpy(dtype=np.float64)
    y = M[TxID].to_numpy(dtype=np.float64)
    cond = condition.reindex(M.index).to_numpy()
    assert not pd.isna(cond).any(), f"{TxID}: train samples missing from {a.coldata}"
    k = len(sel)
    assert sel <= set(feats), f"{TxID}: saved selected features not in the training matrix"

    rng = np.random.default_rng([a.seed, int(TxID[4:])])
    ever = set()
    for b in range(a.n_pairs):
        half = np.zeros(len(y), dtype=bool)
        for lvl in np.unique(cond):
            idx = rng.permutation(np.where(cond == lvl)[0])
            half[idx[: len(idx) // 2]] = True
        for h, part in enumerate((half, ~half)):
            s = SelectKBest(f_regression, k=k).fit(X[part], y[part]).get_support(indices=True)
            ever.update(s)
            f_sel.write(f"{TxID}\t{2 * b + h}\t{','.join(map(str, s))}\n")
    for n, c in enumerate(feats):
        if n in ever or c in sel:
            f_cand.write(f"{TxID}\t{len(feats)}\t{k}\t{n}\t{c[len(TxID) + 1:]}\t{int(c in sel)}\n")
    n_done += 1
    if i % 50 == 0:
        print(f"shard {shard_i}: {i}/{len(sel_files)}", flush=True)
f_sel.close()
f_cand.close()
print(f"shard {shard_i}: {n_done} genes done")
