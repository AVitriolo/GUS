"""Feature-selector benchmark, step 3: per gene x method, on the FULL-train selections.

  n_cpg, n_regions   size of the selection (regions = co-methylated regions from benchmark_selectors.py)
  test_r2            ridge (RidgeCV, standardised) refit on train with ONLY the selected CpGs, R2 on test
                     (0 selected -> predict the train mean)
  jaccard_<M>        region-level Jaccard overlap with method M's selection
  n_consensus        regions selected by >= --min_methods of the methods present (same value for every method of a gene)

Reads bench_inputs/selections shards in <bench_dir>; writes evaluate_genes.tsv there.
"""
import argparse
import glob
import itertools
import os

import h5py
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.metrics import r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

p = argparse.ArgumentParser()
p.add_argument("--bench_dir", required=True)
p.add_argument("--h5", required=True)
p.add_argument("--min_methods", type=int, default=2)
a = p.parse_args()
METHODS = ["SKB", "ENET", "SUSIE", "BORUTA"]          # those present in the selections are used


def read(pattern, names):
    return pd.concat([pd.read_csv(f, sep="\t", header=None, names=names, dtype=str, keep_default_na=False)
                      for f in sorted(glob.glob(os.path.join(a.bench_dir, pattern)))])


inp = read("bench_inputs_shard*.tsv", ["TxID", "what", "v"]).set_index(["TxID", "what"])["v"]
sel = read("bench_selections*_shard*.tsv", ["TxID", "method", "b", "idx"])
sel = sel[sel["b"] == "full"].set_index(["TxID", "method"])["idx"]
rows = {n: np.loadtxt(os.path.join(a.bench_dir, f"bench_{n}_rows.txt"), dtype=int) for n in ["train", "test"]}


def vec(s, dtype=float):
    return np.array(s.split(","), dtype=dtype) if s else np.array([], dtype=dtype)


out = []
with h5py.File(a.h5, "r") as h5:
    d = h5["M"]
    for TxID in inp.index.get_level_values(0).unique():
        cols = vec(inp[(TxID, "cols")], int)
        region = vec(inp[(TxID, "region")], int)
        y = {n: vec(inp[(TxID, f"y_{n}")]) for n in ["train", "test"]}
        block = d[:, cols.min():cols.max() + 1][:, cols - cols.min()]
        X = {n: block[rows[n]].astype(np.float64) for n in ["train", "test"]}

        regs, res = {}, {}
        for m in METHODS:
            if (TxID, m) not in sel.index:
                continue
            idx = vec(sel[(TxID, m)], int)
            regs[m] = set(region[idx])
            if len(idx):
                model = make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(-3, 4, 15)))
                pred = model.fit(X["train"][:, idx], y["train"]).predict(X["test"][:, idx])
            else:
                pred = np.full(len(y["test"]), y["train"].mean())
            res[m] = {"TxID": TxID, "method": m, "n_cpg": len(idx), "n_regions": len(regs[m]),
                      "test_r2": r2_score(y["test"], pred)}
        counts = pd.Series([r for s in regs.values() for r in s]).value_counts()
        n_cons = int((counts >= a.min_methods).sum()) if len(counts) else 0
        for m1, m2 in itertools.permutations(regs, 2):
            u = regs[m1] | regs[m2]
            res[m1][f"jaccard_{m2}"] = len(regs[m1] & regs[m2]) / len(u) if u else np.nan
        for m in res:
            res[m]["n_consensus"] = n_cons
            res[m]["n_in_consensus"] = sum(counts.get(r, 0) >= a.min_methods for r in regs[m])
            out.append(res[m])
pd.DataFrame(out).to_csv(os.path.join(a.bench_dir, "evaluate_genes.tsv"), sep="\t", index=False)
print(pd.DataFrame(out).groupby("method")[["n_cpg", "n_regions", "test_r2"]].median())
