"""Co-methylated regions for region-level stability selection (stability_selection_regions.R).

For each gene in candidates_shard<i>.tsv (output of stability_selection_kbest.py), the gene's
candidate CpGs (2Mb bed) are sorted by position and split into regions: a new region starts when
the next CpG is more than max_gap bp away OR its Pearson correlation with the previous CpG
(beta, TRAIN samples, all conditions) is below the threshold. Expression is not used.
CpGs with undefined correlation (constant in train) start a new region.

Writes, per threshold, region_map_corr<thr>_shard<i>.tsv with one row per CpG that was selected in
any half-sample or by the model: TxID, n_regions, CpGID, region, region start, region end, region n CpGs.
With --bg_per_gene N, also bg_regions_corr<thr>_shard<i>.tsv: N random candidate regions per gene
(TxID, chr, start, end, n CpGs), the pool for the motif-enrichment background (stable_regions_motif_sets.py).
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--in_dir", required=True, help="stability selection dir with candidates_shard<i>.tsv")
p.add_argument("--cand_dir", required=True)
p.add_argument("--h5", required=True)
p.add_argument("--positions", required=True, help="h5_cpg_positions.tsv (chr, pos in h5 column order)")
p.add_argument("--coldata", required=True, help="meth_colData.tsv (row order = h5 rows)")
p.add_argument("--train_samples", required=True)
p.add_argument("--thresholds", default="0.3,0.5,0.7")
p.add_argument("--max_gap", type=int, default=300)
p.add_argument("--prefix", default="hg38_v38_pc_50000_25_10")
p.add_argument("--shard", type=int, required=True)
p.add_argument("--bg_per_gene", type=int, default=0,
               help="also write bg_regions_corr<thr>_shard<i>.tsv: this many random candidate regions per gene "
                    "(motif-enrichment background), sampled without regard to selection")
p.add_argument("--skip_maps", action="store_true", help="only write the background regions")
a = p.parse_args()
thresholds = [float(t) for t in a.thresholds.split(",")]

pos = pd.read_csv(a.positions, sep="\t")
chr_cols = {c: (g.index.to_numpy(), g["pos"].to_numpy()) for c, g in pos.groupby("chr")}
del pos
samples = pd.read_csv(a.coldata, sep="\t", index_col=0).index
train = set(open(a.train_samples).read().split())
rows = np.array([i for i, s in enumerate(samples) if s in train])

needed = pd.read_csv(os.path.join(a.in_dir, f"candidates_shard{a.shard}.tsv"), sep="\t", header=None,
                     usecols=[0, 4], names=["TxID", "CpGID"])
needed = needed.groupby("TxID", sort=False)["CpGID"].apply(set)

outs = {} if a.skip_maps else \
    {t: open(os.path.join(a.in_dir, f"region_map_corr{t}_shard{a.shard}.tsv"), "w") for t in thresholds}
bgs = {t: open(os.path.join(a.in_dir, f"bg_regions_corr{t}_shard{a.shard}.tsv"), "w") for t in thresholds} \
    if a.bg_per_gene else {}
with h5py.File(a.h5, "r") as h5:
    d = h5["M"]
    for TxID, need in needed.items():
        cand = pd.read_csv(os.path.join(a.cand_dir, f"{a.prefix}_{TxID}.bed"), sep="\t", header=None,
                           dtype=str)[1].drop_duplicates()
        chrom = cand.str.split("_").str[0].iloc[0]
        cols_c, pos_c = chr_cols[chrom]
        cpos = cand.str.split("_").str[1].astype(int).to_numpy()
        j = np.searchsorted(pos_c, cpos)
        ok = (j < len(pos_c)) & (pos_c[np.minimum(j, len(pos_c) - 1)] == cpos)
        cand, cpos, cols = cand[ok].to_numpy(), cpos[ok], cols_c[j[ok]]
        order = np.argsort(cpos)
        cand, cpos, cols = cand[order], cpos[order], cols[order]

        X = d[:, cols.min():cols.max() + 1][:, cols - cols.min()][rows].astype(np.float64)
        Xc = X - X.mean(0)
        sd = np.sqrt((Xc ** 2).sum(0))
        with np.errstate(invalid="ignore", divide="ignore"):
            r_adj = (Xc[:, :-1] * Xc[:, 1:]).sum(0) / (sd[:-1] * sd[1:])      # corr(CpG i, CpG i+1)
        far = np.diff(cpos) > a.max_gap

        keep = np.isin(cand, list(need))
        for t in thresholds:
            new = np.concatenate([[True], far | ~(r_adj >= t)])               # NaN corr -> break
            region = np.cumsum(new)
            starts = np.flatnonzero(new)
            ends = np.append(starts[1:], len(cand)) - 1
            reg_start, reg_end, reg_n = cpos[starts], cpos[ends], ends - starts + 1
            if t in outs:
                for cg, rg in zip(cand[keep], region[keep]):
                    outs[t].write(f"{TxID}\t{region[-1]}\t{cg}\t{rg}\t{reg_start[rg - 1]}\t{reg_end[rg - 1]}\t{reg_n[rg - 1]}\n")
            if t in bgs:
                rng = np.random.default_rng([int(TxID[4:]), int(t * 100)])
                for rg in rng.choice(region[-1], size=min(a.bg_per_gene, region[-1]), replace=False):
                    bgs[t].write(f"{TxID}\t{chrom}\t{reg_start[rg]}\t{reg_end[rg]}\t{reg_n[rg]}\n")
for f in list(outs.values()) + list(bgs.values()):
    f.close()
print(f"shard {a.shard}: {len(needed)} genes done")
