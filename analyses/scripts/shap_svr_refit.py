"""SHAP values of the original SVR models vs stability (genestd 2Mb v2 beta, k-capped run).

The fitted models were deleted after each run, so each gene's pipeline is refit with its saved
hyperparameters, exactly as model_selection_split_kbest_repeated_kcapped.py builds it:
StandardScaler (all candidate CpGs) -> SelectKBest(f_regression, k) -> SVR(C, epsilon, gamma, kernel),
on the train samples. The refit selection is checked against the saved sel_feats (same_selection).
SHAP: shap.KernelExplainer on the SVR (the explainer the pipeline uses for non-tree models), with a
k-means background of 20 and nsamples_base + 2k coalition samples (the pipeline's nsamples=100 is too
few for k up to 129), on n_explain random train samples; mean |SHAP| per selected CpG.
l1_reg=False by default: shap's default "auto" runs an AIC lasso first whenever few coalitions are
sampled and sets the SHAP of every CpG it drops to exactly 0 (29% of CpGs in k = 129 models), i.e.
among correlated neighbours one keeps the credit and the others get none.

Joined with stability selection (stability_selection_kbest.py / _regions.R): selection frequency and
whether the CpG lies in a stable co-methylated region. One row per gene x selected CpG.
"""
import argparse
import os

import h5py
import numpy as np
import pandas as pd
import shap
from sklearn.feature_selection import SelectKBest, f_regression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

p = argparse.ArgumentParser()
p.add_argument("--genes", required=True)
p.add_argument("--sel_dir", required=True)
p.add_argument("--hyper_dir", required=True)
p.add_argument("--stab_dir", required=True, help="stability selection dir (cpg_freq.tsv, region_freq_<tag>.tsv)")
p.add_argument("--tag", default="corr0.5")
p.add_argument("--cand_dir", required=True)
p.add_argument("--h5", required=True)
p.add_argument("--positions", required=True)
p.add_argument("--coldata", required=True)
p.add_argument("--counts_train", required=True)
p.add_argument("--prefix", default="hg38_v38_pc_50000_25_10")
p.add_argument("--n_explain", type=int, default=50)
p.add_argument("--nsamples_base", type=int, default=2048)
p.add_argument("--l1_reg", default="False", help='"False" (no regularisation) or a shap l1_reg value, e.g. "auto"')
p.add_argument("--shard", default="0/1")
p.add_argument("--out_dir", required=True)
a = p.parse_args()
shard_i, n_shards = map(int, a.shard.split("/"))
os.makedirs(a.out_dir, exist_ok=True)

genes = [g.strip() for g in open(a.genes) if g.strip()][shard_i::n_shards]
gset = set(genes)
freq = pd.read_csv(os.path.join(a.stab_dir, "cpg_freq.tsv"), sep="\t", usecols=["TxID", "CpGID", "freq", "stable"])
freq = freq[freq["TxID"].isin(gset)].set_index(["TxID", "CpGID"])
reg = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t")
reg = reg[reg["stable"] & reg["TxID"].isin(gset)].groupby("TxID")

pos = pd.read_csv(a.positions, sep="\t")
chr_cols = {c: (g.index.to_numpy(), g["pos"].to_numpy()) for c, g in pos.groupby("chr")}
del pos
row_of = {s: i for i, s in enumerate(pd.read_csv(a.coldata, sep="\t", index_col=0).index)}
counts = pd.read_csv(a.counts_train, sep="\t", index_col=0)
counts = counts[[s for s in counts.columns if s in row_of]]
rows = np.array([row_of[s] for s in counts.columns])

# resumable: genes already complete in this shard's file (all k rows written) are skipped, the rest appended
out_path = os.path.join(a.out_dir, f"shap_shard{shard_i}.tsv")
done = set()
if os.path.exists(out_path):
    prev = pd.read_csv(out_path, sep="\t", on_bad_lines="skip")
    ok = prev.groupby("TxID").agg(n=("CpGID", "size"), k=("k", "first"))
    done = set(ok.index[ok["n"] == ok["k"]])
    prev[prev["TxID"].isin(done)].to_csv(out_path, sep="\t", index=False)      # drop a gene cut mid-write
    out_f = open(out_path, "a")
else:
    out_f = open(out_path, "w")
    out_f.write("TxID\tCpGID\tk\tmean_abs_shap\tshap_rank\tfreq\tstable_cpg\tin_stable_region\tsame_selection\n")
genes = [g for g in genes if g not in done]
print(f"shard {shard_i}: {len(done)} genes already done, {len(genes)} to do", flush=True)
with h5py.File(a.h5, "r") as h5:
    d = h5["M"]
    for gi, TxID in enumerate(genes):
        hp = pd.read_csv(os.path.join(a.hyper_dir, f"{a.prefix}_3_{TxID}"), sep="\t", header=None, index_col=0)[1]
        k = int(hp["feature_selection__k"])
        saved = {l.strip().replace(f"{TxID}_", "", 1) for l in open(os.path.join(a.sel_dir, f"{a.prefix}_3_{TxID}"))
                 if l.strip()}
        cand = pd.read_csv(os.path.join(a.cand_dir, f"{a.prefix}_{TxID}.bed"), sep="\t", header=None,
                           dtype=str)[1].drop_duplicates()
        cols_c, pos_c = chr_cols[cand.str.split("_").str[0].iloc[0]]
        cpos = cand.str.split("_").str[1].astype(int).to_numpy()
        j = np.searchsorted(pos_c, cpos)
        ok = (j < len(pos_c)) & (pos_c[np.minimum(j, len(pos_c) - 1)] == cpos)
        cand, cols = cand[ok].to_numpy(), cols_c[j[ok]]
        order = np.argsort(cols)
        cand, cols = cand[order], cols[order]
        X = d[:, cols.min():cols.max() + 1][:, cols - cols.min()][rows].astype(np.float64)
        y = counts.loc[TxID].to_numpy(dtype=np.float64)

        pipe = Pipeline([("scaling", StandardScaler()), ("feature_selection", SelectKBest(f_regression, k=k)),
                         ("model", SVR(C=float(hp["model__C"]), epsilon=float(hp["model__epsilon"]),
                                       gamma=float(hp["model__gamma"]), kernel=hp["model__kernel"]))]).fit(X, y)
        idx = pipe.named_steps["feature_selection"].get_support(indices=True)
        sel = cand[idx]
        same = len(set(sel) & saved) / max(len(saved), 1)

        Xs = pipe[:-1].transform(X)
        rng = np.random.default_rng(int(TxID[4:]))
        expl_rows = rng.choice(len(Xs), size=min(a.n_explain, len(Xs)), replace=False)
        explainer = shap.KernelExplainer(pipe.named_steps["model"].predict, shap.kmeans(Xs, 20))
        sv = np.array(explainer.shap_values(Xs[expl_rows], nsamples=2 * k + a.nsamples_base, silent=True,
                                            l1_reg=False if a.l1_reg == "False" else a.l1_reg))
        mabs = np.abs(sv).mean(0)
        rank = (-mabs).argsort().argsort() + 1

        st = reg.get_group(TxID) if TxID in reg.groups else pd.DataFrame(columns=["start", "end"])
        for c, m, r in zip(sel, mabs, rank):
            pc = int(c.split("_")[1])
            in_reg = bool(((st["start"] <= pc) & (st["end"] >= pc)).any())
            fr, sc = freq.loc[(TxID, c)] if (TxID, c) in freq.index else (0.0, False)
            out_f.write(f"{TxID}\t{c}\t{k}\t{m:.6g}\t{r}\t{fr}\t{sc}\t{in_reg}\t{same:.3f}\n")
        out_f.flush()
        print(f"shard {shard_i}: {gi + 1}/{len(genes)} {TxID} k={k} same_selection={same:.2f}", flush=True)
out_f.close()
