"""Mean beta of Normal / Tumour / metastasis per significant MethylDriver N_vs_T region, raw vs ComBat (v2 integration),
for the slide showing that ComBat (batch = dataset, all metastases in one dataset) pins metastases between N and T.
Usage: python combat_group_means_stable.py <tag> <out.tsv>
"""
import sys
import h5py
import numpy as np
import pandas as pd

R = "/home/Prostate_GAS/"
tag, out = sys.argv[1], sys.argv[2]
res = pd.read_csv(f"{R}results/methyldriver_2mb_v2_beta_genestd_stable_corr0.5/full_N_vs_T_{tag}.tsv", sep="\t",
                  usecols=["name", "gene", "q_value"])
sig = res[res.q_value < 0.05]
md = pd.read_csv(f"{R}data/to_methyldriver_w20/md_regions_{tag}.tsv", sep="\t", usecols=["mdID", "CpGs"])
md = md[md.mdID.isin(sig.name)]
long = md.assign(CpGID=md.CpGs.str.split(",")).explode("CpGID")[["mdID", "CpGID"]]
pos = pd.read_csv(f"{R}results/test_replication_2mb_v2_beta_genestd_kcapped/h5_cpg_positions.tsv", sep="\t",
                  dtype={"chr": str, "pos": np.int64})
idx = pd.Series(np.arange(len(pos)), index=(pos.chr + "_" + pos.pos.astype(str)).values)
long["i"] = idx.reindex(long.CpGID).values
long = long.dropna(subset=["i"]).astype({"i": np.int64})
cd = pd.read_csv(f"{R}resources/integrated_wgbs_merged_v2_hdf5/meth_colData.tsv", sep="\t")
grp = {"N": (cd.condition == "Normal").to_numpy(), "T": (cd.condition == "Tumour").to_numpy(),
       "M": (cd.condition == "metastasis").to_numpy()}
ui = np.unique(long.i)
for k, (f, ds) in {"raw": ("integrated_wgbs.h5", "beta"), "combat": ("integrated_wgbs_combat_beta.h5", "M")}.items():
    d = h5py.File(f"{R}resources/integrated_wgbs_merged_v2_hdf5/{f}", "r")[ds]
    X = np.concatenate([d[:, ui[a:a + 2000]] for a in range(0, len(ui), 2000)], axis=1)   # samples x CpGs
    for g, m in grp.items():
        long[f"{g}_{k}"] = pd.Series(100 * np.nanmean(X[m], axis=0), index=ui).reindex(long.i).values
    print(k, "done", flush=True)
o = long.drop(columns=["CpGID", "i"]).groupby("mdID").mean().reset_index().merge(sig, left_on="mdID", right_on="name")
o.to_csv(out, sep="\t", index=False)
print(len(o), "regions")
