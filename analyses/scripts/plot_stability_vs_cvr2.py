"""Per model: % of its selected CpGs that are stable vs the model's CV R2 (genestd 2Mb v2 beta, k-capped).
A  stable as single CpGs (cpg_freq.tsv)
B  in a stable co-methylated region (region_freq_<tag>.tsv)
Points coloured by the model's k (no. selected CpGs); line = median per 0.05 CV R2 bin.
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

p = argparse.ArgumentParser()
p.add_argument("--stab_dir", required=True)
p.add_argument("--perf_dir", required=True)
p.add_argument("--tag", default="corr0.5")
p.add_argument("--out_prefix", required=True)
a = p.parse_args()

cv = {os.path.basename(f).split("_")[-1]: pd.read_csv(f, sep="\t", header=None, index_col=0)[1]["cv_r2"]
      for f in glob.glob(os.path.join(a.perf_dir, "*"))}
c = pd.read_csv(os.path.join(a.stab_dir, "cpg_freq.tsv"), sep="\t", usecols=["TxID", "CpGID", "k", "saved_selected", "stable"])
c = c[c["saved_selected"] == 1].copy()
c["pos"] = c["CpGID"].str.split("_").str[1].astype(int)
r = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t")
r = r[r["stable"]]
m = c.merge(r[["TxID", "start", "end"]], on="TxID", how="left")
m["in"] = (m["pos"] >= m["start"]) & (m["pos"] <= m["end"])
c = c.merge(m.groupby(["TxID", "CpGID"], as_index=False)["in"].any(), on=["TxID", "CpGID"])
g = c.groupby("TxID").agg(k=("k", "first"), pct_cpg=("stable", "mean"), pct_region=("in", "mean"))
g[["pct_cpg", "pct_region"]] *= 100
g["cv_r2"] = g.index.map(cv)
g = g.dropna(subset=["cv_r2"])
g.to_csv(f"{a.out_prefix}.tsv", sep="\t")

groups = [("k ≤ 13", g["k"] <= 13, "#2a78d6"), ("k 18–108", (g["k"] > 13) & (g["k"] < 129), "#1baf7a"),
          ("k = 129", g["k"] >= 129, "#eb6834")]
ink = "#555555"
fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), sharex=True)
edges = np.arange(np.floor(g["cv_r2"].min() / 0.05) * 0.05, g["cv_r2"].max() + 0.05, 0.05)
for ax, col, title in [(axes[0], "pct_cpg", "A  Stable as single CpGs"),
                       (axes[1], "pct_region", "B  In a stable co-methylated region (r ≥ 0.5, ≤ 300 bp)")]:
    for lab, mask, colr in groups:
        ax.scatter(g.loc[mask, "cv_r2"], g.loc[mask, col], s=5, alpha=0.35, color=colr, lw=0,
                   label=f"{lab} ({mask.sum():,} models)")
    b = pd.cut(g["cv_r2"], edges)
    med = g.groupby(b, observed=True)[col].median()
    ax.plot([iv.mid for iv in med.index], med.to_numpy(), color="#222222", lw=2, label="median per 0.05 bin")
    rho = spearmanr(g["cv_r2"], g[col]).statistic
    ax.text(0.02, 0.97, f"Spearman ρ = {rho:.2f}\n{len(g):,} models", transform=ax.transAxes, va="top",
            fontsize=9, color=ink)
    ax.set_xlabel("model CV R² (SVR)")
    ax.set_ylabel("% of the model's selected CpGs that are stable")
    ax.set_ylim(-2, 102)
    ax.set_title(title, fontsize=10, loc="left")
    ax.set_axisbelow(True)
    ax.grid(color="#e5e5e5", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
axes[1].legend(frameon=False, fontsize=8, loc="lower right", markerscale=3)
fig.suptitle("Stability of the selected CpGs vs model performance (PFER ≤ 1)", fontsize=11)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(f"{a.out_prefix}.{ext}", dpi=200)
print(g[["pct_cpg", "pct_region"]].describe().round(1))
