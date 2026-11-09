"""SHAP of the original SVR models (shap_svr_refit.py) vs stability selection, CpG level and region level.

Region level: the |SHAP| of a gene's selected CpGs summed within each co-methylated region
(stability_regions_corr.py map), compared between stable and non-stable regions (region_freq_<tag>.tsv).
Relative SHAP = value / the gene's mean over its CpGs (CpG level) or regions (region level), so genes are
comparable.

A  per gene: median relative |SHAP| of stable vs non-stable units, CpG level and region level
B  size-matched: regions with exactly 1, 2, 3 selected CpGs
C  per-gene Spearman correlation between selection frequency and |SHAP|, CpG vs region level
D  how often a gene's top-SHAP unit is stable, observed vs expected by chance
Writes <out_prefix>.png/.pdf and <out_prefix>_regions.tsv (one row per gene x region).
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

p = argparse.ArgumentParser()
p.add_argument("--shap_dir", required=True)
p.add_argument("--stab_dir", required=True)
p.add_argument("--tag", default="corr0.5")
p.add_argument("--out_prefix", required=True)
a = p.parse_args()

d = pd.concat([pd.read_csv(f, sep="\t", on_bad_lines="skip") for f in glob.glob(os.path.join(a.shap_dir, "shap_shard*.tsv"))])
d = d.dropna(subset=["same_selection"])
d = d[d.groupby("TxID")["CpGID"].transform("size") == d["k"]]                       # complete genes only
d["in_stable_region"] = d["in_stable_region"].astype(str) == "True"
rm = pd.concat([pd.read_csv(f, sep="\t", header=None, names=["TxID", "n_regions", "CpGID", "region", "start", "end", "n_cpg"],
                            usecols=["TxID", "CpGID", "region"])
                for f in glob.glob(os.path.join(a.stab_dir, f"region_map_{a.tag}_shard*.tsv"))])
rf = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t", usecols=["TxID", "region", "freq", "stable"])
d = d.merge(rm, on=["TxID", "CpGID"], how="left")
r = d.groupby(["TxID", "region"]).agg(shap=("mean_abs_shap", "sum"), n_sel=("CpGID", "size")).reset_index() \
    .merge(rf, on=["TxID", "region"], how="left")
r["stable"] = r["stable"].fillna(False).astype(bool)
r["freq"] = r["freq"].fillna(0)
r.to_csv(f"{a.out_prefix}_regions.tsv", sep="\t", index=False)

d["rel"] = d["mean_abs_shap"] / d.groupby("TxID")["mean_abs_shap"].transform("mean")
r["rel"] = r["shap"] / r.groupby("TxID")["shap"].transform("mean")


def paired(df, flag):
    g = df.groupby(["TxID", flag])["rel"].median().unstack().dropna()
    g.columns = ["non-stable", "stable"]
    return g


cpg, reg = paired(d, "in_stable_region"), paired(r, "stable")
n_genes = d["TxID"].nunique()
blue, grey, ink = "#2a78d6", "#b7b7b7", "#555555"
fig, axes = plt.subplots(1, 4, figsize=(20, 4.8), gridspec_kw={"width_ratios": [1.2, 1.1, 1, 0.9]})

# A
ax = axes[0]
data, pos, cols = [], [], []
for i, (g, lab) in enumerate([(cpg, "CpG level"), (reg, "region level\n(sum |SHAP| per region)")]):
    for j, (c, col) in enumerate([("non-stable", grey), ("stable", blue)]):
        data.append(g[c].to_numpy()); pos.append(i * 3 + j); cols.append(col)
    pval = wilcoxon(g["stable"], g["non-stable"]).pvalue
    ax.text(i * 3 + 0.5, 1.02, f"stable > non-stable in {100 * (g['stable'] > g['non-stable']).mean():.0f}% genes\n"
            f"p = {pval:.1e}", transform=ax.get_xaxis_transform(), ha="center", fontsize=8, color=ink)
bp = ax.boxplot(data, positions=pos, widths=0.7, patch_artist=True, showfliers=False,
                medianprops=dict(color="white", lw=2), whiskerprops=dict(color=ink), capprops=dict(color=ink))
for b, c in zip(bp["boxes"], cols):
    b.set(facecolor=c, edgecolor=c)
ax.set_xticks([0.5, 3.5], ["CpG level", "region level\n(sum |SHAP| per region)"], fontsize=9)
ax.set_ylabel("relative |SHAP| (1 = gene's mean unit)")
ax.set_title("A  Stable vs non-stable, per gene", fontsize=10, loc="left", pad=28)
ax.legend([bp["boxes"][0], bp["boxes"][1]], ["non-stable", "stable"], frameon=False, fontsize=8, loc="upper left")

# B
ax = axes[1]
sizes = [1, 2, 3]
for j, (c, col) in enumerate([("non-stable", grey), ("stable", blue)]):
    med, lo, hi = [], [], []
    for n in sizes:
        x = r[r["n_sel"] == n].assign(rel=lambda z: z["shap"] / z.groupby("TxID")["shap"].transform("mean"))
        g = paired(x, "stable")
        med.append(g[c].median()); lo.append(g[c].quantile(0.25)); hi.append(g[c].quantile(0.75))
    xs = np.arange(len(sizes)) + (j - 0.5) * 0.35
    ax.bar(xs, med, width=0.33, color=col, label=c)
    ax.errorbar(xs, med, yerr=[np.array(med) - lo, np.array(hi) - med], fmt="none", ecolor=ink, lw=1, capsize=3)
for i, n in enumerate(sizes):
    x = r[r["n_sel"] == n].assign(rel=lambda z: z["shap"] / z.groupby("TxID")["shap"].transform("mean"))
    g = paired(x, "stable")
    ax.text(i, 1.02, f"{len(g)} genes\n{100 * (g['stable'] > g['non-stable']).mean():.0f}% stable >",
            transform=ax.get_xaxis_transform(), ha="center", fontsize=8, color=ink)
ax.set_xticks(range(len(sizes)), [f"{n} selected CpG{'s' if n > 1 else ''}" for n in sizes], fontsize=9)
ax.set_ylabel("relative |SHAP| within size class (median, IQR)")
ax.set_title("B  Same-size regions", fontsize=10, loc="left", pad=28)
ax.legend(frameon=False, fontsize=8, loc="upper right")

# C
ax = axes[2]


def rho(df, f, s):
    return df.groupby("TxID").apply(lambda x: spearmanr(x[f], x[s]).statistic if x[f].nunique() > 1 else np.nan,
                                    include_groups=False).dropna()


rc, rr = rho(d, "freq", "mean_abs_shap"), rho(r, "freq", "shap")
bins = np.linspace(-1, 1, 41)
ax.hist(rc, bins=bins, color=grey, alpha=0.9, label=f"CpG level (median {rc.median():.2f})")
ax.hist(rr, bins=bins, color=blue, alpha=0.75, label=f"region level (median {rr.median():.2f})")
ax.axvline(0, color=ink, lw=0.8)
ax.set_xlabel("Spearman ρ, selection frequency vs |SHAP| (per gene)")
ax.set_ylabel("genes")
ax.set_title("C  Frequency vs SHAP within genes", fontsize=10, loc="left", pad=28)
ax.legend(frameon=False, fontsize=8, loc="upper left")

# D
ax = axes[3]
top_c = d.loc[d.groupby("TxID")["mean_abs_shap"].idxmax()]
top_r = r.loc[r.groupby("TxID")["shap"].idxmax()]
obs = [100 * top_c["in_stable_region"].mean(), 100 * top_r["stable"].mean()]
exp = [100 * d.groupby("TxID")["in_stable_region"].mean().mean(), 100 * r.groupby("TxID")["stable"].mean().mean()]
xs = np.arange(2)
ax.bar(xs - 0.18, exp, width=0.34, color=grey, label="expected by chance")
ax.bar(xs + 0.18, obs, width=0.34, color=blue, label="observed")
for x, v in zip(xs + 0.18, obs):
    ax.text(x, v + 1.5, f"{v:.0f}%", ha="center", fontsize=8, color=ink)
for x, v in zip(xs - 0.18, exp):
    ax.text(x, v + 1.5, f"{v:.0f}%", ha="center", fontsize=8, color=ink)
ax.set_xticks(xs, ["top-SHAP CpG", "top-SHAP region"], fontsize=9)
ax.set_ylim(0, 100)
ax.set_ylabel("% genes where it is stable")
ax.set_title("D  Is the top-SHAP unit stable?", fontsize=10, loc="left", pad=28)
ax.legend(frameon=False, fontsize=8, loc="upper left")

for ax in axes:
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e5e5e5", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"SHAP of the original SVR models vs stability selection ({n_genes} genes; "
             f"regions: co-methylated, {a.tag})", fontsize=11)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(f"{a.out_prefix}.{ext}", dpi=200)
print(f"{n_genes} genes; CpG: {100 * (cpg['stable'] > cpg['non-stable']).mean():.0f}% | region: "
      f"{100 * (reg['stable'] > reg['non-stable']).mean():.0f}% | rho cpg {rc.median():.2f} region {rr.median():.2f} | "
      f"top region stable {obs[1]:.0f}% vs {exp[1]:.0f}%")
