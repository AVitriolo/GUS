"""Distribution of SVR |SHAP| in stable vs non-stable co-methylated regions, all genes pooled.

Input: <shap_dir>/shap_vs_stability_regions.tsv (plot_shap_vs_stability.py: one row per gene x region,
summed |SHAP| of its selected CpGs) and the shap_shard*.tsv files (per CpG).
|SHAP| is divided by the gene's mean over its regions (A, B) or CpGs (C), log10 scale.
A  regions: stable vs non-stable
B  regions with exactly one selected CpG (size-matched)
C  CpGs: in a stable region vs not
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--shap_dir", required=True)
a = p.parse_args()

r = pd.read_csv(os.path.join(a.shap_dir, "shap_vs_stability_regions.tsv"), sep="\t")
r["rel"] = r["shap"] / r.groupby("TxID")["shap"].transform("mean")
d = pd.concat([pd.read_csv(f, sep="\t", on_bad_lines="skip") for f in glob.glob(os.path.join(a.shap_dir, "shap_shard*.tsv"))])
d = d.dropna(subset=["same_selection"])
d = d[d.groupby("TxID")["CpGID"].transform("size") == d["k"]]
d["stable"] = d["in_stable_region"].astype(str) == "True"
d["rel"] = d["mean_abs_shap"] / d.groupby("TxID")["mean_abs_shap"].transform("mean")
r1 = r[r["n_sel"] == 1].copy()
r1["rel"] = r1["shap"] / r1.groupby("TxID")["shap"].transform("mean")

grey, blue, ink = "#9a9a9a", "#2a78d6", "#555555"
fig, axes = plt.subplots(1, 3, figsize=(17, 4.6), sharey=False)
bins = np.linspace(-4, 1.5, 56)
for ax, df, flag, title, unit in [
        (axes[0], r, "stable", "A  Regions (sum of |SHAP| of the region's selected CpGs)", "regions"),
        (axes[1], r1, "stable", "B  Regions with exactly 1 selected CpG (size-matched)", "regions"),
        (axes[2], d, "stable", "C  Single CpGs (in a stable region or not)", "CpGs")]:
    for val, lab, col in [(False, "non-stable", grey), (True, "stable", blue)]:
        x = np.log10(df.loc[df[flag] == val, "rel"].dropna().clip(lower=1e-4))
        ax.hist(x, bins=bins, density=True, histtype="stepfilled", alpha=0.35, color=col)
        ax.hist(x, bins=bins, density=True, histtype="step", lw=1.8, color=col,
                label=f"{lab} (n={len(x):,}; median {10 ** np.median(x):.2f})")
        ax.axvline(np.median(x), color=col, ls="--", lw=1.2)
    ratio = df.loc[df[flag], "rel"].median() / df.loc[~df[flag], "rel"].median()
    ax.text(0.98, 0.95, f"median ratio {ratio:.1f}×", transform=ax.transAxes, ha="right", va="top", fontsize=9,
            color=ink)
    ax.set_xticks([-4, -3, -2, -1, 0, 1], ["0.0001", "0.001", "0.01", "0.1", "1", "10"])
    ax.set_xlabel("relative |SHAP| (1 = gene's mean; log scale)")
    ax.set_ylabel(f"density of {unit}")
    ax.set_title(title, fontsize=10, loc="left")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e5e5e5", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"SVR |SHAP| in stable vs non-stable co-methylated regions ({r['TxID'].nunique():,} genes, pooled)",
             fontsize=11)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(a.shap_dir, f"shap_distribution_stable_vs_nonstable.{ext}"), dpi=200)
