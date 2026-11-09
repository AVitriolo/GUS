"""Per-cohort report of the stability selection (stability_selection_summary.R + stability_selection_regions.R).

Per gene (stability_report_genes.tsv): cv_r2, test_r2, k, n candidates, Nogueira stability phi, selected features
stable as single features, selected features lying in a stable co-methylated region (n and %), stable regions.
Plot (stability_report.png/.pdf): % of each model's selected features in a stable region, and phi, by CV R2 bin.
Usage: python stability_report.py --stab_dir DIR --perf_dir DIR [--tag corr0.5] [--label COHORT]
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
p.add_argument("--stab_dir", required=True)
p.add_argument("--perf_dir", required=True, help="performance files <gene>_SVR.tsv (or <prefix>_<gene>) with cv_r2")
p.add_argument("--tag", default="corr0.5")
p.add_argument("--label", default="")
a = p.parse_args()
S = a.stab_dir


def gene_of(f):
    b = os.path.basename(f)
    return b[:-len("_SVR.tsv")] if b.endswith("_SVR.tsv") else b.split("_")[-1]


perf = {}
for f in glob.glob(os.path.join(a.perf_dir, "*")):
    s = dict(l.rstrip("\n").split("\t")[:2] for l in open(f) if "\t" in l)
    perf[gene_of(f)] = (float(s.get("cv_r2", "nan")), float(s.get("test_r2", "nan")))

cpg = pd.read_csv(os.path.join(S, "cpg_freq.tsv"), sep="\t", dtype={"TxID": str})
genes = pd.read_csv(os.path.join(S, "stability_genes.tsv"), sep="\t", dtype={"TxID": str})
reg = pd.read_csv(os.path.join(S, f"region_freq_{a.tag}.tsv"), sep="\t", dtype={"TxID": str})
rmap = pd.concat([pd.read_csv(f, sep="\t", header=None, dtype={0: str}, usecols=[0, 2, 3], names=["TxID", "CpGID", "region"])
                  for f in glob.glob(os.path.join(S, f"region_map_{a.tag}_shard*.tsv"))]).drop_duplicates(["TxID", "CpGID"])

sel = cpg[cpg["saved_selected"] == 1].merge(rmap, on=["TxID", "CpGID"], how="left") \
    .merge(reg[["TxID", "region", "stable"]].rename(columns={"stable": "region_stable"}), on=["TxID", "region"], how="left")
sel["region_stable"] = sel["region_stable"].fillna(False).astype(bool)
g = sel.groupby("TxID").agg(n_selected=("CpGID", "size"), n_stable_single=("stable", "sum"),
                            n_in_stable_region=("region_stable", "sum"))
g["pct_stable_single"] = 100 * g["n_stable_single"] / g["n_selected"]
g["pct_in_stable_region"] = 100 * g["n_in_stable_region"] / g["n_selected"]
g["n_stable_regions"] = reg[reg["stable"] & (reg["n_saved"] > 0)].groupby("TxID").size().reindex(g.index).fillna(0).astype(int)
g = g.join(genes.set_index("TxID")[["k", "n_cand", "phi"]])
g["cv_r2"] = [perf.get(t, (np.nan, np.nan))[0] for t in g.index]
g["test_r2"] = [perf.get(t, (np.nan, np.nan))[1] for t in g.index]
g.index.name = "gene"
g.to_csv(os.path.join(S, "stability_report_genes.tsv"), sep="\t")

edges = [-np.inf, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, np.inf]
labels = ["<0", "0–0.1", "0.1–0.2", "0.2–0.3", "0.3–0.4", "0.4–0.5", "0.5–0.6", "0.6–0.7", "≥0.7"]
g["bin"] = pd.cut(g["cv_r2"], edges, labels=labels, right=False)
tab = g.groupby("bin", observed=False).agg(models=("k", "size"), median_pct_in_stable_region=("pct_in_stable_region", "median"),
                                           median_pct_stable_single=("pct_stable_single", "median"),
                                           pct_models_with_stable_region=("n_stable_regions", lambda s: 100 * (s > 0).mean()),
                                           median_phi=("phi", "median"))
tab.to_csv(os.path.join(S, "stability_report_bins.tsv"), sep="\t")

blue, ink = "#2a78d6", "#555555"
fig, axes = plt.subplots(1, 2, figsize=(15, 4.8))
for ax, col, ylab, ylim in [(axes[0], "pct_in_stable_region", "% of the model's selected features\nin a stable co-methylated region", (-2, 108)),
                            (axes[1], "phi", "stability index Φ̂ (Nogueira)", (-0.05, 1.05))]:
    data = [g.loc[g["bin"] == l, col].dropna().to_numpy() for l in labels]
    keep = [i for i, d in enumerate(data) if len(d)]
    bp = ax.boxplot([data[i] for i in keep], positions=[i + 1 for i in keep], widths=0.6, patch_artist=True,
                    showfliers=False, medianprops=dict(color="white", lw=2), whiskerprops=dict(color=ink),
                    capprops=dict(color=ink))
    for b in bp["boxes"]:
        b.set(facecolor=blue, edgecolor=blue)
    ax.set_xticks(range(1, len(labels) + 1), [f"{l}\nn={len(d):,}" for l, d in zip(labels, data)], fontsize=8)
    ax.set_xlabel("model CV R² (n = models)")
    ax.set_ylabel(ylab)
    ax.set_ylim(*ylim)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e5e5e5", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
for i, l in enumerate(labels, start=1):
    m = tab.loc[l, "median_pct_in_stable_region"]
    if not np.isnan(m):
        axes[0].text(i, 103, f"{m:.0f}%", ha="center", fontsize=8, color=ink)
axes[0].set_title("A  Selected features in stable regions (PFER ≤ 1)", fontsize=10, loc="left")
axes[1].set_title("B  Stability of each model's selection", fontsize=10, loc="left")
fig.suptitle(f"{a.label + ': ' if a.label else ''}stability selection, {len(g):,} models "
             f"({a.tag.replace('corr', 'regions r ≥ ')}, ≤ 300 bp)", fontsize=11)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(os.path.join(S, f"stability_report.{ext}"), dpi=200)
print(f"{a.label}: {len(g):,} models | selected features {int(g['n_selected'].sum()):,} | stable single "
      f"{100 * g['n_stable_single'].sum() / g['n_selected'].sum():.1f}% | in a stable region "
      f"{100 * g['n_in_stable_region'].sum() / g['n_selected'].sum():.1f}% | models with >= 1 stable region "
      f"{100 * (g['n_stable_regions'] > 0).mean():.1f}%")
print(tab.round(1).to_string())
