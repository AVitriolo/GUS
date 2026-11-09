"""Coherence of the WGBS gene models with TCGA PRAD, split by stability.

WGBS: the footprint files of stable_footprint_files.py (gene, cv_r2, probe; selected CpGs mapped to the HM450
probes whose 50-nt footprint contains them): selected CpGs in a stable co-methylated region, the others, all.
TCGA PRAD: the local same-gene models (sel_feats/<gene>_SVR.tsv) and their candidate probes (the columns of
<data_dir>/<gene>_train.tsv). As plot_wgbs_vs_tcga_feature_coherence.py, per gene:
    observed = % of the WGBS probes also selected by the TCGA model of the same gene
    chance   = (% of the WGBS probes that are TCGA candidates) * (TCGA selected / TCGA candidates)
Genes with >= 1 WGBS probe and a TCGA model. Wilcoxon signed-rank observed vs chance per set and CV R2 bin.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

p = argparse.ArgumentParser()
p.add_argument("--footprint_dir", required=True)
p.add_argument("--tcga_sel_dir", required=True)
p.add_argument("--tcga_data_dir", required=True)
p.add_argument("--out_prefix", required=True)
p.add_argument("--no_title", action="store_true", help="omit the figure title (slide version)")
a = p.parse_args()

SETS = [("selected", "selected (all)", "#5a5a5a"), ("unstable", "selected, not in a stable region", "#86b6ef"),
        ("stable", "selected, in a stable region", "#2a78d6")]
wg = {s: pd.read_csv(os.path.join(a.footprint_dir, f"wgbs_genestd_{s}_footprint.tsv.gz"), sep="\t") for s, _, _ in SETS}
wg = {s: d[d["probe"].str.startswith("cg")] for s, d in wg.items()}
genes = sorted(set().union(*[set(d["gene"]) for d in wg.values()]))

tcga = {}
for g in genes:
    sf = os.path.join(a.tcga_sel_dir, f"{g}_SVR.tsv")
    df = os.path.join(a.tcga_data_dir, f"{g}_train.tsv")
    if not (os.path.exists(sf) and os.path.exists(df)):
        continue
    with open(df) as fh:
        cand = {c for c in fh.readline().rstrip("\n").split("\t") if c.startswith("cg")}
    sel = {l.strip() for l in open(sf) if l.strip()}
    tcga[g] = (cand, sel)
print(f"{len(genes):,} WGBS genes with >= 1 probe; {len(tcga):,} with a TCGA PRAD model")

rows = []
for s, _, _ in SETS:
    for g, x in wg[s].groupby("gene"):
        if g not in tcga:
            continue
        cand, sel = tcga[g]
        pr = set(x["probe"])
        obs = 100 * len(pr & sel) / len(pr)
        exp = 100 * len(pr & cand) / len(pr) * len(sel) / max(1, len(cand))
        rows.append({"set": s, "gene": g, "cv_r2": x["cv_r2"].iloc[0], "n_probes": len(pr), "observed": obs, "chance": exp})
d = pd.DataFrame(rows)
edges = [-np.inf, 0, 0.1, 0.2, 0.3, 0.5, np.inf]
labels = ["<0", "0–0.1", "0.1–0.2", "0.2–0.3", "0.3–0.5", "≥0.5"]
d["bin"] = pd.cut(d["cv_r2"], edges, labels=labels, right=False)
d.to_csv(f"{a.out_prefix}_genes.tsv", sep="\t", index=False)


def summ(x):
    p = wilcoxon(x["observed"], x["chance"]).pvalue if len(x) > 10 and (x["observed"] != x["chance"]).any() else np.nan
    return pd.Series({"genes": len(x), "observed_mean": x["observed"].mean(), "chance_mean": x["chance"].mean(),
                      "ratio": x["observed"].mean() / x["chance"].mean(), "wilcoxon_p": p})


tab = d.groupby(["set", "bin"], observed=False).apply(summ, include_groups=False)
allg = d.groupby("set").apply(summ, include_groups=False)
tab.to_csv(f"{a.out_prefix}_summary.tsv", sep="\t")
print(allg.round(4).to_string())
print(tab.round(3).to_string())

ink = "#555555"
fig, axes = plt.subplots(1, 2, figsize=(14, 4.8))
x = np.arange(len(labels))
w = 0.26
ax = axes[0]
for i, (s, lab, c) in enumerate(SETS):
    t = tab.loc[s]
    ax.bar(x + (i - 1) * w, t["observed_mean"], width=w - 0.02, color=c, label=lab)
    ax.scatter(x + (i - 1) * w, t["chance_mean"], marker="_", s=180, color="#222222", zorder=3,
               label="chance" if i == 0 else None)
ax.set_xticks(x, [f"{l}\nn={int(tab.loc['selected'].loc[l, 'genes']):,}" for l in labels], fontsize=8)
ax.set_xlabel("WGBS model CV R² (n = genes with ≥ 1 probe, all selected)")
ax.set_ylabel("% of the WGBS probes also selected\nby the TCGA PRAD model (mean)")
ax.set_title("A  Coherence with TCGA PRAD (bars) vs chance (—)", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax = axes[1]
for s, lab, c in SETS:
    ax.plot(x, tab.loc[s]["ratio"], color=c, lw=2, marker="o", ms=6, label=lab)
ax.axhline(1, color=ink, lw=0.8, ls="--")
ax.set_xticks(x, labels, fontsize=8)
ax.set_xlabel("WGBS model CV R²")
ax.set_ylabel("observed / chance")
ax.set_title("B  Enrichment over chance", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left")
for ax in axes:
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e5e5e5", lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
if not a.no_title:
    fig.suptitle("Feature coherence of the WGBS gene models with TCGA PRAD (HM450, footprint matching), by stability",
                 fontsize=11)
fig.tight_layout()
for ext in ["png", "pdf"]:
    fig.savefig(f"{a.out_prefix}.{ext}", dpi=200)
