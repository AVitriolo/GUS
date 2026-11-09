"""Stable functional CpGs (stability selection, co-methylated regions) vs the Methodical TMRs (Heery & Schaefer 2025;
resources/TMRs/{N,T,D}_prostate_tmrs.bed, transcript ID first in the name column).
  A  distance of every stable gene-CpG link to the TSS of the linked gene model, with Methodical's search window
     (transcript +/- 5 kb, shown as +/- 5 kb of the TSS) marked
  B  genes with stable CpGs: stable CpG inside a TMR of the same gene (any transcript) / gene has TMRs elsewhere /
     gene has no TMR
Usage: python plot_stable_vs_tmr.py --stab_dir DIR --out_prefix PREFIX [--tag corr0.5]
"""
import argparse
import glob
import gzip
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--stab_dir", required=True)
p.add_argument("--out_prefix", required=True)
p.add_argument("--tag", default="corr0.5")
p.add_argument("--tmr_dir", default="/home/Prostate_GAS/resources/TMRs")
p.add_argument("--tss", default="/home/Prostate_GAS/resources/reference_data/gencode/hg38_v38_pc.bed")
p.add_argument("--gtf", default="/home/Prostate_GAS/resources/reference_data/gencode/gencode.v38.annotation.gtf.gz")
a = p.parse_args()

tx2g = {}
for l in gzip.open(a.gtf, "rt"):
    if "\ttranscript\t" in l:
        tx2g[re.search(r'transcript_id "(ENST\d+)', l).group(1)] = re.search(r'gene_name "([^"]+)', l).group(1)

# stable gene-CpG links: selected CpGs in a stable co-methylated region
cpg = pd.read_csv(os.path.join(a.stab_dir, "cpg_freq.tsv"), sep="\t", usecols=["TxID", "CpGID", "saved_selected"])
cpg = cpg[cpg.saved_selected == 1]
rm = pd.concat([pd.read_csv(f, sep="\t", header=None, usecols=[0, 2, 3], names=["TxID", "CpGID", "region"])
                for f in glob.glob(os.path.join(a.stab_dir, f"region_map_{a.tag}_shard*.tsv"))]).drop_duplicates(["TxID", "CpGID"])
rf = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t", usecols=["TxID", "region", "stable"])
st = cpg.merge(rm, on=["TxID", "CpGID"]).merge(rf[rf.stable], on=["TxID", "region"])[["TxID", "CpGID"]]
st["chr"] = st.CpGID.str.rsplit("_", n=1).str[0]
st["pos"] = st.CpGID.str.rsplit("_", n=1).str[1].astype(int)
st["gene"] = st.TxID.map(tx2g)
tss = pd.read_csv(a.tss, sep="\t", header=None, names=["chr", "s", "e", "TxID"]).drop_duplicates("TxID")
st["dist"] = (st.pos - st.TxID.map(dict(zip(tss.TxID, tss.s)))).abs()

tmr = pd.concat([pd.read_csv(f, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "s", "e", "name"])
                 for f in glob.glob(os.path.join(a.tmr_dir, "*_prostate_tmrs.bed"))]).drop_duplicates(["chr", "s", "e", "name"])
tmr["gene"] = tmr.name.str.split("_").str[0].map(tx2g)
x = st.merge(tmr[["gene", "chr", "s", "e"]], on=["gene", "chr"])
in_own = set(x[(x.pos >= x.s) & (x.pos <= x.e)].gene)
genes = set(st.gene.dropna()); tmr_genes = set(tmr.gene.dropna())
cats = [("stable CpG inside a TMR\nof the same gene", len(genes & in_own), "#1baf7a"),
        ("gene has TMRs, stable CpGs\noutside them", len((genes & tmr_genes) - in_own), "#8fd8bb"),
        ("gene has no TMR", len(genes - tmr_genes), "#c9c8c3")]

d = st.dist.dropna().clip(lower=1)
plt.rcParams.update({"font.size": 9})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(9.6, 3.6), layout="constrained", gridspec_kw={"width_ratios": [1.35, 1]})
bins = np.logspace(0, np.log10(1.05e6), 50)
ax.hist(d, bins=bins, color="#1baf7a", edgecolor="white", linewidth=0.4)
ax.set_xscale("log")
ax.axvspan(1, 5e3, color="#1f5f6b", alpha=0.10, lw=0)
ax.axvline(5e3, color="#1f5f6b", lw=1, ls="--")
ax.text(4.2e3, ax.get_ylim()[1] * 0.95, f"Methodical window\n(±5 kb): {100 * (d <= 5e3).mean():.0f}%",
        ha="right", va="top", fontsize=8, color="#1f5f6b")
ax.text(6e3, ax.get_ylim()[1] * 0.95, f"beyond: {100 * (d > 5e3).mean():.0f}%\n(> 50 kb: {100 * (d > 5e4).mean():.0f}%)",
        ha="left", va="top", fontsize=8, color="#0b0b0b")
ax.set_xlabel("distance of the stable CpG to the TSS of its gene (bp)")
ax.set_ylabel("stable gene–CpG links")
ax.spines[["top", "right"]].set_visible(False)
y = np.arange(len(cats))[::-1]
bx.barh(y, [c[1] for c in cats], color=[c[2] for c in cats], height=0.62)
for yi, (_, n, _) in zip(y, cats):
    bx.text(n + 80, yi, f"{n:,} ({100 * n / len(genes):.0f}%)", va="center", fontsize=8)
bx.set_yticks(y, [c[0] for c in cats], fontsize=8)
bx.tick_params(axis="y", length=0)
bx.set_xlim(0, max(c[1] for c in cats) * 1.35)
bx.set_xlabel(f"genes with stable CpGs (n = {len(genes):,})")
bx.spines[["top", "right"]].set_visible(False)
for e in ("png", "pdf"):
    fig.savefig(f"{a.out_prefix}.{e}", dpi=300)
print(f"links {len(d):,}; <=5kb {100 * (d <= 5e3).mean():.1f}%; >50kb {100 * (d > 5e4).mean():.1f}%;",
      [(c[0].replace(chr(10), ' '), c[1]) for c in cats], "| TMRs", len(tmr), "genes with TMR", len(tmr_genes))
