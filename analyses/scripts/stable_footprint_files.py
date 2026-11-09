"""WGBS footprint files for the TCGA coherence analysis (analyses/tcga_meteor_pancancer/scripts/confirm_analyses.py --wgbs),
split by stability: for every genestd gene model, its selected CpGs in a stable co-methylated region, the others,
and all of them, mapped to the HM450 probes whose footprint contains them (probe_footprint.py, as
wgbs_prostate_functional_footprint.tsv.gz). Columns: gene (symbol), cv_r2, probe; CpGs without a probe are kept
as chr_pos, as in the existing file (confirm_analyses.py keeps only cg probes).

Writes <out_dir>/wgbs_genestd_{stable,unstable,selected}_footprint.tsv.gz
"""
import argparse
import glob
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from probe_footprint import cpgs_to_probes, load_footprints

p = argparse.ArgumentParser()
p.add_argument("--stab_dir", required=True)
p.add_argument("--tag", default="corr0.5")
p.add_argument("--perf_dir", required=True)
p.add_argument("--manifest", default="resources/illumina_probes/HM450.hg38.manifest.tsv.gz")
p.add_argument("--conversion_table", default="resources/conversion_table_gencodeV38.tsv")
p.add_argument("--out_dir", required=True)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

c = pd.read_csv(os.path.join(a.stab_dir, "cpg_freq.tsv"), sep="\t", usecols=["TxID", "CpGID", "saved_selected"])
c = c[c["saved_selected"] == 1].copy()
c["pos"] = c["CpGID"].str.split("_").str[1].astype(int)
r = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t")
r = r[r["stable"]]
m = c.merge(r[["TxID", "start", "end"]], on="TxID", how="left")
m["in"] = (m["pos"] >= m["start"]) & (m["pos"] <= m["end"])
c = c.merge(m.groupby(["TxID", "CpGID"], as_index=False)["in"].any(), on=["TxID", "CpGID"])

conv = pd.read_csv(a.conversion_table, sep=" ", header=None, skiprows=1, usecols=[4, 5], names=["TxID", "gene"])
c["gene"] = c["TxID"].map(dict(zip(conv["TxID"], conv["gene"])))
cv = {os.path.basename(f).split("_")[-1]: pd.read_csv(f, sep="\t", header=None, index_col=0)[1]["cv_r2"]
      for f in glob.glob(os.path.join(a.perf_dir, "*"))}
c["cv_r2"] = c["TxID"].map(cv)
c = c.dropna(subset=["gene", "cv_r2"])
dup = c.groupby("gene")["TxID"].nunique()
if (dup > 1).any():
    print(f"WARNING: {int((dup > 1).sum())} genes with > 1 model; keeping the model with the highest cv_r2")
    best = c.drop_duplicates("TxID").sort_values("cv_r2", ascending=False).drop_duplicates("gene")["TxID"]
    c = c[c["TxID"].isin(set(best))]

hit = cpgs_to_probes(c["CpGID"].unique(), load_footprints(a.manifest))
for name, mask in [("stable", c["in"]), ("unstable", ~c["in"]), ("selected", c["in"] | ~c["in"])]:
    rows = [(g, v, pr) for g, v, k in zip(c.loc[mask, "gene"], c.loc[mask, "cv_r2"], c.loc[mask, "CpGID"])
            for pr in hit.get(k, [k])]
    out = pd.DataFrame(rows, columns=["gene", "cv_r2", "probe"]).drop_duplicates()
    out.to_csv(os.path.join(a.out_dir, f"wgbs_genestd_{name}_footprint.tsv.gz"), sep="\t", index=False,
               compression="gzip")
    cg = out[out["probe"].str.startswith("cg")]
    print(f"{name}: {len(out):,} rows, {cg['probe'].nunique():,} distinct probes, "
          f"{cg['gene'].nunique():,} genes with >= 1 probe")
