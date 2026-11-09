"""Stable-only selected features, in the sel_feats format (one file per model, lines <TxID>_<chr>_<pos>): for every
model, its selected CpGs lying in a stable co-methylated region (stability_selection_regions.R, region_freq_<tag>.tsv).
Models without any are not written. Input of the MethylDriver / gimme w20 pipeline (sel_feats_to_regions_w20.R).
Usage: python stable_sel_feats.py --stab_dir DIR --out_dir DIR [--tag corr0.5] [--prefix hg38_v38_pc_50000_25_10_3]
"""
import argparse
import os

import pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--stab_dir", required=True)
p.add_argument("--out_dir", required=True)
p.add_argument("--tag", default="corr0.5")
p.add_argument("--prefix", default="hg38_v38_pc_50000_25_10_3")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

c = pd.read_csv(os.path.join(a.stab_dir, "cpg_freq.tsv"), sep="\t", usecols=["TxID", "CpGID", "saved_selected"])
c = c[c["saved_selected"] == 1].copy()
c["pos"] = c["CpGID"].str.split("_").str[1].astype(int)
r = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t")
r = r[r["stable"]]
m = c.merge(r[["TxID", "start", "end"]], on="TxID")
m = m[(m["pos"] >= m["start"]) & (m["pos"] <= m["end"])].drop_duplicates(["TxID", "CpGID"]).sort_values(["TxID", "pos"])
for tx, x in m.groupby("TxID"):
    with open(os.path.join(a.out_dir, f"{a.prefix}_{tx}"), "w") as fh:
        fh.write("".join(f"{tx}_{cg}\n" for cg in x["CpGID"]))
print(f"{m['TxID'].nunique():,} models, {len(m):,} stable selected CpGs ({m['CpGID'].nunique():,} distinct) -> {a.out_dir}")
