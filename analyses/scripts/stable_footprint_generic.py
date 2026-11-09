"""Stable-feature footprint file for any WGBS cohort's stability selection (stability_core.py / prostate scripts output):
every model's selected CpGs lying in a stable co-methylated region (cpg_freq.tsv + region_map_<tag>_shard*.tsv +
region_freq_<tag>.tsv), mapped to the HM450 probes whose footprint contains them (probe_footprint.py); CpGs without a
probe kept as chr_pos. Columns gene, cv_r2, probe: the --first input of analyses/tcga_meteor_pancancer/scripts/saturation_plot.py.
Models are named as their performance files (<gene>_SVR.tsv). --pos_offset shifts the CpG positions before matching
(GBM Wu 2020: +1, its positions are 1 bp lower than the Illumina / prostate convention).
Usage: python stable_footprint_generic.py --stab_dir DIR --perf_dir DIR --out FILE.tsv.gz [--pos_offset 0] [--tag corr0.5]
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
p.add_argument("--perf_dir", required=True)
p.add_argument("--out", required=True)
p.add_argument("--tag", default="corr0.5")
p.add_argument("--pos_offset", type=int, default=0)
p.add_argument("--manifest", default="/home/Prostate_GAS/resources/illumina_probes/HM450.hg38.manifest.tsv.gz")
a = p.parse_args()

cpg = pd.read_csv(os.path.join(a.stab_dir, "cpg_freq.tsv"), sep="\t", usecols=["TxID", "CpGID", "saved_selected"],
                  dtype={"TxID": str})
cpg = cpg[cpg["saved_selected"] == 1]
rm = pd.concat([pd.read_csv(f, sep="\t", header=None, usecols=[0, 2, 3], names=["TxID", "CpGID", "region"], dtype={0: str})
                for f in glob.glob(os.path.join(a.stab_dir, f"region_map_{a.tag}_shard*.tsv"))]).drop_duplicates(["TxID", "CpGID"])
rf = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t", usecols=["TxID", "region", "stable"],
                 dtype={"TxID": str})
st = cpg.merge(rm, on=["TxID", "CpGID"]).merge(rf[rf["stable"]], on=["TxID", "region"])[["TxID", "CpGID"]]

cv = {}
for f in glob.glob(os.path.join(a.perf_dir, "*_SVR.tsv")):
    s = dict(l.rstrip("\n").split("\t")[:2] for l in open(f) if "\t" in l)
    cv[os.path.basename(f)[:-len("_SVR.tsv")]] = float(s["cv_r2"])
st["cv_r2"] = st["TxID"].map(cv)
st = st.dropna(subset=["cv_r2"])
st["key"] = st["CpGID"].str.rsplit("_", n=1).str[0] + "_" + (st["CpGID"].str.rsplit("_", n=1).str[1].astype(int) + a.pos_offset).astype(str)
hit = cpgs_to_probes(st["key"].unique(), load_footprints(a.manifest))
rows = [(g, v, pr) for g, v, k in zip(st["TxID"], st["cv_r2"], st["key"]) for pr in hit.get(k, [k])]
out = pd.DataFrame(rows, columns=["gene", "cv_r2", "probe"]).drop_duplicates()
out.to_csv(a.out, sep="\t", index=False, compression="gzip")
cg = out[out["probe"].str.startswith("cg")]
print(f"{st['TxID'].nunique():,} models, {len(st):,} stable selected CpGs -> {len(out):,} rows; "
      f"{cg['probe'].nunique():,} distinct probes on the array, {(~out['probe'].str.startswith('cg')).sum():,} CpGs off the array")
