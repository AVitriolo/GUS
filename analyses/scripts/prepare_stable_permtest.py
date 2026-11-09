#!/usr/bin/env python
# =============================================================================
# prepare_stable_permtest.py
# Inputs of permtest_cpg_localization.R for the stable CpGs, from
# stable_regions_annotation.py (cpg_sets_annotation.tsv.gz):
#   universe : distinct candidate CpGs (n_cand sampled per gene with >= 1 stable
#              region) + the stable CpGs themselves (so that every point has a
#              category in the universe; they are candidates too)
#   points   : stable_cpg_points.tsv.gz    every distinct stable CpG
#              stable_region_points.tsv.gz one point per distinct stable region
#                                          (its middle stable CpG)
# Written twice: <out_dir>/genomic (cCRE / gene-structure category) and
# <out_dir>/repeats (RepeatMasker class; "In a repeat" is n - "No repeat").
# =============================================================================
import argparse, os
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--ann_dir", required=True)      # stable_regions_annotation.py out_dir
ap.add_argument("--stab_dir", required=True)     # region_freq_<tag>.tsv
ap.add_argument("--tag", default="corr0.5")
ap.add_argument("--out_dir", required=True)
a = ap.parse_args()

GENOMIC = ["Promoter (PLS)", "Proximal enhancer (pELS)", "Distal enhancer (dELS)", "Exon", "Intron", "Intergenic"]
REPEATS = ["No repeat", "SINE", "LINE", "LTR", "DNA", "Simple_repeat", "Other"]

long = pd.read_csv(os.path.join(a.ann_dir, "cpg_sets_annotation.tsv.gz"), sep="\t")
long = long[long["set"].isin(["candidates", "stable"])]
long["repClass"] = long["repClass"].where(long["repClass"].isin(REPEATS), "Other")
long["chr"] = long["CpGID"].str.rsplit("_", n=1).str[0]
long["pos"] = long["CpGID"].str.rsplit("_", n=1).str[1].astype(np.int64)
u = long.drop_duplicates("CpGID")
st = long[long["set"] == "stable"]
print(f"universe {len(u)} CpGs ({long[long.set == 'candidates'].CpGID.nunique()} candidates); "
      f"stable CpGs {st.CpGID.nunique()} ({st.CpGID.isin(long[long.set == 'candidates'].CpGID).mean():.1%} among sampled candidates)")

# middle stable CpG of each stable region (regions shared by several genes count once)
reg = pd.read_csv(os.path.join(a.stab_dir, f"region_freq_{a.tag}.tsv"), sep="\t").query("stable")
m = st[["TxID", "CpGID", "chr", "pos"]].merge(reg[["TxID", "chr", "start", "end"]], on=["TxID", "chr"])
m = m[(m["pos"] >= m["start"]) & (m["pos"] <= m["end"])].sort_values("pos")
mid = m.groupby(["chr", "start", "end"], as_index=False).apply(lambda g: g.iloc[(len(g) - 1) // 2],
                                                               include_groups=False)
mid = mid.drop_duplicates("CpGID")
print(f"stable regions {len(reg[['chr', 'start', 'end']].drop_duplicates())}; region points {len(mid)}")

for sub, col, cats in [("genomic", "category", GENOMIC), ("repeats", "repClass", REPEATS)]:
    d = os.path.join(a.out_dir, sub)
    os.makedirs(d, exist_ok=True)
    code = {c: i for i, c in enumerate(cats)}
    pd.DataFrame({"cat": range(len(cats)), "category": cats}).to_csv(os.path.join(d, "categories.tsv"), sep="\t", index=False)
    pd.DataFrame({"chr": u["chr"], "pos": u["pos"], "cat": u[col].map(code)}).sort_values(["chr", "pos"]) \
        .to_csv(os.path.join(d, "candidate_categories.tsv.gz"), sep="\t", index=False)
    s = st.drop_duplicates("CpGID")
    pd.DataFrame({"chr": s["chr"], "pos": s["pos"], "cat": s[col].map(code)}).sort_values(["chr", "pos"]) \
        .to_csv(os.path.join(d, "stable_cpg_points.tsv.gz"), sep="\t", index=False)
    mc = mid.merge(u[["CpGID", col]], on="CpGID")
    pd.DataFrame({"chr": mc["chr"], "pos": mc["pos"], "cat": mc[col].map(code)}).sort_values(["chr", "pos"]) \
        .to_csv(os.path.join(d, "stable_region_points.tsv.gz"), sep="\t", index=False)
    print(f"wrote {d}")
