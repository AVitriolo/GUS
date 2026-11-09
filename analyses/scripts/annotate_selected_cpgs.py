#!/usr/bin/env python
# =============================================================================
# annotate_selected_cpgs.py
# Genomic and repeat annotation of the selected CpGs (w20 regions) compared
# with the candidate CpGs the same models could choose from.
#   sets     : selected   = CpGs selected by the models with cv_r2 > mean
#              all_sel    = CpGs selected by any model (table only)
#              candidates = all candidate CpGs of the models with cv_r2 > mean
#   genomic  : one category per CpG, first that applies:
#              promoter (ENCODE PLS) > proximal enhancer (pELS) > distal
#              enhancer (dELS) > exon > intron > intergenic
#              (exon / intron / intergenic from GENCODE, all genes: exon = any
#              exon of any transcript, intron = in a gene but in no exon);
#              non-exclusive overlaps are also reported
#   repeats  : RepeatMasker (UCSC hg38 rmsk): class and family of the repeat
#              overlapping the CpG (nested repeats: the one with the lowest
#              divergence), "no repeat" otherwise
# A CpG is its C (0-based [pos - 1, pos)).
# =============================================================================
import argparse, gzip, os, subprocess, tempfile
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--features", required=True)      # TxID CpGID ... (sel_feats_to_regions_w20.R)
ap.add_argument("--cvr2", required=True)          # TxID cv_r2 kept
ap.add_argument("--candidates", required=True)    # chr <TAB> pos (candidates of the kept transcripts)
ap.add_argument("--ccre", required=True)          # ENCODE cCREs, class in col 6 (PLS / pELS / dELS)
ap.add_argument("--gtf", required=True)
ap.add_argument("--rmsk", required=True)          # UCSC rmsk.txt.gz
ap.add_argument("--out_dir", required=True)
a = ap.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
log = lambda *x: print(*x, flush=True)

# ---- CpG sets ----
f = pd.read_csv(a.features, sep="\t", usecols=["TxID", "CpGID"])
kept = set(pd.read_csv(a.cvr2, sep="\t").query("kept == True").TxID)
def to_df(ids):
    ids = pd.Series(sorted(set(ids)))
    return pd.DataFrame({"CpGID": ids, "chr": ids.str.replace(r"_\d+$", "", regex=True),
                         "pos": ids.str.extract(r"_(\d+)$")[0].astype(np.int64)})
sets = {"selected": to_df(f.CpGID[f.TxID.isin(kept)]), "all_sel": to_df(f.CpGID)}
c = pd.read_csv(a.candidates, sep="\t", header=None, names=["chr", "pos"], dtype={"chr": str, "pos": np.int64})
c["CpGID"] = c.chr + "_" + c.pos.astype(str)
sets["candidates"] = c
for k, v in sets.items(): log(f"{k}: {len(v)} CpGs")

# ---- interval sets (merged per chromosome) and point queries ----
def merged(df):
    out = {}
    for ch, g in df.groupby("chr"):
        g = g.sort_values("start"); s = g.start.to_numpy(); e = g.end.to_numpy()
        e_run = np.maximum.accumulate(e)
        new = np.r_[True, s[1:] > e_run[:-1]]
        grp = np.cumsum(new) - 1
        out[ch] = (s[new], np.maximum.reduceat(e, np.flatnonzero(new)))
    return out
def hits(points, ivs):
    res = np.zeros(len(points), bool)
    for ch, g in points.groupby("chr"):
        if ch not in ivs: continue
        s, e = ivs[ch]; p = g.pos.to_numpy() - 1
        i = np.searchsorted(s, p, side="right") - 1
        res[g.index] = (i >= 0) & (p < e[np.clip(i, 0, None)])
    return res

cc = pd.read_csv(a.ccre, sep="\t", header=None, usecols=[0, 1, 2, 5], names=["chr", "start", "end", "cls"])
ann = {k: merged(cc[cc.cls == k]) for k in ["PLS", "pELS", "dELS"]}
gtf = pd.read_csv(a.gtf, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4], names=["chr", "feature", "start", "end"])
gtf["start"] -= 1                                          # GTF 1-based closed -> BED
ann["exon"] = merged(gtf[gtf.feature == "exon"])
ann["gene"] = merged(gtf[gtf.feature == "gene"])
log("annotations loaded")

order = [("PLS", "Promoter (PLS)"), ("pELS", "Proximal enhancer (pELS)"), ("dELS", "Distal enhancer (dELS)"),
         ("exon", "Exon"), ("intron", "Intron"), ("intergenic", "Intergenic")]
rows, nonex = [], []
for k, d in sets.items():
    d = d.reset_index(drop=True)
    h = {x: hits(d, ann[x]) for x in ["PLS", "pELS", "dELS", "exon", "gene"]}
    h["intron"] = h["gene"] & ~h["exon"]; h["intergenic"] = ~h["gene"]
    cat = np.full(len(d), "", dtype=object)
    for x, lab in order:
        cat[(cat == "") & h[x]] = lab
    d["category"] = cat; sets[k] = d
    for x, lab in order:
        rows.append([k, lab, int((cat == lab).sum()), 100 * (cat == lab).mean()])
        nonex.append([k, lab, int(h[x].sum()), 100 * h[x].mean()])
gen = pd.DataFrame(rows, columns=["set", "category", "n", "pct"])
gen.to_csv(os.path.join(a.out_dir, "genomic_annotation_counts.tsv"), sep="\t", index=False, float_format="%.4g")
pd.DataFrame(nonex, columns=["set", "category", "n_overlapping", "pct_overlapping"]).to_csv(
    os.path.join(a.out_dir, "genomic_annotation_overlaps_nonexclusive.tsv"), sep="\t", index=False, float_format="%.4g")
log(gen.pivot(index="category", columns="set", values="pct").round(2).to_string())

# ---- repeats: bedtools intersect of the CpGs with RepeatMasker ----
tmp = tempfile.mkdtemp(dir=a.out_dir)
rm_bed = os.path.join(tmp, "rmsk.bed")
rm = pd.read_csv(a.rmsk, sep="\t", header=None, usecols=[2, 5, 6, 7, 10, 11, 12],
                 names=["milliDiv", "chr", "start", "end", "repName", "repClass", "repFamily"])
rm.sort_values(["chr", "start"])[["chr", "start", "end", "repName", "repClass", "repFamily", "milliDiv"]].to_csv(rm_bed, sep="\t", header=False, index=False)
del rm
rep_rows, fam_rows = [], []
for k, d in sets.items():
    pb = os.path.join(tmp, f"{k}.bed")
    pd.DataFrame({"chr": d.chr, "s": d.pos - 1, "e": d.pos, "id": np.arange(len(d))}).sort_values(["chr", "s"]).to_csv(pb, sep="\t", header=False, index=False)
    out = subprocess.run(f"bedtools intersect -a {pb} -b {rm_bed} -wa -wb -sorted", shell=True, check=True,
                         capture_output=True, text=True).stdout
    x = pd.read_csv(pd.io.common.StringIO(out), sep="\t", header=None, usecols=[3, 7, 8, 9, 10],
                    names=["id", "repName", "repClass", "repFamily", "milliDiv"]) if out else pd.DataFrame(columns=["id", "repName", "repClass", "repFamily", "milliDiv"])
    x = x.sort_values(["id", "milliDiv"]).drop_duplicates("id")
    d["repClass"] = "No repeat"; d["repFamily"] = "No repeat"; d["repName"] = ""
    d.loc[x.id.to_numpy(), ["repClass", "repFamily", "repName"]] = x[["repClass", "repFamily", "repName"]].to_numpy()
    for cl, n in d.repClass.value_counts().items(): rep_rows.append([k, cl, n, 100 * n / len(d)])
    for (cl, fa), n in d.groupby(["repClass", "repFamily"]).size().items(): fam_rows.append([k, cl, fa, n, 100 * n / len(d)])
    log(f"{k}: {100 * (d.repClass != 'No repeat').mean():.1f}% in a repeat")
subprocess.run(f"rm -r {tmp}", shell=True)
pd.DataFrame(rep_rows, columns=["set", "repClass", "n", "pct"]).to_csv(os.path.join(a.out_dir, "repeat_class_counts.tsv"), sep="\t", index=False, float_format="%.4g")
pd.DataFrame(fam_rows, columns=["set", "repClass", "repFamily", "n", "pct"]).to_csv(os.path.join(a.out_dir, "repeat_family_counts.tsv"), sep="\t", index=False, float_format="%.4g")

sets["selected"][["CpGID", "category", "repClass", "repFamily", "repName"]].to_csv(
    os.path.join(a.out_dir, "selected_cpgs_annotation.tsv.gz"), sep="\t", index=False, compression="gzip")
log("done")
