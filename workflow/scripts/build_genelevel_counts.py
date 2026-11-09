"""Gene-level expression input (config level: gene).

  representative TxID of a gene = its selected TxID (ff_selected_TxIDs, the transcript-level model list) with
                                  the highest mean in the transcript-level ComBat-corrected TRAIN counts (ties:
                                  the first in the counts-file order; the original hand-made build broke the
                                  4 exact ties of v38 / 25_10 the other way: ANKRD30A, AP1B1, HCK, OR51F2)
  gene counts                   = raw kallisto counts of ALL the gene's transcripts summed (tximport principle;
                                  no transcript lengths, so no lengthScaledTPM)
  rows                          = the representative TxID for genes that have one, else GENE_<symbol> (non-word
                                  characters -> _), so the files plug into filter_counts_NTM_train_test.R and the
                                  per-TxID rules unchanged

Usage: python build_genelevel_counts.py --counts_NT F --counts_M F --conversion_table F --selected_TxIDs F
           --counts_train_tx F --out_NT F --out_M F --out_representative F
"""
import argparse
import re

import pandas as pd

p = argparse.ArgumentParser()
for a in ("counts_NT", "counts_M", "conversion_table", "selected_TxIDs", "counts_train_tx",
          "out_NT", "out_M", "out_representative"):
    p.add_argument(f"--{a}", required=True)
a = p.parse_args()

conv = pd.read_csv(a.conversion_table, sep=" ", usecols=["transcript_id", "gene_name"])
tx2gene = dict(zip(conv.transcript_id, conv.gene_name))
selected = set(open(a.selected_TxIDs).read().split())

mean_train = pd.read_csv(a.counts_train_tx, sep="\t", index_col=0).mean(axis=1)
cand = mean_train[mean_train.index.isin(selected)].to_frame("mean_train")
cand["gene"] = cand.index.map(tx2gene)
cand = cand.dropna(subset=["gene"])
# ties (transcripts with identical counts) -> the first one in the counts-file order
rep = cand.loc[cand.groupby("gene").mean_train.idxmax()].rename_axis("TxID").reset_index()[["TxID", "gene", "mean_train"]]
rep.sort_values("TxID").to_csv(a.out_representative, sep="\t", index=False)
print(f"representative TxIDs: {len(rep):,} genes")
row_of = dict(zip(rep.gene, rep.TxID))

for src, out in ((a.counts_NT, a.out_NT), (a.counts_M, a.out_M)):
    x = pd.read_csv(src, sep="\t", index_col=0)
    g = x.index.map(tx2gene)
    x = x[g.notna()].groupby(g[g.notna()]).sum()
    x.index = [row_of.get(s, "GENE_" + re.sub(r"\W", "_", s)) for s in x.index]
    x.index.name = "transcript_id"
    x.to_csv(out, sep="\t")
    print(f"{src}: {len(x):,} genes ({x.index.str.startswith('ENST').sum():,} named by their representative TxID) -> {out}")
