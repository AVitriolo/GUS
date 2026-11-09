#!/usr/bin/env python
# =============================================================================
# nonredundant_motifs.py
# Groups redundant motifs with the criterion of gimme's own non-redundant
# selection (gimmemotifs.comparison.select_nonredundant_motifs: MotifComparer,
# "partial" / "seqcor" / "mean", similarity >= 0.7), greedily from the most
# enriched motif down: the first motif of a group is its representative.
# Unlike gimme, all motifs are grouped (not only those with >= 2-fold
# enrichment) and the RFE step is not run.
# Motifs longer than --max_len are excluded (they cannot fit in the windows;
# gimme then reports spurious hits with very negative scores).
# =============================================================================
import argparse, logging
import pandas as pd
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import pfmfile_location
from gimmemotifs.comparison import MotifComparer

ap = argparse.ArgumentParser()
ap.add_argument("--report", required=True)          # gimme.roc.report.txt (gimme_known_enrichment.py)
ap.add_argument("--pfm", default="gimme.vertebrate.v5.0")
ap.add_argument("--max_len", type=int, default=20)
ap.add_argument("--min_similarity", type=float, default=0.7)
ap.add_argument("--output", required=True)
a = ap.parse_args()
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("nonredundant")

motifs = {m.id: m for m in read_motifs(pfmfile_location(a.pfm))}
rep = pd.read_csv(a.report, sep="\t")
rep["motif_id"] = rep["Motif"].str.rsplit("_", n=1).str[0]
rep["length"] = rep["motif_id"].map(lambda i: len(motifs[i]))
long_ = rep["length"] > a.max_len
log.info(f"{len(rep)} motifs | excluded, longer than {a.max_len} bp: {long_.sum()}")

order = rep[~long_].sort_values(["Enr. at 1% FPR", "q-value (BH)"], ascending=[False, True])["motif_id"].tolist()
mc = MotifComparer()
cluster, sim = {}, {}
remaining = order[:]
while remaining:
    r = remaining[0]
    scores = mc.get_all_scores([motifs[r]], [motifs[m] for m in remaining], "partial", "seqcor", "mean")[r]
    members = [m for m in remaining if m == r or scores.get(m, [0])[0] >= a.min_similarity]
    for m in members:
        cluster[m], sim[m] = r, (1.0 if m == r else scores[m][0])
    remaining = [m for m in remaining if m not in set(members)]
    if len(set(cluster.values())) % 50 == 0:
        log.info(f"{len(set(cluster.values()))} groups, {len(remaining)} motifs left")

rep["representative"] = rep["motif_id"].map(cluster)
rep["similarity_to_rep"] = rep["motif_id"].map(sim)
rep["is_representative"] = rep["motif_id"] == rep["representative"]
rep["group_size"] = rep["representative"].map(rep["representative"].value_counts())
rep["excluded_too_long"] = long_
cols = ["Motif", "motif_id", "length", "representative", "is_representative", "similarity_to_rep", "group_size",
        "excluded_too_long", "Enr. at 1% FPR", "q-value (BH)", "ROC AUC", "factors_direct", "factors_indirect"]
rep[cols].to_csv(a.output, sep="\t", index=False, float_format="%.4g")
sig = (rep["q-value (BH)"] < 0.05) & (rep["Enr. at 1% FPR"] > 1)
log.info(f"groups: {rep['is_representative'].sum()} from {(~long_).sum()} motifs | "
         f"enriched (q < 0.05): {sig.sum()} motifs -> {(sig & rep['is_representative']).sum()} non-redundant, "
         f"{(sig & long_).sum()} too long")
