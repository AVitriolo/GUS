#!/usr/bin/env python
# =============================================================================
# gimme_known_enrichment.py
# Known-motif enrichment, foreground vs custom background, with the same
# gimmemotifs calls as `gimme motifs --known -b <bg> -s 0`
# (commands/motifs.py): scan_to_file(score_table=True) -> best z-scored,
# GC-normalised score per motif and sequence; calc_stats_iterator -> stats.
# Differences: the non-redundant motif selection (RFE with cross-validated
# logistic regressions, very slow with ~10^6 sequences) and the HTML report are
# skipped; BH-adjusted p-values across motifs are added.
# Regions are used as given (no resizing). Score tables are reused if present.
# =============================================================================
import argparse, os, logging
import numpy as np
import pandas as pd
from gimmemotifs.scanner import scan_to_file
from gimmemotifs.stats import calc_stats_iterator
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import pfmfile_location

ap = argparse.ArgumentParser()
ap.add_argument("--fg", required=True)
ap.add_argument("--bg", required=True)
ap.add_argument("--genome", default="GRCh38.p14")
ap.add_argument("--pfm", default="gimme.vertebrate.v5.0")
ap.add_argument("--outdir", required=True)
ap.add_argument("--ncpus", type=int, default=12)
a = ap.parse_args()
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("enrichment")
os.makedirs(a.outdir, exist_ok=True)
pfm = pfmfile_location(a.pfm)
log.info(f"motifs: {pfm}")

tables = {"fg": os.path.join(a.outdir, "input.motif.score.txt"),
          "bg": os.path.join(a.outdir, "background.motif.score.txt")}
for key, infile in (("fg", a.fg), ("bg", a.bg)):
    if os.path.exists(tables[key]):
        log.info(f"reusing {tables[key]}"); continue
    log.info(f"scanning {infile}")
    tmp = tables[key] + ".part"
    scan_to_file(infile, pfm, filepath_or_buffer=tmp, score_table=True,
                 genome=a.genome, zscore=True, gcnorm=True, ncpus=a.ncpus)
    os.rename(tmp, tables[key])

n_fg = pd.read_csv(tables["fg"], comment="#", sep="\t", usecols=[0]).shape[0]
n_bg = pd.read_csv(tables["bg"], comment="#", sep="\t", usecols=[0]).shape[0]
log.info(f"sequences: {n_fg} foreground, {n_bg} background; calculating stats")

stats = ["phyper_at_fpr", "roc_auc", "pr_auc", "enr_at_fpr", "recall_at_fdr", "matches_at_fpr"]
rows = []
for motif_stats in calc_stats_iterator(motifs=pfm, fg_table=tables["fg"], bg_table=tables["bg"],
                                       stats=stats, ncpus=a.ncpus):
    for m, s in motif_stats.items():
        rows.append({"Motif": m,
                     "# matches": s["matches_at_fpr"][0], "% matches input": 100 * s["matches_at_fpr"][0] / n_fg,
                     "# matches background": s["matches_at_fpr"][1], "% matches background": 100 * s["matches_at_fpr"][1] / n_bg,
                     "P-value": s["phyper_at_fpr"], "ROC AUC": s["roc_auc"], "PR AUC": s["pr_auc"],
                     "Enr. at 1% FPR": s["enr_at_fpr"], "Recall at 10% FDR": s["recall_at_fdr"]})
df = pd.DataFrame(rows)

# BH across motifs
p = df["P-value"].to_numpy(); o = np.argsort(p); r = np.empty_like(o); r[o] = np.arange(1, len(p) + 1)
q = p * len(p) / r
q_sorted = np.minimum.accumulate(q[o][::-1])[::-1]
df["q-value (BH)"] = np.minimum(q_sorted[np.argsort(o)], 1)
df["log10 P-value"] = -np.log10(np.maximum(df["P-value"], 1e-300))

# motif -> factors, from the database annotation (direct / indirect or predicted)
direct, indirect = {}, {}
for m in read_motifs(pfm):
    f = m.factors if hasattr(m.factors, "keys") else {}
    d = sorted(set(f.get("direct", [])))
    direct[str(m)] = ",".join(d)                    # report names are str(motif): <id>_<consensus>
    indirect[str(m)] = ",".join(sorted({x for k, v in f.items() if k != "direct" for x in v} - set(d)))
df["factors_direct"] = df["Motif"].map(direct)
df["factors_indirect"] = df["Motif"].map(indirect)

df = df.sort_values(["q-value (BH)", "ROC AUC"], ascending=[True, False])
out = os.path.join(a.outdir, "gimme.roc.report.txt")
df.to_csv(out, sep="\t", index=False, float_format="%.4g")
log.info(f"{len(df)} motifs -> {out} | q < 0.05 and enrichment > 1: {((df['q-value (BH)'] < 0.05) & (df['Enr. at 1% FPR'] > 1)).sum()}")
