"""Checks behind the WGBS vs TCGA comparisons, redone on the pan-cancer runs.

Candidates of a TCGA gene model are rebuilt exactly as in 03_make_gene_datasets.py (all probes of the cohort matrix within
+-WINDOW_BP of the gene's TSS anchor), because the per-gene tables are not kept.

  model_stats.tsv   per cohort: candidate probes and selected probes (k) per model, k / candidates, k / cap
                    (cap = per-fold training size, the k-capped grid limit), for all models and cv R2 >= 0.1 / 0.3
  coherence.tsv     per cohort, WGBS gene models vs the TCGA model of the same gene: % of the WGBS probes (footprint
                    matching) also selected by TCGA, against chance = (% of WGBS probes that are TCGA candidates) *
                    (TCGA k / TCGA candidates); Wilcoxon signed-rank test, as plot_wgbs_vs_tcga_feature_coherence.py.
                    WGBS models: all, cv R2 >= 0.1, cv R2 >= 0.3 (TCGA models: any cv R2)
  core_enrichment.tsv  WGBS probes in the TCGA core (probes selected in every cohort at cv R2 >= threshold), against
                    the core's share of the probes present in every cohort matrix (thresholds applied to both sides)

Usage: python scripts/confirm_analyses.py <WORK_DIR> --wgbs=FILE [--out_suffix=SUFFIX] [COHORT ...]
  FILE: gene, cv_r2, probe TSV (wgbs_prostate_functional_footprint.tsv.gz)
Writes <WORK_DIR>/confirm<SUFFIX>/
"""
import sys, os, glob, math
import numpy as np, pandas as pd
from scipy.stats import wilcoxon

work = sys.argv[1]; SC = os.environ.get("METH_SCALE", "beta")
opt = {a.split("=", 1)[0][2:]: a.split("=", 1)[1] for a in sys.argv[2:] if a.startswith("--")}
cohorts = [a for a in sys.argv[2:] if not a.startswith("--")] or sorted(os.path.basename(os.path.dirname(os.path.dirname(p))) for p in glob.glob(f"{work}/*/results_{SC}/performance"))
PIPE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
conf = dict(l.split("#")[0].strip().split("=", 1) for l in open(f"{PIPE}/config/pipeline.conf") if "=" in l.split("#")[0])
WIN, TEST, CV = int(conf["WINDOW_BP"]), float(conf["TEST_SIZE"]), int(conf["CV"])
out = f"{work}/confirm" + opt.get("out_suffix", ""); os.makedirs(out, exist_ok=True)

# TSS anchors (same rule as 03_make_gene_datasets.py: probe with the smallest |distToTSS| for the gene)
ga = pd.read_csv(f"{PIPE}/resources/HM450.hg19.gene_association.tsv.gz", sep="\t", dtype=str)
cols = ["geneNames", "transcriptTypes", "transcriptIDs", "distToTSS"]
for c in cols: ga[c] = ga[c].str.split(";")
ex = ga.explode(cols).rename(columns={"geneNames": "geneName"})
ex["geneName"] = ex["geneName"].replace({"": pd.NA, "NA": pd.NA}); ex["distToTSS"] = pd.to_numeric(ex["distToTSS"], errors="coerce")
ex = ex.dropna(subset=["geneName", "distToTSS"]); ex["ad"] = ex.distToTSS.abs()
anc = ex.sort_values("ad", kind="stable").drop_duplicates("geneName")
anchor = {g: (c, int(b)) for g, c, b in zip(anc.geneName, anc.CpG_chrm, anc.CpG_beg)}

w = pd.read_csv(opt["wgbs"], sep="\t"); w = w[w.probe.str.startswith("cg")]
wg = {g: (x.cv_r2.iloc[0], set(x.probe)) for g, x in w.groupby("gene")}
print(f"WGBS: {len(wg):,} gene models with >= 1 footprint probe", flush=True)

def rd(f): return float(dict(l.rstrip("\n").split("\t") for l in open(f) if l.strip())["cv_r2"])
stats, coh, avail, sel = [], [], {}, {}
for c in cohorts:
    rng = pd.read_csv(f"{work}/{c}/{c}_data_meth_range.csv", usecols=["seqnames", "start", "probe.id"])
    avail[c] = set(rng["probe.id"])
    byc = {ch: (g.start.values[np.argsort(g.start.values)], g["probe.id"].values[np.argsort(g.start.values)]) for ch, g in rng.groupby("seqnames")}
    n = len(open(f"{work}/{c}/need_gene_exp_{c}.csv").readline().strip().split(",")) - 1
    cap = int((n - math.ceil(TEST * n)) * (CV - 1) / CV)
    M = {}
    for f in glob.glob(f"{work}/{c}/results_{SC}/performance/*_SVR.tsv"):
        g = os.path.basename(f)[:-8]; sf = f"{work}/{c}/results_{SC}/sel_feats/{g}_SVR.tsv"
        if not os.path.exists(sf) or g not in anchor: continue
        ch, pos = anchor[g]
        if ch not in byc: continue
        st, ids = byc[ch]; cand = set(ids[np.searchsorted(st, pos - WIN, "left"):np.searchsorted(st, pos + WIN, "right")])
        s = {l.strip() for l in open(sf) if l.strip()}
        M[g] = (rd(f), cand, s)
    sel[c] = M
    d = pd.DataFrame([(cv, len(ca), len(s), len(s - ca)) for cv, ca, s in M.values()], columns=["cv", "n_cand", "k", "not_cand"])
    for lab, x in [("all", d), ("cv>=0.1", d[d.cv >= .1]), ("cv>=0.3", d[d.cv >= .3])]:
        stats.append(dict(cohort=c, models=lab, n_models=len(x), n_samples=n, k_cap=cap, median_candidates=x.n_cand.median(),
                          median_k=x.k.median(), mean_k=x.k.mean(), median_k_over_candidates=(x.k / x.n_cand).median(),
                          median_k_over_cap=(x.k / cap).median(), selected_not_candidate=int(x.not_cand.sum())))
    for lab, th in [("all", None), ("cv>=0.1", .1), ("cv>=0.3", .3)]:
        o, e = [], []
        for g, (wcv, wp) in wg.items():
            if (th is not None and wcv < th) or g not in M: continue
            _, ca, s = M[g]
            o.append(100 * len(wp & s) / len(wp)); e.append(100 * len(wp & ca) / len(wp) * len(s) / max(1, len(ca)))
        o, e = np.array(o), np.array(e)
        p = wilcoxon(o, e).pvalue if len(o) > 10 and np.any(o != e) else np.nan
        coh.append(dict(cohort=c, wgbs_models=lab, n_comparable=len(o), observed_mean=o.mean(), chance_mean=e.mean(),
                        ratio=o.mean() / e.mean(), wilcoxon_p=p))
    print(f"{c}: {len(M):,} models, cap {cap}, median k {d.k.median():.0f} of {d.n_cand.median():.0f} candidates; "
          f"coherence {coh[-3]['observed_mean']:.1f}% vs chance {coh[-3]['chance_mean']:.1f}% (n={coh[-3]['n_comparable']:,})", flush=True)

S = pd.DataFrame(stats); S.to_csv(f"{out}/model_stats.tsv", sep="\t", index=False, float_format="%.3f")
C = pd.DataFrame(coh); C.to_csv(f"{out}/coherence.tsv", sep="\t", index=False, float_format="%.4g")
uni = set.intersection(*avail.values()); core = []
for th in (0.1, 0.3):
    k = set.intersection(*({p for cv, _, s in sel[c].values() if cv >= th for p in s} for c in cohorts)) & uni
    wp = {p for cv, ps in wg.values() if cv >= th for p in ps} & uni
    sh = len(wp & k)
    core.append(dict(cv_threshold=th, probes_in_every_matrix=len(uni), tcga_core=len(k), core_share_of_probes=100 * len(k) / len(uni),
                     wgbs_probes_in_every_matrix=len(wp), wgbs_in_core=sh, wgbs_share_in_core=100 * sh / max(1, len(wp)),
                     enrichment=(sh / max(1, len(wp))) / (len(k) / len(uni))))
K = pd.DataFrame(core); K.to_csv(f"{out}/core_enrichment.tsv", sep="\t", index=False, float_format="%.3f")
for t in (S, C, K): print(t.to_string(index=False, float_format="%.3g"))
print(f"wrote {out}/model_stats.tsv, coherence.tsv, core_enrichment.tsv")
