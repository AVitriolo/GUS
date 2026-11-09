"""Saturation and coherence with CpG +-50 bp as the unit, across prostate WGBS, GBM WGBS and the TCGA cohorts.

Every selected CpG (WGBS CpGs; for TCGA the probe's target CpG, hg38) is widened to +-50 bp. At each cv R2 threshold the
windows of ALL cohorts are merged into regions (CpGs <= 100 bp apart fall in one region), and a cohort's catalogue is the
set of regions it hits ("regions") or of (gene, region) pairs ("pairs").
  saturation_regions50.{pdf,png}, saturation_regions50_counts.tsv: as saturation_plot.py; the WGBS cohorts stay first
      (prostate, then GBM) and only the TCGA cohorts are shuffled in the random-order panel
  region_sizes.tsv: size of the merged regions
  coherence_regions50.tsv: per TCGA cohort and WGBS cohort, same gene: share of WGBS-selected CpGs that have a
      TCGA-selected probe within 50 bp, among those with a TCGA candidate probe within 50 bp; chance = 1-(1-k/n)^m
      (m = TCGA candidate probes within 50 bp, k/n = TCGA selected/candidates for that gene). TCGA models: any cv R2.
TCGA candidates are rebuilt as in 03_make_gene_datasets.py (probes within +-WINDOW_BP of the TSS anchor).

Usage: python scripts/region_analyses.py <WORK_DIR> --wgbs="Prostate WGBS=FILE:342" --wgbs="GBM WGBS=FILE:64"
           --probes=hm450_probe_hg38.tsv.gz [--out_suffix=SUFFIX] [COHORT ...]
  WGBS FILE: gene, cv_r2, cpg (chr_pos, hg38, 1-based C as Illumina CpG_beg + 1)
Writes <WORK_DIR>/regions50<SUFFIX>/
"""
import sys, os, glob, math
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

work = sys.argv[1]; SC = os.environ.get("METH_SCALE", "beta"); FLANK = 50; TH = (0.1, 0.3); N_PERM = 500
wgbs_args = [a.split("=", 1)[1] for a in sys.argv[2:] if a.startswith("--wgbs=")]
opt = {a.split("=", 1)[0][2:]: a.split("=", 1)[1] for a in sys.argv[2:] if a.startswith("--") and not a.startswith("--wgbs=")}
cohorts = [a for a in sys.argv[2:] if not a.startswith("--")] or sorted(os.path.basename(os.path.dirname(os.path.dirname(p))) for p in glob.glob(f"{work}/*/results_{SC}/performance"))
PIPE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
conf = dict(l.split("#")[0].strip().split("=", 1) for l in open(f"{PIPE}/config/pipeline.conf") if "=" in l.split("#")[0])
WIN = int(conf["WINDOW_BP"])
out = f"{work}/regions50" + opt.get("out_suffix", ""); os.makedirs(out, exist_ok=True)

p2c = pd.read_csv(opt["probes"], sep="\t"); p2c = dict(zip(p2c.Probe_ID, p2c.cpg))
def split(k): c, p = k.rsplit("_", 1); return c, int(p)

# models: cohort -> list of (gene, cv, [(chr, pos)])
models, nsamp, tcga_sel = {}, {}, {}
for a in wgbs_args:
    lab, f = a.split("=", 1)
    f, n = f.rsplit(":", 1) if ":" in f else (f, 0)
    d = pd.read_csv(f, sep="\t"); models[lab] = [(g, x.cv_r2.iloc[0], [split(k) for k in x.cpg]) for g, x in d.groupby("gene")]; nsamp[lab] = int(n)
    print(f"{lab}: {len(models[lab]):,} gene models", flush=True)
WG = list(models)
rd = lambda f: float(dict(l.rstrip("\n").split("\t") for l in open(f) if l.strip())["cv_r2"])
for c in cohorts:
    rows, sel = [], {}
    for f in glob.glob(f"{work}/{c}/results_{SC}/performance/*_SVR.tsv"):
        g = os.path.basename(f)[:-8]; sf = f"{work}/{c}/results_{SC}/sel_feats/{g}_SVR.tsv"
        if not os.path.exists(sf): continue
        ps = [l.strip() for l in open(sf) if l.strip()]; sel[g] = (rd(f), ps)
        rows.append((g, sel[g][0], [split(p2c[p]) for p in ps if p in p2c]))
    models[c] = rows; tcga_sel[c] = sel
    nsamp[c] = len(open(f"{work}/{c}/need_gene_exp_{c}.csv").readline().strip().split(",")) - 1
    print(f"{c}: {len(rows):,} gene models", flush=True)
ALL = WG + cohorts

# ---------------- saturation on merged +-50 bp regions ----------------
def regions(th):
    pts = {}
    for c in ALL:
        for g, cv, cp in models[c]:
            if cv >= th:
                for ch, p in cp: pts.setdefault(ch, []).append(p)
    rid, base, sizes = {}, 0, []
    for ch in sorted(pts):
        a = np.unique(np.array(pts[ch])); brk = np.r_[True, np.diff(a) > 2 * FLANK]; r = np.cumsum(brk) - 1
        rid[ch] = (a, r + base); n = r.max() + 1
        st = a[brk]; en = a[np.r_[brk[1:], True]]; sizes.append(en - st + 2 * FLANK + 1); base += n
    return rid, base, np.concatenate(sizes)
def lookup(rid, ch, p): a, r = rid[ch]; return r[np.searchsorted(a, p)]
order = WG + sorted(cohorts, key=lambda c: -nsamp[c]); ci = {c: k for k, c in enumerate(ALL)}
rng = np.random.default_rng(0); rows, size_rows = [], []
fig, ax = plt.subplots(2, 4, figsize=(20, 9))
gidx = {g: i for i, g in enumerate(sorted({g for c in ALL for g, _, _ in models[c]}))}
for j, (kind, th) in enumerate([(k, t) for k in ("regions", "pairs") for t in TH]):
    rid, NR, sizes = regions(th)
    if kind == "regions": size_rows.append(dict(cv_threshold=th, n_regions=NR, median_bp=np.median(sizes), p90_bp=np.percentile(sizes, 90), max_bp=sizes.max(), share_101bp=np.mean(sizes == 2 * FLANK + 1)))
    ids, bits = [], []
    for c in ALL:
        e = [lookup(rid, ch, p) if kind == "regions" else gidx[g] * NR + lookup(rid, ch, p) for g, cv, cp in models[c] if cv >= th for ch, p in cp]
        e = np.unique(np.array(e, dtype=np.int64)); ids.append(e); bits.append(np.full(len(e), 1 << ci[c], dtype=np.int64))
    ids, bits = np.concatenate(ids), np.concatenate(bits); o = np.argsort(ids, kind="stable"); ids, bits = ids[o], bits[o]
    um, cnt = np.unique(np.bitwise_or.reduceat(bits, np.r_[0, np.flatnonzero(np.diff(ids)) + 1]), return_counts=True)
    def cumulative(ordr):
        rank = np.array([ordr.index(c) for c in ALL]); first = np.full(len(um), len(ALL))
        for k in range(len(ALL)): first = np.where((um >> k) & 1, np.minimum(first, rank[k]), first)
        return np.array([cnt[first <= s].sum() for s in range(len(ALL))]), np.array([cnt[first == s].sum() for s in range(len(ALL))])
    curves = np.array([cumulative(WG + list(rng.permutation(cohorts)))[0] for _ in range(N_PERM)]); x = np.arange(1, len(ALL) + 1)
    a = ax[0, j]; a.fill_between(x, np.percentile(curves, 5, 0), np.percentile(curves, 95, 0), color="#1baf7a", alpha=.25, lw=0, label=f"5–95% of random TCGA orders ({' + '.join(WG)} first)")
    a.plot(x, curves.mean(0), "o-", color="#1f5f6b", ms=4, label="mean"); a.set_xlabel("number of cohorts"); a.set_ylabel(f"distinct functional {kind}")
    a.set_title(f"{kind} (CpG ±{FLANK} bp), cv R² ≥ {th}: random orders", loc="left", fontsize=10); a.legend(frameon=False, fontsize=7)
    cum, new = cumulative(order); own = {c: int(cnt[((um >> ci[c]) & 1).astype(bool)].sum()) for c in ALL}; a = ax[1, j]
    for i, c in enumerate(order):
        a.bar(i, cum[i] - new[i], color="#e1e0d9", width=.7); a.bar(i, new[i], bottom=cum[i] - new[i], color="#1baf7a", width=.7)
        rows.append(dict(catalogue=kind, cv_threshold=th, step=i + 1, cohort=c, n_samples=nsamp[c], cohort_size=own[c], new=int(new[i]), cumulative=int(cum[i]), new_share_of_cohort=new[i] / max(1, own[c])))
    a.set_xticks(range(len(order)), [f"{c}\n({nsamp[c]})" for c in order], fontsize=7, rotation=90)
    a.set_ylabel(f"cumulative functional {kind}"); a.set_title(f"{kind}, cv R² ≥ {th}: WGBS first, then TCGA by sample size (green = new)", loc="left", fontsize=10)
    print(f"saturation {kind} cv>={th}: {NR:,} regions in total", flush=True)
for a in ax.flat: a.spines[["top", "right"]].set_visible(False)
fig.suptitle(f"Saturation of functional CpG regions (CpG ±{FLANK} bp, merged across cohorts): {' + '.join(WG)} + {len(cohorts)} TCGA cohorts", fontsize=12)
plt.tight_layout(); plt.savefig(f"{out}/saturation_regions50.pdf"); plt.savefig(f"{out}/saturation_regions50.png", dpi=90)
pd.DataFrame(rows).to_csv(f"{out}/saturation_regions50_counts.tsv", sep="\t", index=False)
pd.DataFrame(size_rows).to_csv(f"{out}/region_sizes.tsv", sep="\t", index=False)

# ---------------- coherence (same gene), +-50 bp ----------------
ga = pd.read_csv(f"{PIPE}/resources/HM450.hg19.gene_association.tsv.gz", sep="\t", dtype=str)
cols = ["geneNames", "transcriptTypes", "transcriptIDs", "distToTSS"]
for c in cols: ga[c] = ga[c].str.split(";")
ex = ga.explode(cols).rename(columns={"geneNames": "geneName"}); ex["geneName"] = ex["geneName"].replace({"": pd.NA, "NA": pd.NA})
ex["distToTSS"] = pd.to_numeric(ex["distToTSS"], errors="coerce"); ex = ex.dropna(subset=["geneName", "distToTSS"]); ex["ad"] = ex.distToTSS.abs()
anc = ex.sort_values("ad", kind="stable").drop_duplicates("geneName"); anchor = {g: (c, int(b)) for g, c, b in zip(anc.geneName, anc.CpG_chrm, anc.CpG_beg)}
wsel = {w: {g: (cv, cp) for g, cv, cp in models[w]} for w in WG}
coh = []
for c in cohorts:
    rng_ = pd.read_csv(f"{work}/{c}/{c}_data_meth_range.csv", usecols=["seqnames", "start", "probe.id"])
    byc = {ch: (x.start.values[np.argsort(x.start.values)], x["probe.id"].values[np.argsort(x.start.values)]) for ch, x in rng_.groupby("seqnames")}
    per = {w: {lab: [] for lab in ("all", "cv>=0.1", "cv>=0.3")} for w in WG}
    for g, (tcv, ps) in tcga_sel[c].items():
        if g not in anchor or anchor[g][0] not in byc: continue
        ch19, pos19 = anchor[g]; st, ids = byc[ch19]
        cand = ids[np.searchsorted(st, pos19 - WIN, "left"):np.searchsorted(st, pos19 + WIN, "right")]
        cpos = {}; [cpos.setdefault(split(p2c[p])[0], []).append(split(p2c[p])[1]) for p in cand if p in p2c]
        spos = {}; [spos.setdefault(split(p2c[p])[0], []).append(split(p2c[p])[1]) for p in ps if p in p2c]
        n, k = len(cand), len(ps)
        if n == 0: continue
        cpos = {ch: np.sort(v) for ch, v in cpos.items()}; spos = {ch: np.sort(v) for ch, v in spos.items()}
        for w in WG:
            if g not in wsel[w]: continue
            wcv, wcp = wsel[w][g]; obs, exp_ = [], []
            for ch, p in wcp:
                a = cpos.get(ch)
                if a is None: continue
                m = np.searchsorted(a, p + FLANK, "right") - np.searchsorted(a, p - FLANK, "left")
                if m == 0: continue
                s = spos.get(ch, np.array([], int)); sh = np.searchsorted(s, p + FLANK, "right") - np.searchsorted(s, p - FLANK, "left")
                obs.append(sh > 0); exp_.append(1 - (1 - k / n) ** m)
            if not obs: continue
            for lab, th in (("all", None), ("cv>=0.1", .1), ("cv>=0.3", .3)):
                if th is None or wcv >= th: per[w][lab].append((100 * np.mean(obs), 100 * np.mean(exp_), len(obs)))
    for w in WG:
        for lab, v in per[w].items():
            if not v: continue
            o, e, nn = map(np.array, zip(*v)); p = wilcoxon(o, e).pvalue if len(o) > 10 and np.any(o != e) else np.nan
            coh.append(dict(tcga_cohort=c, wgbs_cohort=w, wgbs_models=lab, n_genes=len(o), comparable_cpgs=int(nn.sum()),
                            observed_mean=o.mean(), chance_mean=e.mean(), ratio=o.mean() / e.mean(), wilcoxon_p=p))
    print(f"coherence {c} done", flush=True)
C = pd.DataFrame(coh); C.to_csv(f"{out}/coherence_regions50.tsv", sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(size_rows).to_string(index=False)); print(C.to_string(index=False, float_format="%.3g"))
print(f"wrote {out}/")
