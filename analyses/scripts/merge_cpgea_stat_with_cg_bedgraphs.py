#!/usr/bin/env python
# Fill the CpGs missing from cpgea_wgbs_with_coverage_hg38 (CPGEA *_mC_Identification.stat:
# only CpGs called significantly methylated, >=5 reads) with the same sample's value from
# cpgea_wgbs_hg38 (CPGEA *_WGBS.CG.bg: all covered CpGs incl. 0%, no read counts).
# Output rows/samples follow cpgea_wgbs_with_coverage_hg38 (R order); filled cells have cov = NA.
# source: 0 = no value, 1 = .stat (with coverage), 2 = WGBS.CG bedGraph (no coverage).
import argparse, time, h5py, numpy as np, pandas as pd
from multiprocessing import Pool

ap = argparse.ArgumentParser()
ap.add_argument("--stat_h5", required=True)
ap.add_argument("--stat_row_index", required=True)
ap.add_argument("--cg_h5", required=True)
ap.add_argument("--out_h5", required=True)
ap.add_argument("--block", type=int, default=200_000)
ap.add_argument("--chunk_cpgs", type=int, default=2_000)
ap.add_argument("--workers", type=int, default=12)
ap.add_argument("--gzip", type=int, default=4)
a = ap.parse_args()

NA_INT = np.iinfo(np.int32).min
R = pd.read_csv(a.stat_row_index).h5_row.to_numpy(np.int64) - 1

def runs(idx):
    br = np.where(np.diff(idx) != 1)[0]
    s = np.r_[0, br + 1]; e = np.r_[br + 1, len(idx)]
    return [(int(i), int(j), int(idx[i]), int(idx[j - 1]) + 1) for i, j in zip(s, e)]

def work(lo_hi):
    lo, hi = lo_hi
    with h5py.File(a.stat_h5, "r") as fs, h5py.File(a.cg_h5, "r") as fc:
        n = fs["assay001"].shape[0]
        b = np.empty((n, hi - lo)); c = np.empty((n, hi - lo), np.int32)
        for i, j, s, e in runs(R[lo:hi]):
            b[:, i:j] = fs["assay001"][:, s:e]; c[:, i:j] = fs["assay002"][:, s:e]
        g = fc["beta"][:, lo:hi]
    src = np.where(np.isfinite(b), 1, 0).astype(np.uint8)
    fill = ~np.isfinite(b) & np.isfinite(g)
    b[fill] = g[fill]; c[fill] = NA_INT; src[fill] = 2
    both = np.isfinite(g) & (src == 1)
    stats = np.array([fill.sum(), (fill & (g == 0)).sum(), both.sum(),
                      (np.abs(b[both] - g[both]) <= 0.005).sum(), (np.isfinite(b) & ~np.isfinite(g)).sum()], np.int64)
    return lo, hi, b, c, src, stats

with h5py.File(a.stat_h5, "r") as fs:
    ns, n = fs["assay001"].shape
with h5py.File(a.cg_h5, "r") as fc:
    assert fc["beta"].shape == (ns, n), fc["beta"].shape
assert len(R) == n
blocks = [(lo, min(lo + a.block, n)) for lo in range(0, n, a.block)]
print(f"{ns} samples x {n:,} CpGs, {len(blocks)} blocks, {a.workers} workers", flush=True)

t0 = time.time(); tot = np.zeros(5, np.int64)
kw = dict(chunks=(ns, a.chunk_cpgs), compression="gzip", compression_opts=a.gzip)
with h5py.File(a.out_h5, "w") as fo:
    db = fo.create_dataset("beta", (ns, n), np.float64, fillvalue=np.nan, **kw)
    dc = fo.create_dataset("cov", (ns, n), np.int32, fillvalue=NA_INT, **kw)
    dsrc = fo.create_dataset("source", (ns, n), np.uint8, fillvalue=0, **kw)
    with Pool(a.workers) as p:
        for k, (lo, hi, b, c, src, st) in enumerate(p.imap(work, blocks), 1):
            db[:, lo:hi] = b; dc[:, lo:hi] = c; dsrc[:, lo:hi] = src; tot += st
            if k % 10 == 0 or k == len(blocks):
                print(f"  {k}/{len(blocks)} blocks, {time.time()-t0:.0f}s", flush=True)
print(f"filled values: {tot[0]:,} (of which 0%: {tot[1]:,}) | values in both: {tot[2]:,}, identical within 0.005: "
      f"{100*tot[3]/max(tot[2],1):.2f}% | values only in .stat: {tot[4]:,} | {time.time()-t0:.0f}s", flush=True)
print("done", flush=True)
