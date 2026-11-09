"""Map CpG IDs (chr_pos, 1-based) to the HM450 probes whose 50-nt footprint contains them.

Footprint of a probe = the 50 nt its oligo covers on hg38 (mapPos_A .. mapPos_A + 49) together with the target CpG
(manifest CpG_beg + 1 and the G next to it). A CpG inside the footprints of several probes (probes < 50 bp apart)
maps to all of them.
"""
import numpy as np
import pandas as pd


def load_footprints(manifest):
    m = pd.read_csv(manifest, sep="\t", usecols=["CpG_chrm", "CpG_beg", "Probe_ID", "mapPos_A"]).dropna()
    m = m[m.Probe_ID.str.startswith("cg")].copy()
    t, a = m.CpG_beg.astype(int) + 1, m.mapPos_A.astype(int)
    m["lo"], m["hi"], m["target"] = np.minimum(a, t), np.maximum(a + 49, t + 1), t
    return {c: (g.lo.values, g.hi.values, g.Probe_ID.values) for c, g in m.sort_values("lo").groupby("CpG_chrm")}


def cpgs_to_probes(keys, fp, max_span=60):
    """keys: iterable of 'chr_pos'. Returns {key: [probe IDs whose footprint contains it]} (keys without a probe omitted)."""
    out = {}
    for k in keys:
        c, p = k.rsplit("_", 1)
        if c not in fp:
            continue
        lo, hi, ids = fp[c]; p = int(p)
        i, j = np.searchsorted(lo, p - max_span, "left"), np.searchsorted(lo, p, "right")
        hit = [ids[x] for x in range(i, j) if hi[x] >= p]
        if hit:
            out[k] = hit
    return out
