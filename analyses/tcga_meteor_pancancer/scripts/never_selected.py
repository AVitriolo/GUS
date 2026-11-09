"""Probes that no gene model ever selected, per cohort and across all cohorts.

Probes available in a cohort = the probes of its preprocessed matrix (<COHORT>_data_meth_range.csv, same rows as the
methylation matrix: probes kept after the REMP imputation). A probe is selected when a gene model with cv R2 >= threshold
has it among its selected features; "any" = every model, whatever its cv R2.

Usage: python scripts/never_selected.py <WORK_DIR> [--out_suffix=SUFFIX] [COHORT ...]     (default: every cohort with results)
Writes <WORK_DIR>/never_selected/never_selected.tsv (counts) and never_selected_probes_<threshold>.txt
(probes available in >= 1 cohort and never selected in any cohort)
"""
import sys, os, glob
import pandas as pd

work = sys.argv[1]; SC = os.environ.get("METH_SCALE", "beta")
cohorts = [a for a in sys.argv[2:] if not a.startswith("--")] or sorted(os.path.basename(os.path.dirname(os.path.dirname(p))) for p in glob.glob(f"{work}/*/results_{SC}/performance"))
TH = (None, 0.1, 0.3)
out = f"{work}/never_selected" + next((a[len("--out_suffix="):] for a in sys.argv[2:] if a.startswith("--out_suffix=")), ""); os.makedirs(out, exist_ok=True)

avail, models = {}, {}
for c in cohorts:
    avail[c] = set(pd.read_csv(f"{work}/{c}/{c}_data_meth_range.csv", usecols=["probe.id"])["probe.id"])
    rows = []
    for f in glob.glob(f"{work}/{c}/results_{SC}/performance/*_SVR.tsv"):
        sf = f"{work}/{c}/results_{SC}/sel_feats/{os.path.basename(f)}"
        if not os.path.exists(sf): continue
        v = dict(l.rstrip("\n").split("\t") for l in open(f) if l.strip())
        rows.append((float(v["cv_r2"]), {l.strip() for l in open(sf) if l.strip()}))
    models[c] = rows
    print(f"{c}: {len(avail[c]):,} probes in the matrix, {len(rows):,} gene models", flush=True)

res = []
union_avail = set().union(*avail.values()); inter_avail = set.intersection(*avail.values())
for th in TH:
    lab = "any" if th is None else str(th)
    sel = {c: set().union(*(s for cv, s in models[c] if th is None or cv >= th)) for c in cohorts}
    for c in cohorts:
        n = len(avail[c] - sel[c])
        res.append(dict(cv_threshold=lab, scope=c, n_probes=len(avail[c]), n_selected=len(avail[c] & sel[c]),
                        n_never_selected=n, pct_never_selected=100 * n / len(avail[c])))
    sel_all = set().union(*sel.values())
    never = union_avail - sel_all
    res.append(dict(cv_threshold=lab, scope="all cohorts (probe in >= 1 matrix)", n_probes=len(union_avail),
                    n_selected=len(union_avail & sel_all), n_never_selected=len(never), pct_never_selected=100 * len(never) / len(union_avail)))
    n2 = len(inter_avail - sel_all)
    res.append(dict(cv_threshold=lab, scope="all cohorts (probe in every matrix)", n_probes=len(inter_avail),
                    n_selected=len(inter_avail & sel_all), n_never_selected=n2, pct_never_selected=100 * n2 / len(inter_avail)))
    with open(f"{out}/never_selected_probes_{lab}.txt", "w") as fh:
        fh.write("\n".join(sorted(never)) + "\n")
R = pd.DataFrame(res); R.to_csv(f"{out}/never_selected.tsv", sep="\t", index=False, float_format="%.2f")
print(R.to_string(index=False, float_format="%.2f"))
print(f"wrote {out}/never_selected.tsv and never_selected_probes_<threshold>.txt")
