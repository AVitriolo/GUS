"""
combat_fit_transform.py

Leakage-safe ComBat: fit batch parameters on TRAIN, apply them to TEST.

Standard ComBat implementations (sva::ComBat, inmoose.pycombat_norm) estimate
and apply in a single call, so they cannot be used for train/test discipline.
This module separates the two steps, following the ComBat algorithm
(Johnson, Li & Rabinovic 2007) with parametric empirical-Bayes priors.

VALIDATION: combat_fit + combat_apply on the SAME data must reproduce
pycombat_norm's output. See validate() below.
"""
import numpy as np


# ---------------------------------------------------------------- priors
def _aprior(gamma_hat):
    m = gamma_hat.mean()
    s2 = gamma_hat.var(ddof=0)
    return (2 * s2 + m ** 2) / s2


def _bprior(gamma_hat):
    m = gamma_hat.mean()
    s2 = gamma_hat.var(ddof=0)
    return (m * s2 + m ** 3) / s2


def _postmean(g_hat, g_bar, n, d_star, t2):
    return (t2 * n * g_hat + d_star * g_bar) / (t2 * n + d_star)


def _postvar(sum2, n, a, b):
    return (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)


def _it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=1e-8):
    """Iterative empirical-Bayes solution (sva::it.sol)."""
    n = (1 - np.isnan(sdat)).sum(axis=1).astype(float)
    g_old = g_hat.copy()
    d_old = d_hat.copy()
    change = 1.0
    it = 0
    while change > conv and it < 1000:
        g_new = _postmean(g_hat, g_bar, n, d_old, t2)
        sum2 = np.nansum((sdat - g_new[:, None]) ** 2, axis=1)
        d_new = _postvar(sum2, n, a, b)
        change = max(
            np.nanmax(np.abs(g_new - g_old) / np.abs(g_old)),
            np.nanmax(np.abs(d_new - d_old) / np.abs(d_old)),
        )
        g_old, d_old = g_new, d_new
        it += 1
    return g_old, d_old


# ---------------------------------------------------------------- fit
def combat_fit(X, batch, covar_mod=None, par_prior=True, ref_batch=None):
    """
    Estimate ComBat parameters from TRAINING data.

    X         : (genes x samples) ndarray
    batch     : (samples,) array of batch labels
    covar_mod : (samples x k) design of covariates to PRESERVE, or None
    ref_batch : label of a batch to use as reference, or None
    """
    X = np.asarray(X, dtype=float)
    batch = np.asarray(batch)
    batches = np.unique(batch)
    n_batch = len(batches)
    n_array = X.shape[1]

    batch_idx = [np.where(batch == b)[0] for b in batches]
    n_batches = np.array([len(i) for i in batch_idx], dtype=float)

    # design = batch indicators (+ covariates to preserve)
    batch_design = np.zeros((n_array, n_batch))
    for i, idx in enumerate(batch_idx):
        batch_design[idx, i] = 1.0
    design = batch_design if covar_mod is None else np.hstack([batch_design, covar_mod])

    # --- standardize ---
    B_hat = np.linalg.lstsq(design, X.T, rcond=None)[0]      # (p x genes)

    if ref_batch is None:
        grand_mean = (n_batches / n_array) @ B_hat[:n_batch, :]
    else:
        r = int(np.where(batches == ref_batch)[0][0])
        grand_mean = B_hat[r, :]

    resid = X - (design @ B_hat).T
    var_pooled = (resid ** 2).sum(axis=1) / float(n_array)

    stand_mean = np.outer(grand_mean, np.ones(n_array))
    if covar_mod is not None:
        # covariate contribution is PRESERVED: add it back into the standardisation mean
        tmp = design.copy()
        tmp[:, :n_batch] = 0.0
        stand_mean = stand_mean + (tmp @ B_hat).T

    sd = np.sqrt(var_pooled)[:, None]
    sd[sd == 0] = np.nan
    s_data = (X - stand_mean) / sd

    # --- batch effect parameters ---
    gamma_hat = np.linalg.lstsq(batch_design, s_data.T, rcond=None)[0]   # (n_batch x genes)
    # NOTE: ddof=0 — matches sva/inmoose (numpy .var default)
    delta_hat = np.vstack([np.nanvar(s_data[:, idx], axis=1, ddof=0) for idx in batch_idx])

    if par_prior:
        gamma_bar = np.nanmean(gamma_hat, axis=1)
        t2 = np.nanvar(gamma_hat, axis=1, ddof=0)   # ddof=0 to match sva/inmoose
        a_prior = np.array([_aprior(d[~np.isnan(d)]) for d in delta_hat])
        b_prior = np.array([_bprior(d[~np.isnan(d)]) for d in delta_hat])

        gamma_star, delta_star = [], []
        for i in range(n_batch):
            g, d = _it_sol(
                s_data[:, batch_idx[i]], gamma_hat[i], delta_hat[i],
                gamma_bar[i], t2[i], a_prior[i], b_prior[i],
            )
            gamma_star.append(g)
            delta_star.append(d)
        gamma_star = np.vstack(gamma_star)
        delta_star = np.vstack(delta_star)
    else:
        gamma_star, delta_star = gamma_hat, delta_hat

    # if a reference batch is used, its own effect is not adjusted
    if ref_batch is not None:
        r = int(np.where(batches == ref_batch)[0][0])
        gamma_star[r, :] = 0.0
        delta_star[r, :] = 1.0

    return {
        "batches": batches,
        "grand_mean": grand_mean,
        "var_pooled": var_pooled,
        "gamma_star": gamma_star,
        "delta_star": delta_star,
        "n_batch": n_batch,
        "ref_batch": ref_batch,
    }


# ---------------------------------------------------------------- apply
def combat_apply(X, batch, fit, covar_mod=None):
    """
    Apply TRAIN-fitted ComBat parameters to data (train or test).

    covar_mod must be the covariate design for THESE samples (same columns
    as used at fit time), or None.
    """
    X = np.asarray(X, dtype=float)
    batch = np.asarray(batch)
    n_array = X.shape[1]

    stand_mean = np.outer(fit["grand_mean"], np.ones(n_array))
    if covar_mod is not None:
        # NOTE: preserving covariates on new data requires the fitted covariate
        # coefficients; supported only when covar_mod is None in this version.
        raise NotImplementedError(
            "covar_mod at apply-time not supported; use covar_mod=None (mod=NULL)."
        )

    sd = np.sqrt(fit["var_pooled"])[:, None]
    sd[sd == 0] = np.nan
    s_data = (X - stand_mean) / sd

    bayesdata = s_data.copy()
    for i, b in enumerate(fit["batches"]):
        idx = np.where(batch == b)[0]
        if len(idx) == 0:
            continue
        bayesdata[:, idx] = (
            (bayesdata[:, idx] - fit["gamma_star"][i][:, None])
            / np.sqrt(fit["delta_star"][i])[:, None]
        )

    return bayesdata * sd + stand_mean


# ---------------------------------------------------------------- validation
def validate(seed=0, n_genes=500, n_per_batch=(60, 40), verbose=True):
    """
    Validate against inmoose.pycombat_norm: fitting and applying on the SAME
    data must reproduce pycombat_norm's corrected matrix.
    """
    import pandas as pd
    from inmoose.pycombat import pycombat_norm

    rng = np.random.default_rng(seed)
    nA, nB = n_per_batch
    n = nA + nB

    base = rng.normal(7, 2, size=(n_genes, 1))
    X = base + rng.normal(0, 1, size=(n_genes, n))
    # inject a batch effect into batch B (additive + multiplicative)
    X[:, nA:] += rng.normal(1.5, 0.4, size=(n_genes, 1))
    X[:, nA:] *= rng.normal(1.3, 0.1, size=(n_genes, 1))

    batch = np.array(["A"] * nA + ["B"] * nB)
    df = pd.DataFrame(X, index=[f"g{i}" for i in range(n_genes)],
                      columns=[f"s{i}" for i in range(n)])

    ref = np.asarray(pycombat_norm(df, list(batch)))

    fit = combat_fit(X, batch, covar_mod=None, par_prior=True)
    mine = combat_apply(X, batch, fit, covar_mod=None)

    diff = np.abs(mine - ref)
    max_abs = np.nanmax(diff)
    rel = max_abs / np.ptp(ref)          # relative to the data's dynamic range
    corr = np.corrcoef(mine.ravel(), ref.ravel())[0, 1]

    if verbose:
        print(f"  genes={n_genes}  samples={n}  batches=A:{nA}/B:{nB}")
        print(f"  max |diff|      : {max_abs:.3e}")
        print(f"  relative to range: {rel:.3e}")
        print(f"  correlation     : {corr:.12f}")
    return rel, corr


if __name__ == "__main__":
    print("Validating combat_fit + combat_apply against inmoose.pycombat_norm")
    print("(fit and apply on the SAME data must reproduce pycombat_norm)\n")
    ok = True
    for seed in range(3):
        print(f"[seed {seed}]")
        rel, corr = validate(seed=seed)
        # tolerance: EB iteration converges to ~1e-6 abs; require 1e-5 relative
        passed = (rel < 1e-5) and (corr > 1 - 1e-12)
        print(f"  -> {'PASS' if passed else 'FAIL'}\n")
        ok &= passed
    print("OVERALL:", "PASS - matches ComBat to EB convergence tolerance"
          if ok else "FAIL - do NOT use")