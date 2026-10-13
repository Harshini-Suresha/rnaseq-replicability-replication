"""Fast differential expression for two-group designs with optional covariates (numpy only).

limma-trend-style pipeline: median-of-ratios size factors -> log2 normalised CPM with a small prior count ->
per-gene linear model (covariates + condition) -> empirical-Bayes variance moderation with a mean-dependent prior
(Smyth 2004) -> moderated t. Calls can use a TREAT-style threshold test (H0: |log2FC| <= thr; cf. DESeq2
lfcThreshold, edgeR glmTreat) or an ordinary test followed by a fold-change filter. BH-adjusted p-values.

This is a substitute for the paper's DESeq2/edgeR (R is not available in the build environment); its agreement with
the authors' own results is quantified in report/validation.md.
"""
import numpy as np, pandas as pd
from scipy import stats, special


def size_factors(counts):  # counts: samples x genes
    c = counts.to_numpy(float)
    pos = (c > 0).all(axis=0)
    if pos.sum() < 100:
        lib = c.sum(1)
        return lib / np.exp(np.log(lib).mean())
    logc = np.log(c[:, pos])
    return np.exp(np.median(logc - logc.mean(0), axis=1))


def trigamma_inverse(x):
    x = np.asarray(x, float)
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        dif = tri * (1 - tri / x) / special.polygamma(2, y)
        y = y + dif
        if np.max(-dif / y) < 1e-8:
            break
    return y


def fit_f_dist_trend(s2, df, covariate, deg=3):
    """Prior df d0 and trended prior variance s0sq (method of moments, as limma fitFDist with trend)."""
    z = np.log(np.maximum(s2, 1e-300))
    X = np.vander((covariate - covariate.mean()) / (covariate.std() + 1e-12), deg + 1)
    beta, *_ = np.linalg.lstsq(X, z, rcond=None)
    fit = X @ beta
    shift = special.digamma(df / 2) - np.log(df / 2)
    e = z - fit - shift
    emean, evar = e.mean(), e.var(ddof=1) - special.polygamma(1, df / 2)
    if evar > 0:
        d0 = 2 * trigamma_inverse(evar)
        s0sq = np.exp(emean + shift - (special.digamma(d0 / 2) - np.log(d0 / 2)) + fit)
    else:
        d0, s0sq = np.inf, np.exp(emean + fit + shift)
    return d0, s0sq


def bh(p):
    p = np.asarray(p); n = len(p)
    o = np.argsort(p); ranked = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1)
    return out


def fit(counts, is_b, covars=None, min_total=10, prior=0.5):
    """counts: samples x genes DataFrame. is_b: bool array (True = group B; logFC = B - A).
    covars: optional DataFrame (samples x k) of numeric covariates, z-scored here.
    Returns (DataFrame[log2FoldChange, se], dft)."""
    keep = counts.sum(0) >= min_total
    cts = counts.loc[:, keep]
    sf = size_factors(cts)
    norm = cts.to_numpy(float) / sf[:, None]
    lib = norm.sum(1).mean()
    Y = np.log2((norm + prior) / (lib + 2 * prior) * 1e6)
    n = Y.shape[0]
    cols = [np.ones(n)]
    if covars is not None and covars.shape[1]:
        z = (covars - covars.mean()) / covars.std(ddof=0).replace(0, 1)
        cols += [z.to_numpy(float)[:, j] for j in range(z.shape[1])]
    cols.append(np.asarray(is_b, float))
    X = np.column_stack(cols)
    p = X.shape[1]
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ Y
    R = Y - X @ B
    df = n - p
    s2 = (R ** 2).sum(0) / df
    d0, s0sq = fit_f_dist_trend(s2, df, Y.mean(0))
    s2m = (d0 * s0sq + df * s2) / (d0 + df) if np.isfinite(d0) else s0sq
    se = np.sqrt(s2m * XtXi[-1, -1])
    dft = df + (d0 if np.isfinite(d0) else 1e6)
    return pd.DataFrame({"log2FoldChange": B[-1], "se": se}, index=cts.columns), float(dft)


def call(fitted, dft, lfc_thr=1.0, fdr=0.05, threshold_test=True):
    """DEG calls. threshold_test=True: TREAT-style H0 |logFC|<=thr then BH. False: ordinary test, BH, |logFC|>thr."""
    lfc, se = fitted.log2FoldChange.to_numpy(), fitted.se.to_numpy()
    if threshold_test and lfc_thr > 0:
        tt = (np.abs(lfc) - lfc_thr) / se
        p = np.where(tt > 0, np.minimum(2 * stats.t.sf(tt, dft), 1.0), 1.0)
    else:
        p = 2 * stats.t.sf(np.abs(lfc / se), dft)
    q = bh(p)
    return pd.DataFrame({"log2FoldChange": lfc, "padj": q, "pvalue": p,
                         "deg": (q < fdr) & (np.abs(lfc) > lfc_thr)}, index=fitted.index)


def run(counts, is_b, covars=None, min_total=10, fdr=0.05, lfc_thr=1.0, threshold_test=True):
    f, dft = fit(counts, is_b, covars, min_total)
    return call(f, dft, lfc_thr, fdr, threshold_test)
