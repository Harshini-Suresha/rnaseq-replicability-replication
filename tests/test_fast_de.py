"""Sanity tests for the fast DE engine. Run: python -m pytest tests/ -q   (or: python tests/test_fast_de.py)"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow"))
import fast_de


def simulate(n=8, genes=3000, n_true=300, lfc=2.0, seed=1):
    rng = np.random.default_rng(seed)
    mu = rng.lognormal(5, 1.5, genes)
    eff = np.zeros(genes); idx = rng.choice(genes, n_true, replace=False)
    eff[idx] = rng.choice([-lfc, lfc], n_true)
    lib = rng.lognormal(0, 0.2, 2 * n)
    isb = np.r_[np.zeros(n, bool), np.ones(n, bool)]
    m = mu[None, :] * lib[:, None] * (2.0 ** (eff[None, :] * isb[:, None]))
    r = 10.0
    counts = rng.negative_binomial(r, r / (r + m))
    return pd.DataFrame(counts, columns=[f"g{i}" for i in range(genes)]), isb, set(f"g{i}" for i in idx)


def test_null_calibration():
    c, isb, _ = simulate(n_true=0)
    res = fast_de.run(c, isb, lfc_thr=0.0, threshold_test=False, fdr=0.05)
    assert res.deg.sum() <= 0.01 * len(res), "too many false positives under the null"


def test_power_and_precision():
    c, isb, truth = simulate(n=8, lfc=2.0)
    res = fast_de.run(c, isb, lfc_thr=1.0, fdr=0.05)
    called = set(res.index[res.deg])
    assert len(called) > 0.5 * len(truth)
    assert len(called & truth) / len(called) > 0.9


def test_sign_convention():
    c, isb, truth = simulate(n=10, lfc=2.0)
    res = fast_de.run(c, isb, lfc_thr=0.0, threshold_test=False)
    top = res.loc[list(truth)]
    assert np.sign(top.log2FoldChange).abs().mean() == 1.0
    # group B has the higher counts when lfc > 0: logFC = B - A
    g = res.index[res.padj < 1e-6]
    assert (res.loc[g, "log2FoldChange"].abs() > 0.5).all()


def test_covariate_runs():
    c, isb, _ = simulate()
    cov = pd.DataFrame({"x": np.random.default_rng(0).normal(size=len(c))})
    res = fast_de.run(c, isb, cov, lfc_thr=1.0)
    assert res.shape[0] > 0 and np.isfinite(res.log2FoldChange).all()


def test_bh_monotone():
    p = np.array([0.001, 0.01, 0.02, 0.5, 0.9])
    q = fast_de.bh(p)
    assert (np.diff(q[np.argsort(p)]) >= -1e-12).all() and (q >= p).all()


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v(); print("ok", k)
