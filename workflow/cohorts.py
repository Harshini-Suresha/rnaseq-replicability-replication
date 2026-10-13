"""Deterministic random-cohort planner (paper design: independent random cohorts, sampled without replacement)."""
import numpy as np, pandas as pd


def plan(cfg):
    """One cohort per (N, trial): N samples from each group drawn without replacement from the full dataset.
    Cohorts at the same N are independent draws and may share samples (as in Degen & Medo 2025). The seed depends
    only on (seed, N, trial), so adding sizes or trials never changes existing cohorts."""
    meta = pd.read_csv(cfg["meta"], index_col=0)
    ga = meta.index[meta[cfg["condition_col"]] == cfg["group_a"]].to_numpy()
    gb = meta.index[meta[cfg["condition_col"]] == cfg["group_b"]].to_numpy()
    out = []
    for n in cfg["cohort_sizes"]:
        for t in range(cfg["trials"]):
            rng = np.random.default_rng([cfg["seed"], n, t])
            out.append({"N": int(n), "trial": t,
                        "samples": list(rng.choice(ga, n, replace=False)) + list(rng.choice(gb, n, replace=False))})
    return out
