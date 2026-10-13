"""PyDESeq2 anchor (slow: ~1 min per cohort at N=5, ~4 min at N=10). Runs the paper-style pipeline -- Wald test with
formal threshold H0: |log2FC| <= 1 (lfc_null=1, greaterAbs), BH FDR < 5 %, covariates tumour purity + age -- on the
first K cohorts of selected sizes, and compares with the fast engine on the *same* cohorts.
Reference for precision/recall: the authors' full-data DESeq2 result (data/upstream_BSLA.deseq2.lfc1.csv).
Resumable: finished cohorts are stored in results/pydeseq2_anchor_sets.json and skipped on re-run.
Usage: python workflow/pydeseq2_crosscheck.py config/config.yaml --sizes 3 5 7 10 15 --k 10
Output: results/pydeseq2_anchor.csv (one row per cohort and engine), results/pydeseq2_anchor_pairs.csv,
        results/pydeseq2_anchor_sets.json (DEG sets per cohort and engine)
"""
import argparse, json, os, sys, time, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, yaml
import fast_de
from cohorts import plan

warnings.filterwarnings("ignore")


def pydeseq2_call(counts, meta, cfg, n_cpus):
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    cts = counts.loc[:, counts.sum() >= cfg["min_total_count"]]
    cov = meta[cfg["covariates"]]
    md = ((cov - cov.mean()) / cov.std(ddof=0)).copy()
    cond = cfg["condition_col"]
    md[cond] = meta[cond].values
    design = "~" + "+".join(cfg["covariates"] + [cond])
    dds = DeseqDataSet(counts=cts, metadata=md, design=design, quiet=True, n_cpus=n_cpus)
    dds.deseq2()
    st = DeseqStats(dds, contrast=[cond, cfg["group_b"], cfg["group_a"]], lfc_null=cfg["primary_lfc"],
                    alt_hypothesis="greaterAbs", quiet=True, n_cpus=n_cpus)
    st.summary()
    r = st.results_df
    return set(r.index[(r.padj < cfg["primary_fdr"])])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("config")
    ap.add_argument("--sizes", type=int, nargs="+", default=[5, 10]); ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--cpus", type=int, default=2)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    counts = pd.read_csv(cfg["counts"], index_col=0).T
    meta = pd.read_csv(cfg["meta"], index_col=0).loc[counts.index]
    cond = cfg["condition_col"]
    ref_f, ref_dft = fast_de.fit(counts, (meta[cond] == cfg["group_b"]).to_numpy(), meta[cfg["covariates"]], cfg["min_total_count"])
    ref_deg = set(ref_f.index[fast_de.call(ref_f, ref_dft, cfg["primary_lfc"], cfg["primary_fdr"]).deg])
    au = pd.read_csv(cfg["authors_reference"], index_col=0); auth = set(au.index[au.FDR < cfg["authors_reference_fdr"]])
    sp = "results/pydeseq2_anchor_sets.json"
    store = json.load(open(sp)) if os.path.exists(sp) else {}
    t0 = time.time()
    todo = sorted([c for c in plan(cfg) if c["N"] in a.sizes and c["trial"] < a.k], key=lambda c: (c["trial"], c["N"]))
    for c in todo:
        key = f"{c['N']}_{c['trial']}"
        if key in store:
            continue
        ids = c["samples"]; m = meta.loc[ids]
        f, dft = fast_de.fit(counts.loc[ids], (m[cond] == cfg["group_b"]).to_numpy(), m[cfg["covariates"]], cfg["min_total_count"])
        fc = fast_de.call(f, dft, cfg["primary_lfc"], cfg["primary_fdr"])
        d_fast = sorted(fc.index[fc.deg]); d_py = sorted(pydeseq2_call(counts.loc[ids], m, cfg, a.cpus))
        store[key] = {"fast": d_fast, "pydeseq2": d_py}
        json.dump(store, open(sp, "w"))
        print(f"N={c['N']} trial={c['trial']} fast={len(d_fast)} pydeseq2={len(d_py)}  {time.time()-t0:.0f}s", flush=True)
    rows, pr = [], []
    for key, v in store.items():
        n, t = map(int, key.split("_"))
        for eng in ("fast", "pydeseq2"):
            d = set(v[eng])
            rows.append(dict(N=n, trial=t, engine=eng, n_deg=len(d),
                             precision_authors=len(d & auth) / len(d) if d else np.nan, recall_authors=len(d & auth) / len(auth),
                             precision_own_ref=len(d & ref_deg) / len(d) if d else np.nan, recall_own_ref=len(d & ref_deg) / len(ref_deg)))
    pd.DataFrame(rows).sort_values(["N", "trial", "engine"]).to_csv("results/pydeseq2_anchor.csv", index=False)
    for n in sorted({int(k.split("_")[0]) for k in store}):
        ts = sorted(int(k.split("_")[1]) for k in store if int(k.split("_")[0]) == n)
        for eng in ("fast", "pydeseq2"):
            ss = [set(store[f"{n}_{t}"][eng]) for t in ts]
            for i in range(len(ss)):
                for j in range(i + 1, len(ss)):
                    u = len(ss[i] | ss[j]); pr.append(dict(N=n, engine=eng, jaccard=len(ss[i] & ss[j]) / u if u else 0.0))
    pd.DataFrame(pr).to_csv("results/pydeseq2_anchor_pairs.csv", index=False)
    print("done", flush=True)
