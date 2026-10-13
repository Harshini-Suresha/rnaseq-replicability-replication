"""Main analysis (fast engine). For each cohort size N: 100 random cohorts -> DE -> DEG calls under each
(lfc, fdr) setting -> metrics against the full-data reference.

Definitions (following Degen & Medo 2025, PLOS Comput Biol):
  reference ("truth")  the same pipeline run on all 252 samples
  precision / recall   |DEG_cohort & DEG_ref| / |DEG_cohort|  and  / |DEG_ref|   (precision undefined if no calls)
  replicability        Jaccard index between the DEG sets of two cohorts, over all cohort pairs at that N
                       (0 if either set is empty); summarised by the median
  spearman_truth       Spearman correlation of cohort logFC with full-data logFC
Outputs (results/): metrics_by_cohort.csv, replicability_by_N.csv, pair_jaccard.csv.gz, gene_recurrence.csv.gz, reference_full_data.csv.gz
Usage: python workflow/run_sweep.py config/config.yaml [--workers 2]
"""
import argparse, os, sys, time, warnings
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd, yaml
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fast_de
from cohorts import plan

warnings.filterwarnings("ignore")
G = {}


def init(cfg_path):
    cfg = yaml.safe_load(open(cfg_path))
    G["cfg"] = cfg
    G["counts"] = pd.read_csv(cfg["counts"], index_col=0).T
    G["meta"] = pd.read_csv(cfg["meta"], index_col=0).loc[G["counts"].index]
    f, dft = fit_cohort(list(G["counts"].index))
    G["ref_fit"] = f
    G["ref"] = {(l, q): fast_de.call(f, dft, l, q) for l in cfg["lfc_grid"] for q in cfg["fdr_grid"]}
    a = pd.read_csv(cfg["authors_reference"], index_col=0)
    G["authors_deg"] = set(a.index[a.FDR < cfg["authors_reference_fdr"]])


def fit_cohort(ids):
    cfg = G["cfg"]
    m = G["meta"].loc[ids]
    isb = (m[cfg["condition_col"]] == cfg["group_b"]).to_numpy()
    cov = m[cfg["covariates"]] if cfg.get("covariates") else None
    return fast_de.fit(G["counts"].loc[ids], isb, cov, cfg["min_total_count"])


def work(item):
    cfg = G["cfg"]
    f, dft = fit_cohort(item["samples"])
    cm = f.index.intersection(G["ref_fit"].index)
    sp = spearmanr(f.loc[cm].log2FoldChange, G["ref_fit"].loc[cm].log2FoldChange)[0]
    rows, degsets = [], {}
    for (l, q), ref in G["ref"].items():
        c = fast_de.call(f, dft, l, q)
        d = set(c.index[c.deg]); r = set(ref.index[ref.deg])
        tp = len(d & r)
        row = dict(N=item["N"], trial=item["trial"], lfc=l, fdr=q, n_deg=len(d), n_ref=len(r),
                   precision=tp / len(d) if d else np.nan, recall=tp / len(r) if r else np.nan, spearman_truth=sp)
        if (l, q) == (cfg["primary_lfc"], cfg["primary_fdr"]):
            ta = len(d & G["authors_deg"])
            row.update(precision_authors=ta / len(d) if d else np.nan, recall_authors=ta / len(G["authors_deg"]))
        rows.append(row)
        degsets[(l, q)] = d
    return item["N"], item["trial"], rows, degsets


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("config"); ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    init(a.config)
    cfg = G["cfg"]
    os.makedirs("results", exist_ok=True)
    items = plan(cfg)
    t0 = time.time()
    with ProcessPoolExecutor(a.workers, initializer=init, initargs=(a.config,)) as ex:
        out = list(ex.map(work, items, chunksize=8))
    print(f"fitted {len(out)} cohorts in {time.time() - t0:.0f}s")
    pd.DataFrame([r for o in out for r in o[2]]).to_csv("results/metrics_by_cohort.csv", index=False)
    genes = G["ref_fit"].index
    gidx = {g: i for i, g in enumerate(genes)}
    rep_rows, pair_rows = [], []
    for n in cfg["cohort_sizes"]:
        sets_n = [o[3] for o in out if o[0] == n]
        for k in sets_n[0]:
            X = np.zeros((len(sets_n), len(genes)), dtype=np.float32)
            for i, s in enumerate(sets_n):
                X[i, [gidx[g] for g in s[k]]] = 1
            sz = X.sum(1); inter = X @ X.T; union = sz[:, None] + sz[None, :] - inter
            with np.errstate(invalid="ignore", divide="ignore"):
                J = np.where(union > 0, inter / union, 0.0)
            v = J[np.triu_indices(len(sets_n), 1)]
            rep_rows.append(dict(N=n, lfc=k[0], fdr=k[1], pairs=len(v), rep_median=np.median(v), rep_mean=v.mean(),
                                 rep_q05=np.quantile(v, .05), rep_q25=np.quantile(v, .25),
                                 rep_q75=np.quantile(v, .75), rep_q95=np.quantile(v, .95)))
            if k == (cfg["primary_lfc"], cfg["primary_fdr"]):
                pair_rows.append(pd.DataFrame({"N": n, "jaccard": v.round(4), "min_genes": np.minimum(sz[:, None], sz[None, :])[np.triu_indices(len(sets_n), 1)].astype(int)}))
                if n == cfg["cohort_sizes"][0]:
                    rec = pd.DataFrame(index=genes)
                rec[f"N{n}"] = X.mean(0)
    pd.DataFrame(rep_rows).to_csv("results/replicability_by_N.csv", index=False)
    pd.concat(pair_rows).to_csv("results/pair_jaccard.csv.gz", index=False)  # every cohort pair at the primary setting
    prim = (cfg["primary_lfc"], cfg["primary_fdr"])
    rec.insert(0, "authors_deseq2_deg", [int(g in G["authors_deg"]) for g in genes])
    rec.insert(0, "full_data_deg", G["ref"][prim].deg.astype(int).to_numpy())
    rec.insert(0, "full_data_log2FC", G["ref_fit"].log2FoldChange.round(4).to_numpy())
    rec.round(3).to_csv("results/gene_recurrence.csv.gz")
    ref = pd.DataFrame({"log2FoldChange": G["ref_fit"].log2FoldChange, "se": G["ref_fit"].se})
    for (l, q), r in G["ref"].items():
        ref[f"deg_lfc{l}_fdr{q}"] = r.deg.astype(int)
    ref.round(5).to_csv("results/reference_full_data.csv.gz")
    print("done")
