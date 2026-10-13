"""Generate report/results_tables.md from results/* so every number in the write-up comes from the pipeline output.
Usage: python workflow/make_report.py config/config.yaml"""
import json, os, sys
import numpy as np, pandas as pd, yaml

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
L, P = float(cfg["primary_lfc"]), float(cfg["primary_fdr"])
mc = pd.read_csv("results/metrics_by_cohort.csv"); rp = pd.read_csv("results/replicability_by_N.csv")
mc = mc[(mc.lfc == L) & (mc.fdr == P)]; rp = rp[(rp.lfc == L) & (rp.fdr == P)].set_index("N")
pub = pd.read_csv("data/published_degen_medo_2025_results.csv", index_col=0).loc["BSLA"]
f = lambda x: "–" if pd.isna(x) else f"{x:.2f}"
out = []

out += ["## Table 1. Fast engine, 100 random cohorts per size (primary setting: |log2FC| > 1, 5 % FDR)\n",
        "Precision and recall are against the authors' full-data DESeq2 result; medians over cohorts. "
        "Precision is undefined for cohorts that call no genes and is computed over the rest.\n",
        "| N per group | Median genes called | Cohorts with no calls | Precision | Recall | Replicability (median Jaccard) | Log2FC rank agreement with full data |",
        "|---|---|---|---|---|---|---|"]
for n, h in mc.groupby("N"):
    out.append(f"| {n} | {h.n_deg.median():.0f} | {(h.n_deg == 0).mean() * 100:.0f} % | {f(h.precision_authors.median())} | "
               f"{f(h.recall_authors.median())} | {f(rp.loc[n, 'rep_median'])} | {f(h.spearman_truth.median())} |")

sp = "results/pydeseq2_anchor_sets.json"
if os.path.exists(sp):
    st = json.load(open(sp))
    au = pd.read_csv(cfg["authors_reference"], index_col=0); auth = set(au.index[au.FDR < cfg["authors_reference_fdr"]])
    out += ["", "## Table 2. PyDESeq2 anchor (covariates, formal |log2FC| > 1 test, 5 % FDR) vs the fast engine on the same cohorts\n",
            "Medians over the cohorts run; replicability is the median Jaccard over all pairs of those cohorts. "
            "Reference for precision and recall: the authors' full-data DESeq2 result.\n",
            "| N per group | Cohorts | PyDESeq2 genes called | PyDESeq2 precision | PyDESeq2 recall | PyDESeq2 replicability | Fast engine replicability, same cohorts | Fast engine precision, same cohorts |",
            "|---|---|---|---|---|---|---|---|"]
    for n in sorted({int(k.split("_")[0]) for k in st}):
        ts = sorted(int(k.split("_")[1]) for k in st if int(k.split("_")[0]) == n)
        res = {}
        for eng in ("pydeseq2", "fast"):
            ss = [set(st[f"{n}_{t}"][eng]) for t in ts]
            pj = [len(ss[i] & ss[j]) / len(ss[i] | ss[j]) if (ss[i] | ss[j]) else 0.0 for i in range(len(ss)) for j in range(i + 1, len(ss))]
            pr = [len(s & auth) / len(s) for s in ss if s]
            rc = [len(s & auth) / len(auth) for s in ss]
            res[eng] = dict(nd=np.median([len(s) for s in ss]), pj=np.median(pj) if pj else np.nan,
                            pr=np.median(pr) if pr else np.nan, rc=np.median(rc))
        out.append(f"| {n} | {len(ts)} | {res['pydeseq2']['nd']:.0f} | {f(res['pydeseq2']['pr'])} | {f(res['pydeseq2']['rc'])} | "
                   f"{f(res['pydeseq2']['pj'])} | {f(res['fast']['pj'])} | {f(res['fast']['pr'])} |")

out += ["", "## Table 3. Authors' published values for this data set (their results table, row BSLA)\n",
        "| N per group | Precision | Recall | Replicability | Log2FC rank agreement with full data (mean) |", "|---|---|---|---|---|",
        f"| 5 | {f(pub.Prec_N5)} | {f(pub.Rec_N5)} | {f(pub.Rep_N5)} | {f(pub.Spear_Truth_N5_mean)} |",
        f"| 10 | {f(pub.Prec_N10)} | {f(pub.Rec_N10)} | {f(pub.Rep_N10)} | {f(pub.Spear_Truth_N10_mean)} |"]

rows = pd.read_csv("results/metrics_by_cohort.csv")
rows = rows[rows.N == 10]
sens = rows.groupby(["lfc", "fdr"]).agg(n_ref=("n_ref", "first"), n_deg=("n_deg", "median"), prec=("precision", "median"), rec=("recall", "median")).reset_index()
rr = pd.read_csv("results/replicability_by_N.csv"); rr = rr[rr.N == 10].set_index(["lfc", "fdr"]).rep_median
out += ["", "## Table 4. Sensitivity of the N = 10 result to the calling threshold (fast engine; reference = this repo's full-data call at the same setting)\n",
        "| Fold-change threshold | FDR | Reference genes | Median genes called | Precision | Recall | Replicability |", "|---|---|---|---|---|---|---|"]
for _, r in sens.iterrows():
    out.append(f"| {'none' if r.lfc == 0 else '|log2FC| > ' + str(r.lfc)} | {r.fdr * 100:.0f} % | {r.n_ref:.0f} | {r.n_deg:.0f} | {f(r.prec)} | {f(r.rec)} | {f(rr.loc[(r.lfc, r.fdr)])} |")
open("report/results_tables.md", "w").write("\n".join(out) + "\n")
print("\n".join(out))
