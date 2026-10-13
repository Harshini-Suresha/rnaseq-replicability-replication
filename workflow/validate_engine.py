"""Validate the fast engine and the PyDESeq2 anchor against the authors' own results shipped in their repo.
Writes report/validation.md and report/validation.json. Usage: python workflow/validate_engine.py config/config.yaml"""
import json, os, sys
import numpy as np, pandas as pd, yaml
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fast_de

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
counts = pd.read_csv(cfg["counts"], index_col=0).T
meta = pd.read_csv(cfg["meta"], index_col=0).loc[counts.index]
cov = meta[cfg["covariates"]]
isb = (meta[cfg["condition_col"]] == cfg["group_b"]).to_numpy()
out, L = {}, ["# Engine validation\n"]

# 1. full data vs authors' DESeq2
full = fast_de.run(counts, isb, cov, cfg["min_total_count"], cfg["primary_fdr"], cfg["primary_lfc"])
up = pd.read_csv(cfg["authors_reference"], index_col=0)
c = full.index.intersection(up.index); ok = up.loc[c, "logFC"].notna()
sp = spearmanr(full.loc[c].log2FoldChange[ok], up.loc[c].logFC[ok])[0]
A = set(full.index[full.deg]); B = set(up.index[up.FDR < cfg["authors_reference_fdr"]])
out["full_vs_authors_deseq2"] = dict(genes=len(c), spearman_logfc=round(sp, 3), deg_engine=len(A), deg_authors=len(B),
                                      shared=len(A & B), jaccard=round(len(A & B) / len(A | B), 3))
L += ["## 1. Full data (126 vs 126): fast engine vs the authors' DESeq2 (`BSLA.deseq2.lfc1.csv`)\n",
      f"Both use covariates (tumour purity, age) and a formal |log2FC| > 1 threshold test at 5 % FDR. "
      f"logFC Spearman **{sp:.3f}** over {len(c)} genes. DEGs: engine **{len(A)}**, authors **{len(B)}**, shared **{len(A & B)}** "
      f"(Jaccard **{out['full_vs_authors_deseq2']['jaccard']}**).\n"]

# 2. authors' N=5 edgeR cohort
c5 = pd.read_csv("data/upstream_BSLA.N5.csv", index_col=0).T
m5 = pd.read_csv("data/upstream_BSLA.N5.meta.csv", index_col=0).loc[c5.index]
e = pd.read_csv("data/upstream_BSLA.N5.edgeR.lfc1.csv", index_col=0)
r5 = fast_de.run(c5, (m5.Condition == cfg["group_b"]).to_numpy(), m5[cfg["covariates"]], cfg["min_total_count"],
                 cfg["primary_fdr"], cfg["primary_lfc"])
c = r5.index.intersection(e.index)
sp5 = spearmanr(r5.loc[c].log2FoldChange, e.loc[c].logFC)[0]
A = set(r5.index[r5.deg]); B = set(e.index[e.FDR < 0.05])
out["n5_vs_authors_edger"] = dict(spearman_logfc=round(sp5, 3), deg_engine=len(A), deg_edger=len(B), shared=len(A & B))
L += ["## 2. The authors' N = 5 example cohort: fast engine vs their edgeR glmTreat result (`BSLA.N5.original.edgeR.lfc1.csv`)\n",
      f"logFC Spearman **{sp5:.3f}**. DEGs: engine **{len(A)}**, edgeR **{len(B)}**, shared **{len(A & B)}**. "
      "The fast engine is much more conservative than edgeR at this cohort size: it calls few genes, almost all of which edgeR also calls.\n"]

# 3. PyDESeq2 anchor
sp_path = "results/pydeseq2_anchor_sets.json"
if os.path.exists(sp_path):
    st = json.load(open(sp_path))
    rows = []
    for key, v in st.items():
        n, t = map(int, key.split("_"))
        f_, p_ = set(v["fast"]), set(v["pydeseq2"])
        rows.append(dict(N=n, trial=t, fast=len(f_), pydeseq2=len(p_), shared=len(f_ & p_),
                         frac_fast_in_pydeseq2=len(f_ & p_) / len(f_) if f_ else np.nan))
    d = pd.DataFrame(rows).sort_values(["N", "trial"])
    g = d.groupby("N").agg(cohorts=("trial", "count"), fast_median=("fast", "median"), pydeseq2_median=("pydeseq2", "median"),
                           frac_fast_in_pydeseq2=("frac_fast_in_pydeseq2", "median"))
    out["anchor_engine_agreement"] = {int(k): {kk: (None if pd.isna(vv) else round(float(vv), 3)) for kk, vv in r.items()}
                                      for k, r in g.iterrows()}
    L += ["## 3. Same cohorts, two engines: fast engine vs PyDESeq2 (covariates + formal |log2FC| > 1 test, 5 % FDR)\n",
          "| N per group | cohorts | median DEGs, fast | median DEGs, PyDESeq2 | median share of fast calls also called by PyDESeq2 |",
          "|---|---|---|---|---|"]
    L += [f"| {int(n)} | {int(r.cohorts)} | {r.fast_median:.0f} | {r.pydeseq2_median:.0f} | {r.frac_fast_in_pydeseq2:.2f} |" for n, r in g.iterrows()]
    L.append("")
json.dump(out, open("report/validation.json", "w"), indent=1)
open("report/validation.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
