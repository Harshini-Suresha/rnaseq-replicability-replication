"""Assemble docs/index.html: a self-contained static dashboard (no server, works on GitHub Pages).
Reads results/*, data/published_degen_medo_2025_results.csv and report/validation.json, writes the data as JSON
into workflow/dashboard_template.html -> docs/index.html. Safe to re-run while the PyDESeq2 anchor is still filling.
Usage: python workflow/build_dashboard.py config/config.yaml"""
import json, os, sys
import numpy as np, pandas as pd, yaml

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
r4 = lambda x: None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), 4)


def qs(x, qlist):
    x = pd.Series(x).dropna()
    return [r4(x.quantile(q)) for q in qlist] if len(x) else [None] * len(qlist)


mc = pd.read_csv("results/metrics_by_cohort.csv")
rp = pd.read_csv("results/replicability_by_N.csv")
fast = {}
for (l, q), g in mc.groupby(["lfc", "fdr"]):
    key = f"{l}|{q}"
    fast[key] = {}
    for n, h in g.groupby("N"):
        rr = rp[(rp.lfc == l) & (rp.fdr == q) & (rp.N == n)].iloc[0]
        e = dict(n_deg=qs(h.n_deg, [.25, .5, .75]), precision=qs(h.precision, [.25, .5, .75]),
                 recall=qs(h.recall, [.25, .5, .75]), spearman=qs(h.spearman_truth, [.25, .5, .75]),
                 rep=[r4(rr.rep_q05), r4(rr.rep_q25), r4(rr.rep_median), r4(rr.rep_q75), r4(rr.rep_q95)],
                 frac_no_calls=r4(h.n_deg.eq(0).mean()), n_ref=int(h.n_ref.iloc[0]))
        if "precision_authors" in h and h.precision_authors.notna().any():
            e["precision_authors"] = qs(h.precision_authors, [.25, .5, .75])
            e["recall_authors"] = qs(h.recall_authors, [.25, .5, .75])
        fast[key][int(n)] = e

anchor = {}
sp = "results/pydeseq2_anchor_sets.json"
if os.path.exists(sp):
    st = json.load(open(sp))
    au = pd.read_csv(cfg["authors_reference"], index_col=0); auth = set(au.index[au.FDR < cfg["authors_reference_fdr"]])
    rf = pd.read_csv("results/reference_full_data.csv.gz", index_col=0)
    own = set(rf.index[rf[f"deg_lfc{float(cfg['primary_lfc'])}_fdr{float(cfg['primary_fdr'])}"] == 1])
    for n in sorted({int(k.split("_")[0]) for k in st}):
        ts = sorted(int(k.split("_")[1]) for k in st if int(k.split("_")[0]) == n)
        anchor[n] = {"k": len(ts)}
        for eng in ("pydeseq2", "fast"):
            ss = [set(st[f"{n}_{t}"][eng]) for t in ts]
            pj = []; pj10 = []
            for i in range(len(ss)):
                for j in range(i + 1, len(ss)):
                    u = len(ss[i] | ss[j]); v_ = r4(len(ss[i] & ss[j]) / u) if u else 0.0; pj.append(v_)
                    if min(len(ss[i]), len(ss[j])) >= 10: pj10.append(v_)
            anchor[n][eng] = dict(n_deg=[len(s) for s in ss],
                                  precision_authors=[r4(len(s & auth) / len(s)) if s else None for s in ss],
                                  recall_authors=[r4(len(s & auth) / len(auth)) for s in ss],
                                  precision_own=[r4(len(s & own) / len(s)) if s else None for s in ss],
                                  recall_own=[r4(len(s & own) / len(own)) for s in ss], rep=pj, rep10=pj10)

pub = pd.read_csv("data/published_degen_medo_2025_results.csv", index_col=0)
published = [dict(code=c, rep5=r4(r.Rep_N5), rep10=r4(r.Rep_N10), prec5=r4(r.Prec_N5), prec10=r4(r.Prec_N10),
                  rec5=r4(r.Rec_N5), rec10=r4(r.Rec_N10), sp5=r4(r.Spear_Truth_N5_mean), sp10=r4(r.Spear_Truth_N10_mean))
             for c, r in pub.iterrows()]

rec = pd.read_csv("results/gene_recurrence.csv.gz", index_col=0)
rd = rec[rec.full_data_deg == 1]
ncols = [c for c in rec.columns if c.startswith("N")]
annot = pd.read_csv("data/gene_annotation.csv", index_col=0) if os.path.exists("data/gene_annotation.csv") else pd.DataFrame(columns=["symbol"])
def sym(g):
    v = annot["symbol"].get(g) if g in annot.index else None
    return v if isinstance(v, str) and v else g
genes = dict(id=list(rd.index), sym=[sym(g) for g in rd.index], lfc=[r4(x) for x in rd.full_data_log2FC], auth=[int(x) for x in rd.authors_deseq2_deg],
             rec={c[1:]: [int(round(x * 100)) for x in rd[c]] for c in ncols})


# ---------- sample embeddings, single genes, per-cohort distributions (for the box / violin / scatter panels) ----------
from scipy.stats import gaussian_kde
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from cohorts import plan


def box_stats(v):
    v = np.asarray(v, float); q1, med, q3 = np.percentile(v, [25, 50, 75]); iqr = q3 - q1
    lo = v[v >= q1 - 1.5 * iqr].min(); hi = v[v <= q3 + 1.5 * iqr].max()
    return [r4(lo), r4(q1), r4(med), r4(q3), r4(hi)]


def density(v, npts=36):
    v = np.asarray(v, float)
    if len(v) < 3 or np.ptp(v) < 1e-9: return None
    try:
        k = gaussian_kde(v)
    except Exception:
        return None
    g = np.linspace(v.min(), v.max(), npts); d = k(g)
    return [[r4(a), r4(b / d.max())] for a, b in zip(g, d)]


emb = pd.read_csv("results/embeddings.csv"); pv = pd.read_csv("results/pca_variance.csv")
est = json.load(open("results/embedding_stats.json"))
sidx = {s_: i for i, s_ in enumerate(emb["sample"])}
cl = {}
for c in plan(cfg):
    if c["trial"] < 20: cl.setdefault(str(c["N"]), {})[str(c["trial"])] = [sidx[x] for x in c["samples"]]
D_emb = dict(sample=list(emb["sample"]), cond=[0 if c == cfg["group_a"] else 1 for c in emb.condition], groups=[cfg["group_a"], cfg["group_b"]],
             purity=[r4(x) for x in emb.purity], age_years=[r4(-x / 365.25) for x in emb.age_days],
             pca=[[r4(a), r4(b)] for a, b in zip(emb.PC1, emb.PC2)], tsne=[[r4(a), r4(b)] for a, b in zip(emb.tsne1, emb.tsne2)],
             umap=[[r4(a), r4(b)] for a, b in zip(emb.umap1, emb.umap2)],
             pc_var=[r4(pv.explained_variance_ratio[0]), r4(pv.explained_variance_ratio[1])], stats=est, cohorts=cl)

tg = pd.read_csv("results/top_gene_expression.csv", index_col=0)
refg = pd.read_csv("results/reference_full_data.csv.gz", index_col=0)
gcols = [c for c in tg.columns if c != "condition"]
cond = tg.condition.values
D_gene = dict(ids=gcols, sym=[sym(g) for g in gcols], lfc=[r4(refg.loc[g, "log2FoldChange"]) for g in gcols], groups=[cfg["group_a"], cfg["group_b"]], per={})
for g in gcols:
    D_gene["per"][g] = [dict(dots=[r4(x) for x in tg.loc[cond == grp, g]], box=box_stats(tg.loc[cond == grp, g]),
                             dens=density(tg.loc[cond == grp, g])) for grp in D_gene["groups"]]

pj = pd.read_csv("results/pair_jaccard.csv.gz")
rngd = np.random.default_rng(3)
Lp, Pp = float(cfg["primary_lfc"]), float(cfg["primary_fdr"])
mcp = mc[(mc.lfc == Lp) & (mc.fdr == Pp)]
dist = {k: {} for k in ("rep", "rep10", "precision", "recall", "n_deg")}
over = {k: {} for k in dist}
for n_, A_ in anchor.items():
    P_ = A_["pydeseq2"]
    over["rep"][str(n_)] = P_["rep"]; over["rep10"][str(n_)] = P_["rep10"]
    over["precision"][str(n_)] = [x for x in P_["precision_authors"] if x is not None]
    over["recall"][str(n_)] = P_["recall_authors"]; over["n_deg"][str(n_)] = P_["n_deg"]
for n in cfg["cohort_sizes"]:
    sets = dict(rep=pj[pj.N == n].jaccard.to_numpy(), rep10=pj[(pj.N == n) & (pj.min_genes >= 10)].jaccard.to_numpy(), precision=mcp[mcp.N == n].precision_authors.dropna().to_numpy(),
                recall=mcp[mcp.N == n].recall_authors.to_numpy(), n_deg=mcp[mcp.N == n].n_deg.to_numpy())
    for k, v in sets.items():
        if len(v) == 0: continue
        dots = rngd.choice(v, 100, replace=False) if len(v) > 100 else v
        dist[k][str(n)] = dict(dots=[r4(x) for x in dots], box=box_stats(v), dens=density(v), n=int(len(v)))

val = json.load(open("report/validation.json")) if os.path.exists("report/validation.json") else {}
D = dict(meta=dict(sizes=cfg["cohort_sizes"], trials=cfg["trials"], primary=[float(cfg["primary_lfc"]), float(cfg["primary_fdr"])],
                   lfc_grid=cfg["lfc_grid"], fdr_grid=cfg["fdr_grid"], n_genes=int(len(rec)), n_per_group=126,
                   n_ref_own=int(rd.shape[0]), n_ref_authors=int(rec.authors_deseq2_deg.sum())),
         fast=fast, anchor=anchor, published=published, genes=genes, validation=val, emb=D_emb, gexpr=D_gene, dist=dist, over=over)
tpl = open("workflow/dashboard_template.html").read()
html = tpl.replace("/*__DATA__*/null", json.dumps(D, separators=(",", ":")))
os.makedirs("docs", exist_ok=True)
open("docs/index.html", "w").write(html)
print(f"docs/index.html written ({len(html)/1024:.0f} KB); anchor sizes: {sorted(anchor)}")
