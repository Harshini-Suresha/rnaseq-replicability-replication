"""Extra static figures: PCA / t-SNE / UMAP, and box-with-dots + violin plots (top genes; cohort-level metrics).
Usage: python workflow/make_figures_extra.py config/config.yaml   (run workflow/embeddings.py and run_sweep.py first)"""
import json, os, sys
import numpy as np, pandas as pd, yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
L, P = float(cfg["primary_lfc"]), float(cfg["primary_fdr"])
CA, CB = "#2a78d6", "#eb6834"          # Basal, LumA (series 1 and 2 of the validated palette)
INK, MUT, BG, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e8e7e3"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": "#cfcec9", "axes.labelcolor": MUT,
                     "xtick.color": MUT, "ytick.color": MUT, "text.color": INK, "font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": GRID, "grid.linewidth": .8, "axes.axisbelow": True})
os.makedirs("figures", exist_ok=True)
rng = np.random.default_rng(11)
emb = pd.read_csv("results/embeddings.csv"); st = json.load(open("results/embedding_stats.json"))
annot = pd.read_csv("data/gene_annotation.csv", index_col=0)
def gname(g):
    sym = annot.symbol.get(g)
    return f"{sym}" if isinstance(sym, str) and sym else g
auth = set(pd.read_csv(cfg["authors_reference"], index_col=0).query("FDR < @cfg['authors_reference_fdr']").index)
anchor = {}
if os.path.exists("results/pydeseq2_anchor_sets.json"):
    _st = json.load(open("results/pydeseq2_anchor_sets.json"))
    for _n in sorted({int(k.split("_")[0]) for k in _st}):
        _ts = sorted(int(k.split("_")[1]) for k in _st if int(k.split("_")[0]) == _n)
        _ss = [set(_st[f"{_n}_{t}"]["pydeseq2"]) for t in _ts]
        anchor[_n] = dict(rep=[len(_ss[i] & _ss[j]) / len(_ss[i] | _ss[j]) for i in range(len(_ss)) for j in range(i + 1, len(_ss)) if _ss[i] | _ss[j]],
                          rep_min10=[len(_ss[i] & _ss[j]) / len(_ss[i] | _ss[j]) for i in range(len(_ss)) for j in range(i + 1, len(_ss)) if min(len(_ss[i]), len(_ss[j])) >= 10],
                          precision=[len(x & auth) / len(x) for x in _ss if x], recall=[len(x & auth) / len(auth) for x in _ss], n_deg=[len(x) for x in _ss])
col = emb.condition.map({cfg["group_a"]: CA, cfg["group_b"]: CB})
seq = LinearSegmentedColormap.from_list("blue", ["#cde2fb", "#6da7ec", "#2a78d6", "#184f95", "#0d366b"])
pv = pd.read_csv("results/pca_variance.csv")

# ---- Fig 4: embeddings -------------------------------------------------------------------------------------------
fig, axs = plt.subplots(2, 3, figsize=(14.5, 9), layout="constrained")
views = [("PC1", "PC2", f"PCA", f"PC1 score, {pv.explained_variance_ratio[0]*100:.0f} % of variance (standardised log2 CPM units)", f"PC2 score, {pv.explained_variance_ratio[1]*100:.0f} % (same units)"),
         ("tsne1", "tsne2", "t-SNE", "t-SNE 1 (arbitrary units)", "t-SNE 2 (arbitrary units)"), ("umap1", "umap2", "UMAP", "UMAP 1 (arbitrary units)", "UMAP 2 (arbitrary units)")]
for j, (x, y, t, xl, yl) in enumerate(views):
    ax = axs[0, j]
    for g, c in ((cfg["group_a"], CA), (cfg["group_b"], CB)):
        m = emb.condition == g
        ax.scatter(emb.loc[m, x], emb.loc[m, y], s=26, color=c, alpha=.8, edgecolor=BG, linewidth=.6, label=g)
    ax.set_title(t, loc="left", fontweight="bold"); ax.set_xlabel(xl); ax.set_ylabel(yl)
    if j == 0: ax.legend(frameon=False, loc="best")
    ax = axs[1, j]
    sc = ax.scatter(emb[x], emb[y], s=26, c=emb.purity, cmap=seq, vmin=.2, vmax=1, edgecolor=BG, linewidth=.4)
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(f"{t}, coloured by tumour purity", loc="left", fontsize=10)
cb = fig.colorbar(sc, ax=axs[1, :], shrink=.8, pad=.015, aspect=30); cb.set_label("Tumour purity (fraction)")
fig.suptitle("All 252 tumours: 2,000 most variable genes, PCA to 50 components, then t-SNE and UMAP on the components.\nPCA scores are in standardised log2 CPM units. t-SNE and UMAP coordinates are arbitrary units: only neighbourhoods carry meaning, not distances or axis values.", x=.01, ha="left", fontsize=10.5, color=MUT)
fig.savefig("figures/fig4_pca_tsne_umap.png", dpi=160); plt.close(fig)

# ---- helpers ----------------------------------------------------------------------------------------------------
def box_dots(ax, groups, colors, labels, jitter=.14, dot=14, positions=None, width=.5, alpha_dot=.55):
    pos = positions if positions is not None else list(range(1, len(groups) + 1))
    bp = ax.boxplot(groups, positions=pos, widths=width, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.6),
                    whiskerprops=dict(color="#8d8c86"), capprops=dict(color="#8d8c86"), boxprops=dict(edgecolor="#8d8c86"))
    for b, c in zip(bp["boxes"], colors): b.set_facecolor(c); b.set_alpha(.22)
    for p, g, c in zip(pos, groups, colors):
        ax.scatter(p + rng.uniform(-jitter, jitter, len(g)), g, s=dot, color=c, alpha=alpha_dot, linewidth=0, zorder=3)
    ax.set_xticks(pos); ax.set_xticklabels(labels)


def violin(ax, groups, colors, labels, positions=None, width=.8):
    pos = positions if positions is not None else list(range(1, len(groups) + 1))
    vp = ax.violinplot([np.asarray(g, float) for g in groups], positions=pos, widths=width, showextrema=False)
    for b, c in zip(vp["bodies"], colors): b.set_facecolor(c); b.set_alpha(.35); b.set_edgecolor(c); b.set_linewidth(1.2)
    for p, g in zip(pos, groups):
        q1, med, q3 = np.percentile(g, [25, 50, 75])
        ax.vlines(p, q1, q3, color=INK, lw=4, zorder=4); ax.scatter([p], [med], s=26, color="white", zorder=5, edgecolor=INK, linewidth=1)
    ax.set_xticks(pos); ax.set_xticklabels(labels)


# ---- Fig 5 / 6: top genes -----------------------------------------------------------------------------------------
tg = pd.read_csv("results/top_gene_expression.csv", index_col=0); ref = pd.read_csv("results/reference_full_data.csv.gz", index_col=0)
genes = [g for g in tg.columns if g != "condition"][:6]
for name, fn in (("fig5_top_genes_box_dots", "box"), ("fig6_top_genes_violin", "violin")):
    fig, axs = plt.subplots(2, 3, figsize=(12.5, 7.4))
    for ax, g in zip(axs.flat, genes):
        a = tg.loc[tg.condition == cfg["group_a"], g].to_numpy(); b = tg.loc[tg.condition == cfg["group_b"], g].to_numpy()
        (box_dots if fn == "box" else violin)(ax, [a, b], [CA, CB], [cfg["group_a"], cfg["group_b"]])
        ax.set_title(f"{gname(g)}" + (f"  ({g})" if gname(g) != g else ""), loc="left", fontsize=10, fontweight="bold")
        ax.text(.98, .02, f"log2FC {ref.loc[g, 'log2FoldChange']:+.1f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=MUT)
        ax.set_ylabel("log2 CPM"); ax.grid(axis="x", visible=False)
    fig.suptitle(f"Six reference genes with the largest standardised effect, 126 tumours per group ({'box with every sample' if fn == 'box' else 'violin with median and IQR'}).\nValues at the bottom of a panel are samples with zero counts. Titles show gene symbol and Ensembl ID.",
                 x=.01, ha="left", fontsize=11, color=MUT)
    fig.tight_layout(rect=(0, 0, 1, .94)); fig.savefig(f"figures/{name}.png", dpi=160); plt.close(fig)

# ---- Fig 7 / 8 / 9: cohort-level metrics by N ----------------------------------------------------------------------
mc = pd.read_csv("results/metrics_by_cohort.csv"); mc = mc[(mc.lfc == L) & (mc.fdr == P)]
pj = pd.read_csv("results/pair_jaccard.csv.gz")
sizes = cfg["cohort_sizes"]
PY = "#eb6834"


def draw_metric(ax, data, extra, kind, ylab, title, ylim=None, pair_dots=False, n_labels=None):
    """data: list (per N) of arrays. extra: dict N -> array of PyDESeq2 values (drawn as orange squares). kind: box | violin."""
    pos = list(range(1, len(sizes) + 1))
    labs = [str(n) if not n_labels else f"{n}\n(n={n_labels[i]})" for i, n in enumerate(sizes)]
    if kind == "box":
        bp = ax.boxplot(data, positions=pos, widths=.6, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.6),
                        whiskerprops=dict(color="#8d8c86"), capprops=dict(color="#8d8c86"), boxprops=dict(edgecolor="#8d8c86"))
        for b in bp["boxes"]: b.set_facecolor(CA); b.set_alpha(.22)
        for i, d in enumerate(data, 1):
            dd = np.random.default_rng(3).choice(d, 100, replace=False) if len(d) > 100 else d
            ax.scatter(i - .12 + rng.uniform(-.14, .14, len(dd)), dd, s=9, color=CA, alpha=.35, linewidth=0, zorder=3)
    else:
        violin(ax, data, [CA] * len(sizes), labs, positions=pos, width=.85)
    for i, n in enumerate(sizes, 1):
        if n in extra and len(extra[n]):
            v = np.asarray(extra[n], float)
            ax.scatter(i + .24 + rng.uniform(-.06, .06, len(v)), v, s=22, marker="s", color=PY, edgecolor=BG, linewidth=.6, alpha=.9, zorder=6)
    ax.set_xticks(pos); ax.set_xticklabels(labs); ax.set_title(title, loc="left", fontweight="bold"); ax.set_ylabel(ylab); ax.grid(axis="x", visible=False)
    if ylim: ax.set_ylim(*ylim)


spec = [("Replicability", "Overlap between two cohorts' gene lists\n(Jaccard index, 0 to 1)", lambda n: pj[pj.N == n].jaccard.to_numpy(), lambda n: anchor.get(n, {}).get("rep", []), (0, 1)),
        ("Precision", "Fraction of called genes\nthat are in the reference", lambda n: mc[mc.N == n].precision_authors.dropna().to_numpy(), lambda n: anchor.get(n, {}).get("precision", []), (0, 1.02)),
        ("Recall", "Fraction of reference genes\nthat are called", lambda n: mc[mc.N == n].recall_authors.to_numpy(), lambda n: anchor.get(n, {}).get("recall", []), (0, .6)),
        ("Genes called", "Genes called per cohort (count)", lambda n: mc[mc.N == n].n_deg.to_numpy(), lambda n: anchor.get(n, {}).get("n_deg", []), None)]
from matplotlib.lines import Line2D
for name, kind in (("fig7_cohort_metrics_box_dots", "box"), ("fig8_cohort_metrics_violin", "violin")):
    fig, axs = plt.subplots(2, 2, figsize=(13.5, 8.8), sharex=True)
    for ax, (t, yl, get, ex, yl_lim) in zip(axs.flat, spec):
        draw_metric(ax, [get(n) for n in sizes], {n: ex(n) for n in sizes}, kind, yl, t, yl_lim)
        if ax in axs[1]: ax.set_xlabel("Samples per group (N)")
    note = ("Blue: fast engine (box from all values, dots = up to 100 per size; replicability uses all 4,950 cohort pairs per N, 100 shown). "
            if kind == "box" else "Blue violins: fast engine, with median (white) and IQR (black). Replicability uses all 4,950 cohort pairs per N. ")
    fig.suptitle(f"|log2FC| > 1, 5 % FDR; precision and recall against the authors' full-data DESeq2.\n{note}Orange squares: PyDESeq2 on 10 cohorts per size (45 pairs for replicability).",
                 x=.01, ha="left", fontsize=10.5, color=MUT)
    fig.legend([Line2D([], [], marker="o", ls="", color=CA, alpha=.6), Line2D([], [], marker="s", ls="", color=PY)], ["Fast engine", "PyDESeq2 (anchor)"],
               loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(.5, -.005))
    fig.tight_layout(rect=(0, .03, 1, .93)); fig.savefig(f"figures/{name}.png", dpi=160); plt.close(fig)

# fig 9: replicability restricted to pairs in which both cohorts call at least 10 genes
MING = 10
sub = [pj[(pj.N == n) & (pj.min_genes >= MING)].jaccard.to_numpy() for n in sizes]
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
ext = {n: anchor.get(n, {}).get("rep_min10", []) for n in sizes}
for ax, kind in zip(axs, ("box", "violin")):
    draw_metric(ax, sub, ext, kind, "Overlap between two cohorts' gene lists\n(Jaccard index, 0 to 1)", "Box with dots" if kind == "box" else "Violin", (0, 1), n_labels=[len(v) for v in sub])
    ax.set_xlabel("Samples per group (N), with number of qualifying pairs")
fig.suptitle(f"Replicability restricted to pairs where both cohorts call at least {MING} genes. This removes the artefact of tiny gene lists coinciding (or being empty), at the cost of keeping only the luckier cohorts at small N.\nBlue: fast engine; orange squares: PyDESeq2 pairs meeting the same rule.",
             x=.01, ha="left", fontsize=10.2, color=MUT)
fig.tight_layout(rect=(0, 0, 1, .9)); fig.savefig("figures/fig9_replicability_min10genes.png", dpi=160); plt.close(fig)
print(sorted(os.listdir("figures")))
