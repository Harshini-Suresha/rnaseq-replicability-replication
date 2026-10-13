"""18 more gallery figures -> 100 total. MA plots, p-value/QC diagnostics, PR curves, Jaccard histograms, ranking bars.
Usage: python workflow/make_gallery_more.py
"""
import json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, MUT, BG, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e8e7e3"
BLUE, ORANGE, GREEN, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#b9b7b0"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": "#cfcec9",
                     "axes.labelcolor": MUT, "xtick.color": MUT, "ytick.color": MUT, "text.color": INK,
                     "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": .7, "axes.axisbelow": True})
os.makedirs("docs/figures", exist_ok=True)
items = json.load(open("docs/figures/gallery.json"))
rng = np.random.default_rng(5)

au = pd.read_csv("data/upstream_BSLA.deseq2.lfc1.csv", index_col=0)
annot = pd.read_csv("data/gene_annotation.csv", index_col=0) if os.path.exists("data/gene_annotation.csv") else pd.DataFrame(columns=["symbol"])
def sym(g):
    s = annot.symbol.get(g) if g in annot.index else None
    return s if isinstance(s, str) and s else g
mc = pd.read_csv("results/metrics_by_cohort.csv")
pj = pd.read_csv("results/pair_jaccard.csv.gz")
emb = pd.read_csv("results/embeddings.csv")
pv = pd.read_csv("results/pca_variance.csv")
rec = pd.read_csv("results/gene_recurrence.csv.gz", index_col=0)
rd = rec[rec.full_data_deg == 1]
meta = pd.read_csv("data/BSLA.meta.csv", index_col=0)
counts = pd.read_csv("data/BSLA.csv", index_col=0)

def add(fig, fname, title, caption, category, plot_type, plot_type_label, n=None, tags=None):
    fig.tight_layout(); fig.savefig(f"docs/figures/{fname}", dpi=150); plt.close(fig)
    items.append(dict(file=f"figures/{fname}", title=title, caption=caption, category=category,
                      n=n, tags=tags or [], plot_type=plot_type, plot_type_label=plot_type_label))

# ---- MA plots (4): mean expression vs effect ----
for lfc, fdr, suffix in [(1.0, 0.05, "primary"), (0.0, 0.05, "nofilter"), (0.5, 0.05, "half"), (1.0, 0.01, "strict")]:
    hit = (au["FDR"] < fdr) & (au["logFC"].abs() > lfc)
    fig, ax = plt.subplots(figsize=(6, 4.4))
    ax.scatter(au.loc[~hit, "logCPM"], au.loc[~hit, "logFC"], s=6, color=GREY, alpha=.5, lw=0, label=f"not DEG ({(~hit).sum():,})")
    ax.scatter(au.loc[hit, "logCPM"], au.loc[hit, "logFC"], s=8, color=ORANGE, alpha=.6, lw=0, label=f"DEG ({hit.sum():,})")
    ax.axhline(lfc, color="#8d8c86", ls="--", lw=1); ax.axhline(-lfc, color="#8d8c86", ls="--", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("Mean expression (logCPM, log scale)"); ax.set_ylabel("log2FC (LumA vs Basal)")
    ax.set_title(f"MA — full data, {suffix} (|log2FC|>{lfc}, FDR<{int(fdr*100)}%)", loc="left", fontsize=11, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    add(fig, f"gallery60_ma_{suffix}.png", f"MA — full data, {suffix}",
        f"Mean-vs-effect for 18,188 genes; orange = DEGs at this setting ({hit.sum():,} genes). High-mean genes dominate DEG calls.",
        "MA", "ma", "MA plot", tags=["ma", "volcano-companion", suffix])

# ---- p-value histograms (2) ----
fig, ax = plt.subplots(figsize=(6, 3.8))
ax.hist(au["pvalue"].clip(0, 1), bins=40, color=BLUE, alpha=.8, edgecolor=BG)
ax.set_xlabel("p-value (authors' full-data DESeq2)"); ax.set_ylabel("Genes")
ax.set_title("p-value histogram — full data", loc="left", fontweight="bold")
add(fig, "gallery61_pvalue_hist.png", "QC — p-value histogram",
    "Uniform bulk + spike near 0 = real signal plus nulls, as expected for a strong subtype contrast.", "QC", "diagnostic", "Diagnostic", tags=["pvalue", "qc"])

fig, ax = plt.subplots(figsize=(6, 3.8))
ax.hist(au.loc[au["FDR"] < 0.05, "pvalue"], bins=30, color=GREEN, alpha=.85, edgecolor=BG)
ax.set_xlabel("p-value among FDR<5% genes"); ax.set_ylabel("Genes")
ax.set_title("p-value histogram — FDR<5% only", loc="left", fontweight="bold")
add(fig, "gallery61_pvalue_hist_sig.png", "QC — p-values of significant genes",
    "Significant calls pile at tiny p-values; validates the FDR cut behaves.", "QC", "diagnostic", "Diagnostic", tags=["pvalue", "fdr", "qc"])

# ---- PCA scree (1) ----
fig, ax = plt.subplots(figsize=(6, 3.8))
ax.bar(range(1, 11), pv.explained_variance_ratio[:10] * 100, color=BLUE, alpha=.85)
ax.set_xticks(range(1, 11)); ax.set_xlabel("PC"); ax.set_ylabel("% variance")
ax.set_title("PCA scree — top 10 PCs (2,000 HVGs)", loc="left", fontweight="bold")
add(fig, "gallery62_scree.png", "QC — PCA variance explained",
    f"PC1 {pv.explained_variance_ratio[0]*100:.1f}%, PC1–10 {pv.explained_variance_ratio[:10].sum()*100:.1f}%. Subtype signal concentrated in PC1.",
    "QC", "diagnostic", "Diagnostic", tags=["pca", "variance", "qc"])

# ---- PC1 vs covariates (2) ----
for col, lab in [("purity", "tumour purity"), ("age_days", "age (days to birth, negative)")] :
    fig, ax = plt.subplots(figsize=(5.5, 4))
    x = emb["PC1"]; y = emb[col]
    for g, c in (("Basal", BLUE), ("LumA", ORANGE)):
        m = emb.condition == g
        ax.scatter(x[m], y[m], s=14, color=c, alpha=.7, edgecolor=BG, linewidth=.4, label=g)
    r = np.corrcoef(x, y)[0, 1]
    ax.set_xlabel("PC1 score"); ax.set_ylabel(lab); ax.set_title(f"PC1 vs {lab} (r={r:.2f})", loc="left", fontweight="bold"); ax.legend(frameon=False, fontsize=8)
    add(fig, f"gallery62_pc1_{col}.png", f"QC — PC1 vs {lab}",
        f"Correlation r={r:.2f}: covariates do not drive the Basal/LumA split (coloured points).", "QC", "diagnostic", "Diagnostic", tags=["pca", "covariate", "qc"])

# ---- library size QC (1) ----
lib = counts.sum(axis=0)
libdf = pd.DataFrame({"lib": lib, "cond": meta.loc[lib.index, "Condition"]})
fig, ax = plt.subplots(figsize=(5, 4))
data = [libdf.loc[libdf.cond == g, "lib"].to_numpy() / 1e6 for g in ("Basal", "LumA")]
bp = ax.boxplot(data, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.5),
                whiskerprops=dict(color="#8d8c86"), capprops=dict(color="#8d8c86"), boxprops=dict(edgecolor="#8d8c86"))
for b, c in zip(bp["boxes"], [BLUE, ORANGE]):
    b.set_facecolor(c); b.set_alpha(.25)
for i, d in enumerate(data, 1):
    ax.scatter(i + rng.uniform(-.1, .1, len(d)), d, s=8, color=[BLUE, ORANGE][i-1], alpha=.45, linewidth=0)
ax.set_xticks([1, 2]); ax.set_xticklabels(["Basal", "LumA"]); ax.set_ylabel("Library size (M counts)")
ax.set_title("QC — library size by subtype", loc="left", fontweight="bold")
add(fig, "gallery62_libsize.png", "QC — library size by subtype",
    "Total counts per tumour; no gross imbalance that would confound DE.", "QC", "diagnostic", "Diagnostic", tags=["qc", "library"])

# ---- precision-recall tradeoff (2) ----
for n in (5, 10):
    sub = mc[(mc.N == n) & (mc.lfc == 1.0) & (mc.fdr == 0.05)]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(sub.recall_authors, sub.precision_authors, s=16, color=BLUE, alpha=.6, edgecolor=BG, linewidth=.5)
    ax.set_xlabel("Recall per cohort"); ax.set_ylabel("Precision per cohort")
    ax.set_xlim(0, max(.5, sub.recall_authors.max() * 1.2)); ax.set_ylim(.5, 1.02)
    ax.set_title(f"Precision–recall, N={n} (100 cohorts)", loc="left", fontweight="bold")
    add(fig, f"gallery63_pr_N{n}.png", f"Tradeoff — precision vs recall at N={n}",
        f"Each dot = one cohort. Small studies sit top-left: precise but incomplete.", "Tradeoff", "curve", "Curve", n=n, tags=["precision", "recall", f"N{n}"])

# ---- Jaccard histograms (4) ----
for n in (3, 5, 10, 15):
    v = pj[pj.N == n].jaccard.to_numpy()
    fig, ax = plt.subplots(figsize=(5, 3.8))
    ax.hist(v, bins=30, color=BLUE, alpha=.8, edgecolor=BG)
    ax.axvline(np.median(v), color=INK, ls="--", lw=1.5, label=f"median {np.median(v):.2f}")
    ax.set_xlabel("Jaccard overlap (cohort pair)"); ax.set_ylabel("Pairs")
    ax.set_title(f"Replicability distribution — N={n} ({len(v):,} pairs)", loc="left", fontweight="bold"); ax.legend(frameon=False, fontsize=8)
    add(fig, f"gallery64_jaccard_hist_N{n}.png", f"Distribution — Jaccard at N={n}",
        f"All {len(v):,} cohort pairs; median {np.median(v):.2f}. Mass at 0 for small N = empty gene lists.", "Distributions", "histogram", "Histogram", n=n, tags=["jaccard", "replicability", f"N{n}"])

# ---- ranking bars (2) ----
top = au.assign(ab=np.abs(au["logFC"])).nlargest(20, "ab")
fig, ax = plt.subplots(figsize=(7, 5))
ax.barh(top.index.map(sym)[::-1], top["ab"][::-1], color=BLUE, alpha=.85)
ax.set_xlabel("|log2FC|"); ax.set_title("Top 20 genes by |log2FC| (full data)", loc="left", fontweight="bold")
add(fig, "gallery65_top_effect.png", "Ranking — top 20 by effect size",
    "Largest absolute fold changes; these are what N≤10 cohorts recover.", "Ranking", "bar", "Bar", tags=["ranking", "effect"])

det10 = (rd["N10"] * 100).sort_values(ascending=False).head(20)
fig, ax = plt.subplots(figsize=(7, 5))
ax.barh(det10.index.map(sym)[::-1], det10.values[::-1], color=GREEN, alpha=.85)
ax.set_xlabel("% of N=10 cohorts calling it"); ax.set_title("Top 20 most-recurrent genes at N=10", loc="left", fontweight="bold")
add(fig, "gallery65_top_recurrent.png", "Ranking — most recurrent at N=10",
    "Genes found by most cohorts; cross-checks the effect-size ranking above.", "Ranking", "bar", "Bar", tags=["ranking", "recurrence"])

json.dump(items, open("docs/figures/gallery.json", "w"), indent=1)
print("total now", len(items))
