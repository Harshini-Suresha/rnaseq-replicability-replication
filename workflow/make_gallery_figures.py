"""Generate ~70 gallery figures from results/ into docs/figures/ + manifest gallery.json.
Usage: python workflow/make_gallery_figures.py config/config.yaml
Light theme, same palette as dashboard. Total: 3 overview + 10 recurrence + 4 distributions + 40 genes + 9 embeddings + 4 thresholds = 70.
"""
import json, os, sys
import numpy as np, pandas as pd, yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
L, P = float(cfg["primary_lfc"]), float(cfg["primary_fdr"])
sizes = cfg["cohort_sizes"]
BLUE, ORANGE, GREEN = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUT, BG, GRID = "#0b0b0b", "#52514e", "#fcfcfb", "#e8e7e3"
plt.rcParams.update({"figure.facecolor": BG, "axes.facecolor": BG, "axes.edgecolor": "#cfcec9",
                     "axes.labelcolor": MUT, "xtick.color": MUT, "ytick.color": MUT, "text.color": INK,
                     "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": .7, "axes.axisbelow": True})

OUT = "docs/figures"
os.makedirs(OUT, exist_ok=True)
manifest = []
rng = np.random.default_rng(11)

def save(fig, fname, title, caption, category, n=None, tags=None):
    path = os.path.join(OUT, fname)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    manifest.append(dict(file=f"figures/{fname}", title=title, caption=caption,
                         category=category, n=n, tags=tags or []))

def slug_gene(g, annot):
    s = annot.symbol.get(g) if g in annot.index else None
    sym = s if isinstance(s, str) and s else g
    safe = "".join(c if c.isalnum() else "_" for c in sym)[:24]
    return sym, safe

# ---------- load ----------
mc = pd.read_csv("results/metrics_by_cohort.csv")
rp = pd.read_csv("results/replicability_by_N.csv")
rec = pd.read_csv("results/gene_recurrence.csv.gz", index_col=0)
rd = rec[rec.full_data_deg == 1]
tg = pd.read_csv("results/top_gene_expression.csv", index_col=0)
ref = pd.read_csv("results/reference_full_data.csv.gz", index_col=0)
emb = pd.read_csv("results/embeddings.csv")
pv = pd.read_csv("results/pca_variance.csv")
pub = pd.read_csv("data/published_degen_medo_2025_results.csv", index_col=0)
pj = pd.read_csv("results/pair_jaccard.csv.gz")
annot = pd.read_csv("data/gene_annotation.csv", index_col=0) if os.path.exists("data/gene_annotation.csv") else pd.DataFrame(columns=["symbol"])
mcp = mc[(mc.lfc == L) & (mc.fdr == P)]
rpp = rp[(rp.lfc == L) & (rp.fdr == P)].set_index("N")

# ---------- A. overview (3) ----------
# A1 metrics by cohort
fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
for ax, col, ttl, ylab in zip(axs,
    ["rep", "precision", "recall"],
    ["Replicability (median Jaccard)", "Precision (vs authors' full DESeq2)", "Recall (vs authors' full DESeq2)"],
    ["Jaccard overlap", "Precision", "Recall"]):
    x = np.array(sizes)
    if col == "rep":
        lo = [rpp.loc[n, "rep_q25"] for n in sizes]; mid = [rpp.loc[n, "rep_median"] for n in sizes]; hi = [rpp.loc[n, "rep_q75"] for n in sizes]
    elif col == "precision":
        lo = [mcp[mcp.N == n].precision_authors.dropna().quantile(.25) for n in sizes]
        mid = [mcp[mcp.N == n].precision_authors.dropna().median() for n in sizes]
        hi = [mcp[mcp.N == n].precision_authors.dropna().quantile(.75) for n in sizes]
    else:
        lo = [mcp[mcp.N == n].recall_authors.quantile(.25) for n in sizes]
        mid = [mcp[mcp.N == n].recall_authors.median() for n in sizes]
        hi = [mcp[mcp.N == n].recall_authors.quantile(.75) for n in sizes]
    ax.fill_between(x, lo, hi, color=BLUE, alpha=.14, lw=0)
    ax.plot(x, mid, color=BLUE, lw=2, marker="o", ms=5, mec=BG, mew=1.2)
    ax.set_title(ttl, loc="left", fontsize=10, fontweight="bold"); ax.set_xlabel("N per group"); ax.set_ylabel(ylab); ax.set_ylim(0, 1.02); ax.set_xticks(sizes)
save(fig, "gallery01_metrics_by_cohort.png", "Overview — metrics by cohort size",
     "Median + IQR over 100 fast-engine cohorts per size, primary setting (|log2FC|>1, FDR 5%).", "Overview", tags=["metrics", "replicability", "precision", "recall"])

# A2 published 18
fig, ax = plt.subplots(figsize=(7.5, 5.5))
d = pub["Rep_N10"].sort_values()
cols = [GREEN if i == "BSLA" else "#a9dcc8" for i in d.index]
ax.barh(d.index, d.values, color=cols, height=.6)
ax.set_title("Published replicability at N=10 (18 data sets)", loc="left", fontweight="bold")
ax.set_xlabel("Median Jaccard"); ax.set_xlim(0, 1)
save(fig, "gallery02_published_18.png", "Overview — 18 published data sets",
     "Authors' results table, sorted. BSLA (green) is the data set reproduced here.", "Overview", tags=["published"])

# A3 embeddings overview
fig, axs = plt.subplots(1, 3, figsize=(12, 3.8))
for ax, x, y, t in zip(axs, ["PC1", "tsne1", "umap1"], ["PC2", "tsne2", "umap2"], ["PCA", "t-SNE", "UMAP"]):
    for g, c in ((cfg["group_a"], BLUE), (cfg["group_b"], ORANGE)):
        m = emb.condition == g
        ax.scatter(emb.loc[m, x], emb.loc[m, y], s=14, color=c, alpha=.75, edgecolor=BG, linewidth=.4, label=g)
    ax.set_title(t, loc="left", fontweight="bold"); ax.set_xlabel(x); ax.set_ylabel(y)
axs[0].legend(frameon=False, fontsize=8)
save(fig, "gallery03_embeddings_overview.png", "Overview — tumour map (PCA / t-SNE / UMAP)",
     "All 252 tumours from 2,000 most variable genes. Basal (blue) vs LumA (orange).", "Overview", tags=["pca", "tsne", "umap"])

# ---------- B. recurrence per N (10) ----------
for n in sizes:
    fig, ax = plt.subplots(figsize=(5.5, 4))
    col = f"N{n}"
    ax.scatter(rd.full_data_log2FC, rd[col] * 100, s=8, color=BLUE, alpha=.45, lw=0)
    ax.axvline(0, color="#cfcec9", lw=1)
    ax.set_title(f"Gene recurrence at N={n} per group", loc="left", fontweight="bold")
    ax.set_xlabel("Full-data log2FC (LumA vs Basal)"); ax.set_ylabel("% of 100 cohorts calling it"); ax.set_ylim(-3, 103)
    det = int((rd[col] > 0).sum()); half = int((rd[col] >= .5).sum())
    save(fig, f"gallery10_recurrence_N{n}.png", f"Gene recurrence — N={n}",
         f"{det} of {len(rd)} reference genes called by ≥1 cohort; {half} by ≥half. Large |log2FC| survives.",
         "Gene recurrence", n=n, tags=["recurrence", f"N{n}"])

# ---------- C. distributions across N (4) ----------
specs = [("rep", "Replicability — distribution by N", "Jaccard", lambda n: pj[pj.N == n].jaccard.to_numpy(), (0, 1)),
         ("precision", "Precision — distribution by N", "Precision", lambda n: mcp[mcp.N == n].precision_authors.dropna().to_numpy(), (0, 1.02)),
         ("recall", "Recall — distribution by N", "Recall", lambda n: mcp[mcp.N == n].recall_authors.to_numpy(), (0, .6)),
         ("n_deg", "Genes called — distribution by N", "Count", lambda n: mcp[mcp.N == n].n_deg.to_numpy(), None)]
for key, ttl, ylab, get, ylim in specs:
    fig, ax = plt.subplots(figsize=(9, 4))
    data = [get(n) for n in sizes]
    bp = ax.boxplot(data, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.4),
                    whiskerprops=dict(color="#8d8c86"), capprops=dict(color="#8d8c86"), boxprops=dict(edgecolor="#8d8c86"))
    for b in bp["boxes"]:
        b.set_facecolor(BLUE); b.set_alpha(.22)
    for i, d in enumerate(data, 1):
        dd = rng.choice(d, min(100, len(d)), replace=False) if len(d) > 100 else d
        ax.scatter(i + rng.uniform(-.18, .18, len(dd)), dd, s=7, color=BLUE, alpha=.35, linewidth=0, zorder=3)
    ax.set_xticks(range(1, len(sizes) + 1)); ax.set_xticklabels([str(s) for s in sizes])
    ax.set_title(ttl, loc="left", fontweight="bold"); ax.set_xlabel("N per group"); ax.set_ylabel(ylab)
    if ylim: ax.set_ylim(*ylim)
    save(fig, f"gallery20_dist_{key}.png", ttl,
         "Fast engine, 100 cohorts per size (replicability: all 4,950 pairs per N, 100 dots shown). Primary setting.",
         "Distributions", tags=[key])

# ---------- D. top 40 genes (40) ----------
genes40 = [c for c in tg.columns if c != "condition"][:40]
for i, g in enumerate(genes40, 1):
    sym, safe = slug_gene(g, annot)
    a = tg.loc[tg.condition == cfg["group_a"], g].to_numpy()
    b = tg.loc[tg.condition == cfg["group_b"], g].to_numpy()
    lfc = float(ref.loc[g, "log2FoldChange"]) if g in ref.index else float("nan")
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    bp = ax.boxplot([a, b], patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.5),
                    whiskerprops=dict(color="#8d8c86"), capprops=dict(color="#8d8c86"), boxprops=dict(edgecolor="#8d8c86"))
    for box, c in zip(bp["boxes"], [BLUE, ORANGE]):
        box.set_facecolor(c); box.set_alpha(.25)
    for k, (v, c) in enumerate(zip([a, b], [BLUE, ORANGE]), 1):
        ax.scatter(k + rng.uniform(-.12, .12, len(v)), v, s=8, color=c, alpha=.5, linewidth=0, zorder=3)
    ax.set_xticks([1, 2]); ax.set_xticklabels([cfg["group_a"], cfg["group_b"]])
    ax.set_title(f"{sym}  (log2FC {lfc:+.1f})", loc="left", fontsize=10, fontweight="bold")
    ax.set_ylabel("log2 CPM"); ax.grid(axis="x", visible=False)
    save(fig, f"gallery30_gene_{i:02d}_{safe}.png", f"Single gene {i:02d} — {sym}",
         f"{sym} ({g}), log2FC {lfc:+.2f} LumA vs Basal. 126 tumours per group, log2 CPM; zeros at panel bottom.",
         "Single genes", tags=["gene", sym, g])

# ---------- E. embeddings (9) ----------
seq = LinearSegmentedColormap.from_list("blue", ["#cde2fb", "#6da7ec", "#2a78d6", "#184f95", "#0d366b"])
views = [("PC1", "PC2", "pca", "PCA"), ("tsne1", "tsne2", "tsne", "t-SNE"), ("umap1", "umap2", "umap", "UMAP")]
for x, y, key, ttl in views:
    # subtype
    fig, ax = plt.subplots(figsize=(5, 4.2))
    for g, c in ((cfg["group_a"], BLUE), (cfg["group_b"], ORANGE)):
        m = emb.condition == g
        ax.scatter(emb.loc[m, x], emb.loc[m, y], s=18, color=c, alpha=.8, edgecolor=BG, linewidth=.5, label=g)
    ax.set_title(f"{ttl} coloured by subtype", loc="left", fontweight="bold"); ax.set_xlabel(x); ax.set_ylabel(y); ax.legend(frameon=False, fontsize=8)
    save(fig, f"gallery70_embed_{key}_subtype.png", f"Map — {ttl} by subtype",
         "252 tumours; Basal blue vs LumA orange.", "Sample map", tags=[key, "subtype"])
    # purity
    fig, ax = plt.subplots(figsize=(5, 4.2))
    sc = ax.scatter(emb[x], emb[y], s=18, c=emb.purity, cmap=seq, vmin=.2, vmax=1, edgecolor=BG, linewidth=.4)
    plt.colorbar(sc, ax=ax, shrink=.85).set_label("Tumour purity")
    ax.set_title(f"{ttl} coloured by purity", loc="left", fontweight="bold"); ax.set_xlabel(x); ax.set_ylabel(y)
    save(fig, f"gallery70_embed_{key}_purity.png", f"Map — {ttl} by purity",
         "Same projection coloured by tumour purity (fraction).", "Sample map", tags=[key, "purity"])
    # age
    fig, ax = plt.subplots(figsize=(5, 4.2))
    age = -emb.age_days / 365.25
    sc = ax.scatter(emb[x], emb[y], s=18, c=age, cmap=seq, edgecolor=BG, linewidth=.4)
    plt.colorbar(sc, ax=ax, shrink=.85).set_label("Age (years)")
    ax.set_title(f"{ttl} coloured by age", loc="left", fontweight="bold"); ax.set_xlabel(x); ax.set_ylabel(y)
    save(fig, f"gallery70_embed_{key}_age.png", f"Map — {ttl} by age",
         "Same projection coloured by age at diagnosis.", "Sample map", tags=[key, "age"])

# ---------- F. thresholds at N=10 (4) ----------
thr = [(0.0, 0.01), (0.0, 0.05), (0.5, 0.01), (0.5, 0.05), (1.0, 0.01), (1.0, 0.05)]
labels = [f"{'none' if l == 0 else f'|log2FC|>{l}'} / {int(f * 100)}%" for l, f in thr]
for key, ttl, ylab, col in [("rep", "Threshold sensitivity — replicability at N=10", "Median Jaccard", "rep_median"),
                            ("precision", "Threshold sensitivity — precision at N=10", "Median precision", None),
                            ("recall", "Threshold sensitivity — recall at N=10", "Median recall", None),
                            ("n_deg", "Threshold sensitivity — genes called at N=10", "Median genes", None)]:
    fig, ax = plt.subplots(figsize=(7, 3.6))
    vals = []
    sub = mc[(mc.N == 10)]
    for l, f in thr:
        h = sub[(sub.lfc == l) & (sub.fdr == f)]
        if key == "rep":
            v = float(rp[(rp.N == 10) & (rp.lfc == l) & (rp.fdr == f)].iloc[0]["rep_median"])
        elif key == "precision":
            v = float(h.precision.median())
        elif key == "recall":
            v = float(h.recall.median())
        else:
            v = float(h.n_deg.median())
        vals.append(v)
    ax.bar(range(len(thr)), vals, color=BLUE, alpha=.8)
    ax.set_xticks(range(len(thr))); ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=8)
    ax.set_title(ttl, loc="left", fontweight="bold"); ax.set_ylabel(ylab)
    save(fig, f"gallery80_thresh_{key}.png", ttl,
         "Fast engine, N=10, 100 cohorts; reference = own full-data call at same setting.", "Thresholds", n=10, tags=["threshold", key])

json.dump(manifest, open(os.path.join(OUT, "gallery.json"), "w"), indent=1)
print(f"wrote {len(manifest)} figures to {OUT}/")
