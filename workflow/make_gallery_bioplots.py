"""Meaningful volcano + heatmap figures for the gallery (non-generic).
Uses authors' full-data DESeq2 (real p-values) + top-gene expression matrix.
Writes 8 PNGs into docs/figures/ and appends to gallery.json.
Usage: python workflow/make_gallery_bioplots.py
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

au = pd.read_csv("data/upstream_BSLA.deseq2.lfc1.csv", index_col=0)
annot = pd.read_csv("data/gene_annotation.csv", index_col=0) if os.path.exists("data/gene_annotation.csv") else pd.DataFrame(columns=["symbol"])
def sym(g):
    s = annot.symbol.get(g) if g in annot.index else None
    return s if isinstance(s, str) and s else g

au["neglog10p"] = -np.log10(au["pvalue"].clip(lower=1e-300))
au["is_deg"] = (au["FDR"] < 0.05) & (au["logFC"].abs() > 1)
items = json.load(open("docs/figures/gallery.json"))

def volcano(fname, title, caption, lfc_thr, fdr_thr, label_n=6, xlim=10):
    d = au.copy()
    d["hit"] = (d["FDR"] < fdr_thr) & (d["logFC"].abs() > lfc_thr)
    fig, ax = plt.subplots(figsize=(6, 4.6))
    ax.scatter(d.loc[~d["hit"], "logFC"], d.loc[~d["hit"], "neglog10p"], s=6, color=GREY, alpha=.5, lw=0, label=f"not DEG ({(~d['hit']).sum():,})")
    ax.scatter(d.loc[d["hit"], "logFC"], d.loc[d["hit"], "neglog10p"], s=8, color=BLUE, alpha=.65, lw=0, label=f"DEG ({d['hit'].sum():,})")
    ax.axvline(lfc_thr, color="#8d8c86", ls="--", lw=1); ax.axvline(-lfc_thr, color="#8d8c86", ls="--", lw=1)
    if fdr_thr < 1:
        ax.axhline(-np.log10(fdr_thr), color="#8d8c86", ls=":", lw=1)
    top = d[d["hit"]].nlargest(label_n, "neglog10p") if d["hit"].sum() else d.nlargest(label_n, "neglog10p")
    for g, r in top.iterrows():
        ax.text(r["logFC"], r["neglog10p"], " " + sym(g), fontsize=7.5, color=INK, va="bottom")
    ax.set_xlim(-xlim, xlim)
    ax.set_xlabel("log2 fold change (LumA vs Basal)"); ax.set_ylabel("-log10 p-value (authors' full-data DESeq2)")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold"); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(f"docs/figures/{fname}", dpi=150); plt.close(fig)
    items.append(dict(file=f"figures/{fname}", title=title, caption=caption, category="Volcano",
                      n=None, tags=["volcano", "DEG", f"FDR{fdr_thr}", f"LFC{lfc_thr}"], plot_type="volcano", plot_type_label="Volcano"))

volcano("gallery40_volcano_primary.png", "Volcano — full data, primary (|log2FC|>1, FDR 5%)",
        "Authors' full-data DESeq2, 18,188 genes; blue = 2,039 DEGs. Labels = top significant DEGs.", 1.0, 0.05)
volcano("gallery40_volcano_nofilter.png", "Volcano — full data, no fold-change filter (FDR 5%)",
        "Same p-values, DEG = FDR<5% only (13,498 genes). Shows how the |log2FC|>1 filter shortens the list.", 0.0, 0.05)
volcano("gallery40_volcano_strict.png", "Volcano — full data, strict (FDR 1%, |log2FC|>1)",
        "DEG = FDR<1% + |log2FC|>1 (1,630 genes). Only the tip of the volcano survives.", 1.0, 0.01)
volcano("gallery40_volcano_half.png", "Volcano — full data, |log2FC|>0.5 (FDR 5%)",
        "DEG = FDR<5% + |log2FC|>0.5 (4,928 genes). Intermediate stringency.", 0.5, 0.05)

# ---------- heatmaps ----------
tg = pd.read_csv("results/top_gene_expression.csv", index_col=0)
genes40 = [c for c in tg.columns if c != "condition"][:40]
genes20 = genes40[:20]
cond = tg["condition"].values
order = np.argsort(cond)  # Basal then LumA (alphabetical) — groups together
cond_o = cond[order]

def heatmap(fname, title, caption, genes, zscore=True):
    mat = tg[genes].to_numpy()[order, :].T  # genes x samples
    labels = [sym(g) for g in genes]
    if zscore:
        mu = mat.mean(axis=1, keepdims=True); sd = mat.std(axis=1, keepdims=True) + 1e-9
        mat = (mat - mu) / sd
        vmin, vmax, clab = -3, 3, "z-score per gene (row)"
        cmap = plt.get_cmap("RdBu_r")
    else:
        vmin, vmax, clab = np.percentile(mat, 2), np.percentile(mat, 98), "log2 CPM"
        cmap = plt.get_cmap("viridis")
    n_g, n_s = mat.shape
    fig, ax = plt.subplots(figsize=(9, max(4, n_g * 0.22 + 1.5)))
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    # condition strip on top
    strip = np.array([0 if c == "Basal" else 1 for c in cond_o])[None, :]
    from matplotlib.colors import ListedColormap
    ax2 = ax.inset_axes([0, 1.005, 1, 0.04])
    ax2.imshow(strip, aspect="auto", cmap=ListedColormap([BLUE, ORANGE]), vmin=0, vmax=1)
    ax2.set_xticks([]); ax2.set_yticks([]); ax2.set_title("Basal (blue) → LumA (orange), n=252", fontsize=8, color=MUT, loc="left")
    ax.set_yticks(range(n_g)); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xticks([]); ax.set_xlabel(f"Samples ordered by subtype (252) — {clab}")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
    cb = plt.colorbar(im, ax=ax, shrink=.8); cb.set_label(clab)
    fig.tight_layout(); fig.savefig(f"docs/figures/{fname}", dpi=150); plt.close(fig)
    items.append(dict(file=f"figures/{fname}", title=title, caption=caption, category="Heatmap",
                      n=None, tags=["heatmap", "top-genes", "expression"], plot_type="heatmap", plot_type_label="Heatmap"))

heatmap("gallery50_heatmap_top20_z.png", "Heatmap — top 20 genes (z-scored)",
        "20 largest-effect genes × 252 tumours, row z-scores. Clean Basal/LumA split — why small cohorts recover these first.", genes20, True)
heatmap("gallery50_heatmap_top40_z.png", "Heatmap — top 40 genes (z-scored)",
        "All 40 gallery genes × 252 tumours, row z-scores. Effect fades down the list.", genes40, True)
heatmap("gallery50_heatmap_top20_raw.png", "Heatmap — top 20 genes (log2 CPM)",
        "Same 20 genes, raw log2 CPM (2nd–98th pct scale). High expressors dominate visually.", genes20, False)
heatmap("gallery50_heatmap_top40_raw.png", "Heatmap — top 40 genes (log2 CPM)",
        "Same 40 genes, raw log2 CPM. Compare with z-scored version above.", genes40, False)

json.dump(items, open("docs/figures/gallery.json", "w"), indent=1)
print("total now", len(items))
