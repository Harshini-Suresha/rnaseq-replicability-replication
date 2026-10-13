"""Sample-level views of the full data set: PCA, t-SNE, UMAP, plus expression of top genes for box/violin plots.

Pipeline: size-factor normalisation -> log2 CPM (prior 0.5, as the DE engine) -> 2,000 most variable genes ->
z-score per gene -> PCA (50 components) -> t-SNE and UMAP on the 50 PCs. All seeded.
Outputs: results/embeddings.csv, results/pca_variance.csv, results/top_gene_expression.csv, results/embedding_stats.json
Usage: python workflow/embeddings.py config/config.yaml
"""
import json, os, sys
import numpy as np, pandas as pd, yaml
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score
import umap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fast_de

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
SEED = 7
counts = pd.read_csv(cfg["counts"], index_col=0).T
meta = pd.read_csv(cfg["meta"], index_col=0).loc[counts.index]
counts = counts.loc[:, counts.sum() >= cfg["min_total_count"]]
sf = fast_de.size_factors(counts)
norm = counts.to_numpy(float) / sf[:, None]
Y = np.log2((norm + 0.5) / (norm.sum(1).mean() + 1.0) * 1e6)
logcpm = pd.DataFrame(Y, index=counts.index, columns=counts.columns)

hv = logcpm.var(axis=0).nlargest(2000).index
Z = logcpm[hv]; Z = (Z - Z.mean()) / Z.std(ddof=0)
pca = PCA(n_components=50, random_state=SEED).fit(Z.to_numpy())
P = pca.transform(Z.to_numpy())
tsne = TSNE(n_components=2, perplexity=30, init="pca", learning_rate="auto", random_state=SEED).fit_transform(P)
um = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.3, random_state=SEED, n_jobs=1).fit_transform(P)

emb = pd.DataFrame({"sample": counts.index, "condition": meta[cfg["condition_col"]].values,
                    "purity": meta["TumorPurity"].values, "age_days": meta["days_to_birth"].values,
                    "PC1": P[:, 0], "PC2": P[:, 1], "PC3": P[:, 2], "tsne1": tsne[:, 0], "tsne2": tsne[:, 1],
                    "umap1": um[:, 0], "umap2": um[:, 1]})
emb.round(4).to_csv("results/embeddings.csv", index=False)
pd.DataFrame({"PC": np.arange(1, 11), "explained_variance_ratio": pca.explained_variance_ratio_[:10]}).round(5).to_csv("results/pca_variance.csv", index=False)

# genes for box / violin plots: the full-data reference DEGs with the largest standardised effect
ref = pd.read_csv("results/reference_full_data.csv.gz", index_col=0)
sel = ref[ref[f"deg_lfc{float(cfg['primary_lfc'])}_fdr{float(cfg['primary_fdr'])}"] == 1]
sel = (sel.log2FoldChange.abs() / sel.se).nlargest(40).index
tg = logcpm[sel].copy(); tg.insert(0, "condition", meta[cfg["condition_col"]].values)  # symbols are in data/gene_annotation.csv
tg.round(3).to_csv("results/top_gene_expression.csv", index_label="sample")

y = emb.condition.values
stats = dict(n_hvg=2000, pc1_var=round(float(pca.explained_variance_ratio_[0]), 4), pc2_var=round(float(pca.explained_variance_ratio_[1]), 4),
             pc1_to_10_var=round(float(pca.explained_variance_ratio_[:10].sum()), 4),
             silhouette_pca2=round(float(silhouette_score(P[:, :2], y)), 3), silhouette_tsne=round(float(silhouette_score(tsne, y)), 3),
             silhouette_umap=round(float(silhouette_score(um, y)), 3),
             pc1_purity_corr=round(float(np.corrcoef(P[:, 0], meta["TumorPurity"])[0, 1]), 3),
             pc1_age_corr=round(float(np.corrcoef(P[:, 0], meta["days_to_birth"])[0, 1]), 3))
json.dump(stats, open("results/embedding_stats.json", "w"), indent=1)
print(json.dumps(stats, indent=1))
