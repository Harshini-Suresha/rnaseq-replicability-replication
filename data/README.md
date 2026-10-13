# Data provenance

Everything in this folder was copied unchanged from the authors' repository
[pdegen/BootstrapSeq](https://github.com/pdegen/BootstrapSeq) (MIT licence, copy in `UPSTREAM_LICENSE_BootstrapSeq`).
The underlying expression data derive from TCGA breast cancer (BRCA), accessed via the Genomic Data Commons.

| File here | Original name in `resources/` | Content |
|---|---|---|
| `BSLA.csv` | `BSLA.csv` | Raw counts, 18,188 genes x 252 samples (126 Basal, 126 LumA) |
| `BSLA.meta.csv` | `BSLA.meta.csv` | Sample metadata: tumour purity, age (`days_to_birth`), condition |
| `upstream_BSLA.deseq2.lfc1.csv` | `BSLA.deseq2.lfc1.csv` | Authors' full-data DESeq2 result (formal threshold log2FC = 1) |
| `upstream_BSLA.N5.csv`, `upstream_BSLA.N5.meta.csv` | `BSLA.N5.csv`, `BSLA.N5.meta.csv` | The authors' example 5-vs-5 cohort |
| `upstream_BSLA.N5.edgeR.lfc1.csv` | `BSLA.N5.original.edgeR.lfc1.csv` | Authors' edgeR glmTreat result on that cohort |
| `published_degen_medo_2025_results.csv` | `degen_medo_results.csv` | Authors' per-data-set results table (N = 5 and N = 10 only) |

Data-set code `BSLA` is the label used in the authors' repository. The paper's main text refers to a breast cancer
Luminal A vs Basal comparison under a different abbreviation; the mapping was not confirmed against the paper's Table 1.

## Gene annotation (added in this repository)

`gene_annotation.csv` maps the Ensembl gene IDs in the count matrix to HGNC-style symbols, biotype and description (16,879 of
18,188 genes have a symbol; the rest are shown by Ensembl ID). It is derived from the GRCh38 table of the R package
[annotables](https://github.com/stephenturner/annotables) (Stephen Turner, GPL-3; underlying annotation from Ensembl) by
`workflow/fetch_gene_annotation.py`. It carries that source's licence and attribution, separate from the MIT licence of the code.
Symbols reflect the annotation release used by that package and can differ from the release used to build the TCGA matrix.
