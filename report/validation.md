# Engine validation

## 1. Full data (126 vs 126): fast engine vs the authors' DESeq2 (`BSLA.deseq2.lfc1.csv`)

Both use covariates (tumour purity, age) and a formal |log2FC| > 1 threshold test at 5 % FDR. logFC Spearman **0.965** over 18188 genes. DEGs: engine **1852**, authors **2039**, shared **1546** (Jaccard **0.659**).

## 2. The authors' N = 5 example cohort: fast engine vs their edgeR glmTreat result (`BSLA.N5.original.edgeR.lfc1.csv`)

logFC Spearman **0.966**. DEGs: engine **48**, edgeR **340**, shared **48**. The fast engine is much more conservative than edgeR at this cohort size: it calls few genes, almost all of which edgeR also calls.

## 3. Same cohorts, two engines: fast engine vs PyDESeq2 (covariates + formal |log2FC| > 1 test, 5 % FDR)

| N per group | cohorts | median DEGs, fast | median DEGs, PyDESeq2 | median share of fast calls also called by PyDESeq2 |
|---|---|---|---|---|
| 3 | 10 | 1 | 423 | 1.00 |
| 5 | 10 | 23 | 365 | 1.00 |
| 7 | 10 | 110 | 519 | 1.00 |
| 10 | 10 | 174 | 465 | 0.99 |
| 15 | 10 | 443 | 778 | 0.95 |

