## Table 1. Fast engine, 100 random cohorts per size (primary setting: |log2FC| > 1, 5 % FDR)

Precision and recall are against the authors' full-data DESeq2 result; medians over cohorts. Precision is undefined for cohorts that call no genes and is computed over the rest.

| N per group | Median genes called | Cohorts with no calls | Precision | Recall | Replicability (median Jaccard) | Log2FC rank agreement with full data |
|---|---|---|---|---|---|---|
| 3 | 0 | 64 % | 0.86 | 0.00 | 0.00 | 0.61 |
| 4 | 1 | 47 % | 1.00 | 0.00 | 0.00 | 0.72 |
| 5 | 22 | 21 % | 0.96 | 0.01 | 0.02 | 0.75 |
| 6 | 54 | 4 % | 0.96 | 0.03 | 0.13 | 0.80 |
| 7 | 96 | 1 % | 0.96 | 0.05 | 0.21 | 0.82 |
| 8 | 128 | 0 % | 0.95 | 0.06 | 0.24 | 0.84 |
| 9 | 175 | 0 % | 0.95 | 0.08 | 0.28 | 0.85 |
| 10 | 200 | 0 % | 0.95 | 0.09 | 0.32 | 0.86 |
| 12 | 322 | 0 % | 0.95 | 0.15 | 0.37 | 0.89 |
| 15 | 444 | 0 % | 0.95 | 0.21 | 0.44 | 0.91 |

## Table 2. PyDESeq2 anchor (covariates, formal |log2FC| > 1 test, 5 % FDR) vs the fast engine on the same cohorts

Medians over the cohorts run; replicability is the median Jaccard over all pairs of those cohorts. Reference for precision and recall: the authors' full-data DESeq2 result.

| N per group | Cohorts | PyDESeq2 genes called | PyDESeq2 precision | PyDESeq2 recall | PyDESeq2 replicability | Fast engine replicability, same cohorts | Fast engine precision, same cohorts |
|---|---|---|---|---|---|---|---|
| 3 | 10 | 423 | 0.71 | 0.14 | 0.14 | 0.00 | 0.82 |
| 5 | 10 | 365 | 0.84 | 0.15 | 0.20 | 0.02 | 0.92 |
| 7 | 10 | 519 | 0.86 | 0.22 | 0.29 | 0.25 | 0.97 |
| 10 | 10 | 465 | 0.91 | 0.21 | 0.33 | 0.33 | 0.95 |
| 15 | 10 | 778 | 0.91 | 0.35 | 0.41 | 0.43 | 0.95 |

## Table 3. Authors' published values for this data set (their results table, row BSLA)

| N per group | Precision | Recall | Replicability | Log2FC rank agreement with full data (mean) |
|---|---|---|---|---|
| 5 | 0.84 | 0.14 | 0.19 | 0.59 |
| 10 | 0.90 | 0.25 | 0.30 | 0.73 |

## Table 4. Sensitivity of the N = 10 result to the calling threshold (fast engine; reference = this repo's full-data call at the same setting)

| Fold-change threshold | FDR | Reference genes | Median genes called | Precision | Recall | Replicability |
|---|---|---|---|---|---|---|
| none | 1 % | 12146 | 1100 | 1.00 | 0.09 | 0.30 |
| none | 5 % | 13498 | 2648 | 0.99 | 0.19 | 0.35 |
| |log2FC| > 0.5 | 1 % | 4332 | 276 | 1.00 | 0.06 | 0.29 |
| |log2FC| > 0.5 | 5 % | 4928 | 600 | 0.99 | 0.12 | 0.32 |
| |log2FC| > 1.0 | 1 % | 1630 | 102 | 1.00 | 0.06 | 0.28 |
| |log2FC| > 1.0 | 5 % | 1852 | 200 | 1.00 | 0.11 | 0.32 |
