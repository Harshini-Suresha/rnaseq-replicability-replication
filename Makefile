# One-command rebuilds. `make all` skips the slow PyDESeq2 anchor (uses stored results/pydeseq2_anchor_sets.json).
CFG=config/config.yaml
.PHONY: all test sweep anchor validate tables figures embeddings dashboard serve

all: sweep embeddings validate tables figures dashboard

test:
	python tests/test_fast_de.py

sweep:        ## ~1 min: 1,000 random cohorts with the fast engine
	python workflow/run_sweep.py $(CFG) --workers 2

anchor:       ## slow (about 1 min per cohort); resumable
	python workflow/pydeseq2_crosscheck.py $(CFG) --sizes 3 5 7 10 15 --k 10 --cpus 2

validate:
	python workflow/validate_engine.py $(CFG)

tables:
	python workflow/make_report.py $(CFG)

annotation:   ## gene symbols (needs raw.githubusercontent.com)
	python workflow/fetch_gene_annotation.py

embeddings:   ## PCA, t-SNE, UMAP of all 252 samples (~20 s)
	python workflow/embeddings.py $(CFG)

figures:
	python workflow/make_figures.py $(CFG)
	python workflow/make_figures_extra.py $(CFG)

dashboard:
	python workflow/build_dashboard.py $(CFG)

serve:        ## open http://localhost:8000
	python -m http.server 8000 --directory docs
