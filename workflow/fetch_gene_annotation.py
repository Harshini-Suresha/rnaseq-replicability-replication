"""Build data/gene_annotation.csv (Ensembl gene ID -> symbol, biotype, description) for the genes in the count matrix.

Source: the GRCh38 table of the R package 'annotables' (github.com/stephenturner/annotables, data/grch38.rda), read with pyreadr.
Network access to raw.githubusercontent.com is needed; the derived CSV is committed, so this only has to be rerun to refresh it.
Usage: python workflow/fetch_gene_annotation.py config/config.yaml
"""
import io, os, sys, tempfile, urllib.request
import pandas as pd, pyreadr, yaml

URL = "https://raw.githubusercontent.com/stephenturner/annotables/master/data/grch38.rda"
cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
genes = pd.read_csv(cfg["counts"], index_col=0, usecols=[0]).index
with tempfile.TemporaryDirectory() as d:
    f = os.path.join(d, "grch38.rda")
    urllib.request.urlretrieve(URL, f)
    a = list(pyreadr.read_r(f).values())[0]
a = a.drop_duplicates("ensgene").set_index("ensgene")
out = a.reindex(genes)[["symbol", "biotype", "description"]]
out["symbol"] = out["symbol"].replace("", pd.NA)
out.to_csv("data/gene_annotation.csv", index_label="ensembl_id")
print(f"{out.symbol.notna().sum()} of {len(out)} genes have a symbol; {out.index.isin(a.index).sum()} found in the table")
