"""Static PNG figures for the report (same palette as the dashboard). Usage: python workflow/make_figures.py config/config.yaml"""
import json, os, sys
import numpy as np, pandas as pd, yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cfg = yaml.safe_load(open(sys.argv[1] if len(sys.argv) > 1 else "config/config.yaml"))
L, P = float(cfg["primary_lfc"]), float(cfg["primary_fdr"])
C = dict(fast="#2a78d6", py="#eb6834", pub="#1baf7a", ink="#0b0b0b", mut="#52514e", grid="#e8e7e3", bg="#fcfcfb")
plt.rcParams.update({"figure.facecolor": C["bg"], "axes.facecolor": C["bg"], "axes.edgecolor": "#cfcec9", "axes.labelcolor": C["mut"],
                     "xtick.color": C["mut"], "ytick.color": C["mut"], "text.color": C["ink"], "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": C["grid"],
                     "grid.linewidth": 0.8, "axes.axisbelow": True, "font.family": "DejaVu Sans"})
os.makedirs("figures", exist_ok=True)
mc = pd.read_csv("results/metrics_by_cohort.csv"); rp = pd.read_csv("results/replicability_by_N.csv")
mc = mc[(mc.lfc == L) & (mc.fdr == P)]; rp = rp[(rp.lfc == L) & (rp.fdr == P)].set_index("N")
pub = pd.read_csv("data/published_degen_medo_2025_results.csv", index_col=0)
au = pd.read_csv(cfg["authors_reference"], index_col=0); auth = set(au.index[au.FDR < cfg["authors_reference_fdr"]])
sizes = cfg["cohort_sizes"]
anchor = {}
if os.path.exists("results/pydeseq2_anchor_sets.json"):
    st = json.load(open("results/pydeseq2_anchor_sets.json"))
    for n in sorted({int(k.split("_")[0]) for k in st}):
        ts = sorted(int(k.split("_")[1]) for k in st if int(k.split("_")[0]) == n)
        ss = [set(st[f"{n}_{t}"]["pydeseq2"]) for t in ts]
        rep = [len(ss[i] & ss[j]) / len(ss[i] | ss[j]) for i in range(len(ss)) for j in range(i + 1, len(ss)) if ss[i] | ss[j]]
        anchor[n] = dict(k=len(ts), rep=rep, prec=[len(s & auth) / len(s) for s in ss if s], rec=[len(s & auth) / len(auth) for s in ss])


def panel(ax, title, ylab, fastq, anch=None, pubv=None):
    x = np.array(sizes)
    ax.fill_between(x, [fastq[n][0] for n in sizes], [fastq[n][2] for n in sizes], color=C["fast"], alpha=.12, lw=0)
    ax.plot(x, [fastq[n][1] for n in sizes], color=C["fast"], lw=2, marker="o", ms=6, mec=C["bg"], mew=1.5, label="Fast engine (100 cohorts per size; band = IQR)")
    if anch:
        ns = sorted(anch)
        for n in ns:
            v = anch[n]
            if v:
                ax.plot([n, n], [min(v), max(v)], color=C["py"], lw=2, alpha=.5)
        ax.plot(ns, [np.median(anch[n]) for n in ns], color=C["py"], lw=2, marker="s", ms=6, mec=C["bg"], mew=1.5, label="PyDESeq2 (few cohorts per size; bar = range)")
    if pubv:
        ax.scatter(list(pubv), list(pubv.values()), marker="D", s=46, color=C["pub"], edgecolor=C["bg"], linewidth=1.5, zorder=5, label="Published (authors' table)")
    ax.set_title(title, loc="left", fontsize=11, fontweight="bold"); ax.set_xlabel("Samples per group (N)"); ax.set_ylabel(ylab)
    ax.set_ylim(0, 1.02); ax.set_xticks(sizes)


q = lambda s: [s.quantile(.25), s.median(), s.quantile(.75)]
fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
panel(axs[0], "Replicability", "Median Jaccard overlap of DEG lists",
      {n: [rp.loc[n, "rep_q25"], rp.loc[n, "rep_median"], rp.loc[n, "rep_q75"]] for n in sizes},
      {n: a["rep"] for n, a in anchor.items()}, {5: pub.loc["BSLA", "Rep_N5"], 10: pub.loc["BSLA", "Rep_N10"]})
panel(axs[1], "Precision", "Share of called genes also in the reference",
      {n: q(mc[mc.N == n].precision_authors.dropna()) for n in sizes}, {n: a["prec"] for n, a in anchor.items()},
      {5: pub.loc["BSLA", "Prec_N5"], 10: pub.loc["BSLA", "Prec_N10"]})
panel(axs[2], "Recall", "Share of reference genes that are called",
      {n: q(mc[mc.N == n].recall_authors) for n in sizes}, {n: a["rec"] for n, a in anchor.items()},
      {5: pub.loc["BSLA", "Rec_N5"], 10: pub.loc["BSLA", "Rec_N10"]})
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(.5, -.02))
fig.tight_layout(rect=(0, .06, 1, 1)); fig.savefig("figures/fig1_metrics_by_cohort_size.png", dpi=160); plt.close(fig)

# fig 2: gene recurrence at N = 5, 10, 15
rec = pd.read_csv("results/gene_recurrence.csv.gz", index_col=0); rd = rec[rec.full_data_deg == 1]
fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), sharey=True)
for ax, n in zip(axs, (5, 10, 15)):
    ax.scatter(rd.full_data_log2FC, rd[f"N{n}"] * 100, s=7, color=C["fast"], alpha=.4, lw=0)
    ax.set_title(f"N = {n} per group", loc="left", fontsize=11, fontweight="bold"); ax.set_xlabel("Full-data log2 fold change (LumA vs Basal)")
    ax.set_ylim(-2, 102); ax.axvline(0, color="#cfcec9", lw=1)
axs[0].set_ylabel("Cohorts that call the gene (%)")
fig.tight_layout(); fig.savefig("figures/fig2_gene_recurrence.png", dpi=160); plt.close(fig)

# fig 3: published replicability across the 18 data sets
fig, axs = plt.subplots(1, 2, figsize=(11, 5.2), sharey=False)
for ax, k, n in zip(axs, ("Rep_N5", "Rep_N10"), (5, 10)):
    d = pub[k].sort_values()
    cols = [C["pub"] if i == "BSLA" else "#a9dcc8" for i in d.index]
    ax.barh(d.index, d.values, color=cols, height=.62)
    mine = rp.loc[n, "rep_median"] if n in rp.index else None
    if mine is not None and "BSLA" in d.index:
        ax.scatter([mine], ["BSLA"], color=C["fast"], s=60, zorder=5, edgecolor=C["bg"], linewidth=1.5, label="This reproduction (fast engine)")
    if n in anchor and anchor[n]["rep"]:
        ax.scatter([np.median(anchor[n]["rep"])], ["BSLA"], marker="s", color=C["py"], s=60, zorder=6, edgecolor=C["bg"], linewidth=1.5, label="This reproduction (PyDESeq2)")
    ax.set_title(f"Replicability at N = {n}", loc="left", fontsize=11, fontweight="bold"); ax.set_xlim(0, 1); ax.grid(axis="y", visible=False)
axs[0].legend(frameon=False, loc="lower right", fontsize=8.5)
fig.tight_layout(); fig.savefig("figures/fig3_published_18_datasets.png", dpi=160); plt.close(fig)
print("figures written:", sorted(os.listdir("figures")))
