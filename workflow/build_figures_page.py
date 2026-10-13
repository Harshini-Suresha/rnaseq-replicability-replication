"""Build docs/figures.html gallery from docs/figures/gallery.json (self-contained, light theme).
Usage: python workflow/build_figures_page.py
"""
import json, os

items = json.load(open("docs/figures/gallery.json"))
cats = sorted({d["category"] for d in items})
ns = sorted({d["n"] for d in items if d.get("n") is not None})
data_js = json.dumps(items, separators=(",", ":"))

html = """<!doctype html>
<html lang="en" data-theme="light">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Figure gallery — RNA-seq Replicability (100 figures)</title>
<meta name="description" content="Browsable gallery of 100 figures from the RNA-seq replicability reproduction, grouped by plot type and biology.">
<style>
:root{color-scheme:light;--surface-1:#fcfcfb;--surface-2:#f4f3f0;--text-primary:#0b0b0b;--text-secondary:#52514e;--text-muted:#6b6a66;--grid:#e8e7e3;--axis:#cfcec9;--s1:#2a78d6;--s2:#eb6834;--s3:#1baf7a}
:root[data-theme="dark"]{color-scheme:dark;--surface-1:#1a1a19;--surface-2:#232321;--text-primary:#f0efec;--text-secondary:#c3c2b7;--text-muted:#9b9a93;--grid:#2c2c2a;--axis:#3d3d3a;--s1:#3987e5;--s2:#d95926;--s3:#199e70}
*{box-sizing:border-box}html{-webkit-text-size-adjust:100%}
body{margin:0;background:var(--surface-1);color:var(--text-primary);font:15px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
main{max-width:1200px;margin:0 auto;padding:28px 16px 64px}
h1{font-size:26px;margin:0 0 6px;letter-spacing:-.01em}.lede{color:var(--text-secondary);max-width:78ch}
.topbar{display:flex;justify-content:space-between;align-items:flex-start;gap:12px}
button,select,input{font:inherit;color:var(--text-primary);background:var(--surface-1);border:1px solid var(--axis);border-radius:8px;padding:6px 10px;min-height:34px}
button{cursor:pointer}button:hover,select:hover{background:var(--surface-2)}
a{color:inherit}.back{display:inline-block;margin:12px 0;font-size:13.5px;color:var(--text-secondary)}
.toolbar{position:sticky;top:0;z-index:10;background:var(--surface-1);padding:12px 0;border-top:1px solid var(--grid);border-bottom:1px solid var(--grid);display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end;margin:16px 0}
.toolbar label{display:flex;flex-direction:column;gap:3px;font-size:12px;color:var(--text-muted)}
.toolbar input[type=search]{min-width:240px}
.count{font-size:13px;color:var(--text-secondary);margin:4px 0 12px}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(270px,1fr));gap:14px}
section.gsec{margin:22px 0}.gsec h2{font-size:17px;margin:0}.gsec .desc{font-size:13px;color:var(--text-secondary);margin:2px 0 10px;max-width:90ch}
#secnav a{border:1px solid var(--axis);border-radius:20px;padding:3px 10px;text-decoration:none;background:var(--surface-1)}
.card{background:var(--surface-1);border:1px solid var(--grid);border-radius:12px;overflow:hidden;display:flex;flex-direction:column}
.card img{width:100%;height:210px;object-fit:contain;object-position:center;padding:8px;background:var(--surface-2);cursor:zoom-in;display:block}
.card .body{padding:10px 12px 12px;display:flex;flex-direction:column;gap:6px;flex:1}
.badge{display:inline-block;font-size:11px;color:var(--text-secondary);background:var(--surface-2);border-radius:20px;padding:2px 9px;width:fit-content}
.card h3{font-size:14px;margin:0;line-height:1.35}.card p{font-size:12.5px;color:var(--text-secondary);margin:0;flex:1}
.row{display:flex;gap:8px}.row button,.row a{font-size:12.5px;min-height:30px;padding:4px 10px;text-decoration:none;border:1px solid var(--axis);border-radius:8px;background:var(--surface-1)}
.empty{padding:30px;text-align:center;color:var(--text-secondary)}
.overlay{position:fixed;inset:0;background:rgba(0,0,0,.72);display:none;align-items:center;justify-content:center;z-index:50;padding:20px}
.overlay.open{display:flex}.viewer{background:var(--surface-1);border-radius:12px;max-width:1020px;width:100%;max-height:92vh;overflow:auto;padding:16px}
.viewer img{width:100%;height:auto;border-radius:8px;background:var(--surface-2)}
.viewer h3{margin:10px 0 4px}.viewer p{color:var(--text-secondary);font-size:13.5px;margin:0 0 10px}
.nav{display:flex;gap:8px;align-items:center;justify-content:space-between;margin-top:6px}
section.prose{margin-top:34px;max-width:82ch}section.prose h2{font-size:18px;margin:0 0 6px}section.prose p{color:var(--text-secondary);font-size:14px}section.prose ul{padding-left:20px;margin:6px 0 12px;color:var(--text-secondary);font-size:14px}section.prose li{margin-bottom:5px}
code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.88em;background:var(--surface-2);padding:1px 5px;border-radius:4px}
pre{background:var(--surface-2);padding:12px 14px;border-radius:10px;overflow:auto;font-size:12.5px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:10px;margin:14px 0}
.stat{background:var(--surface-2);border-radius:10px;padding:10px 12px}.stat .v{font-size:20px;font-weight:650}.stat .k{font-size:12px;color:var(--text-secondary)}
details.data{font-size:12px;color:var(--text-secondary)}details.data summary{cursor:pointer;font-size:12px}
details.data table{width:100%;border-collapse:collapse;font-size:11.5px;margin-top:6px}details.data td{border-bottom:1px solid var(--grid);padding:3px 4px;vertical-align:top}details.data td:first-child{color:var(--text-muted);white-space:nowrap;width:70px}
footer{margin-top:40px;color:var(--text-muted);font-size:12.5px;border-top:1px solid var(--grid);padding-top:14px}
</style>
</head>
<body>
<main>
<div class="topbar"><div>
<h1>Figure gallery — 100 figures</h1>
<p class="lede">Every static figure from the reproduction, grouped two ways: by <strong>biology</strong> (overview, recurrence, distributions, single genes, sample map, thresholds) and by <strong>plot type</strong> (volcano, MA, heatmap, PCA / t-SNE / UMAP scatter, box + dots, violin, dot, bar, line, diagnostics, curves, histograms). Light background; use Theme for dark. Click any figure to enlarge.</p>
<a class="back" href="index.html">&larr; Back to interactive dashboard</a>
</div><button id="theme" aria-label="Toggle colour theme">Theme</button></div>
<div class="toolbar" role="group" aria-label="Gallery filters">
<label>Search<input type="search" id="q" placeholder="e.g. AGR3, UMAP, N=10, precision"></label>
<label>Category (biology)<select id="cat"><option value="">All categories</option></select></label>
<label>Plot type<select id="pt"><option value="">All plot types</option></select></label>
<label>Samples per group<select id="n"><option value="">All N</option></select></label>
<label>Sort<select id="sort"><option value="name">File name</option><option value="cat">Category</option><option value="plot">Plot type</option></select></label>
</div>
<div class="count" id="count" aria-live="polite"></div>
<div class="stats" aria-label="Dataset summary">
<div class="stat"><div class="v">252</div><div class="k">tumours (126 Basal + 126 LumA, TCGA-BRCA, BSLA)</div></div>
<div class="stat"><div class="v">18,188</div><div class="k">genes tested</div></div>
<div class="stat"><div class="v">1,000</div><div class="k">random cohorts (10 sizes × 100 trials)</div></div>
<div class="stat"><div class="v">100</div><div class="k">figures: 4 volcano, 4 MA, 4 heatmap, 44 box, 4 violin, 10 dot, 10 scatter, 7 bar, 1 line, 6 diagnostics, 2 curves, 4 histograms</div></div>
</div>
<nav id="secnav" aria-label="Sections" style="display:flex;gap:8px;flex-wrap:wrap;margin:4px 0 6px;font-size:12.5px"></nav>
<div id="sections"></div>
<div class="empty" id="empty" hidden>No figures match. Clear search or filters.</div>
<section class="prose" aria-labelledby="mqt"><h2 id="mqt">Methodology — how these 100 figures were made (in detail)</h2>
<p><strong>1. Question.</strong> When a differential-expression study has few replicates, how much would a repeat study of the same size find again (replicability), and how much of the full-data truth does one small study recover (precision/recall)? Design follows Degen &amp; Medo (2025): repeatedly draw small random cohorts from a large RNA-seq data set, run DE, and compare cohort-vs-cohort and cohort-vs-full-data.</p>
<p><strong>2. Data provenance.</strong> TCGA breast cancer, Basal vs Luminal A, 126 + 126 tumours, 18,188 genes after <code>min_total_count: 10</code> filtering — the authors' <code>BSLA</code> set, taken from <code>pdegen/BootstrapSeq</code> (MIT), not re-downloaded. Files: <code>data/BSLA.csv</code> (genes × samples counts), <code>data/BSLA.meta.csv</code> (<code>Condition, TumorPurity, days_to_birth</code>), <code>data/upstream_BSLA.deseq2.lfc1.csv</code> (authors' full-data DESeq2: <code>logCPM, logFC, lfcSE, pvalue, FDR</code>), <code>data/published_degen_medo_2025_results.csv</code> (their 18-dataset table, row BSLA), <code>data/gene_annotation.csv</code> (16,879/18,188 symbols). All parameters in <code>config/config.yaml</code> (<code>seed: 20261012</code>).</p>
<p><strong>3. Preprocessing (both engines).</strong> Size-factor normalisation → log2 CPM with prior 0.5 (as the DE engine) for embeddings; raw counts + covariates for modelling. Genes kept per-cohort if summed count ≥ 10.</p>
<p><strong>4. Cohort design (<code>workflow/cohorts.py</code>).</strong> N per group in {3, 4, 5, 6, 7, 8, 9, 10, 12, 15} × 100 trials = 1,000 cohorts, drawn without replacement within a cohort, seeded from <code>(seed, N, trial)</code> so reruns are identical. Same-N cohorts are independent draws and <em>can share samples</em>, as in the paper. Pairwise Jaccard uses all 4,950 pairs per N at the primary setting (<code>results/pair_jaccard.csv.gz</code>).</p>
<p><strong>5. DE models.</strong> (a) <em>Fast numpy engine</em> (<code>workflow/fast_de.py</code>, limma-trend style): weighted least squares on log2 CPM with purity + age covariates, empirical-Bayes moderated t, formal threshold test H0: |log2FC| ≤ L against L ∈ {0, 0.5, 1.0}, BH-FDR ∈ {1%, 5%}. Runs 1,000 cohorts in ~1 min but is conservative below ~8/group (64%/47%/21% no-calls at N = 3/4/5). (b) <em>PyDESeq2 anchor</em> (<code>workflow/pydeseq2_crosscheck.py</code>): Python port of DESeq2, Wald + formal |log2FC|&gt;1 test with covariates, ~1 min/cohort → 10 cohorts/size at N = 3, 5, 7, 10, 15 only (45 pairs/size; <code>results/pydeseq2_anchor_sets.json</code>). R/DESeq2/edgeR were unavailable at build time.</p>
<p><strong>6. References (“truth”).</strong> Two full-252 analyses: authors' DESeq2 (2,039 DEGs at |log2FC|&gt;1, FDR 5%) and own fast-engine full-data call (1,852 DEGs). Primary precision/recall in the gallery use the authors' reference.</p>
<ul>
<li><strong>Metrics (<code>workflow/run_sweep.py</code>).</strong> Per cohort: <code>precision = |DEG_cohort ∩ DEG_ref| / |DEG_cohort|</code> (NaN/undefined if no calls, excluded from medians); <code>recall = |∩| / |DEG_ref|</code>; <code>replicability = |A∩B| / |A∪B|</code> (Jaccard; 0 if either list empty), median over all pairs; <code>spearman_truth</code> = Spearman(logFC_cohort, logFC_full). Summaries: median + IQR (bands/boxes), whiskers 1.5×IQR, violins via Gaussian KDE.</li>
<li><strong>Threshold grid.</strong> 6 settings (LFC × FDR) in <code>results/metrics_by_cohort.csv</code> (6,000 rows); primary L=1.0, P=0.05. Volcano/MA sections show 4 of the 6.</li>
<li><strong>Embeddings (<code>workflow/embeddings.py</code>, seed 7).</strong> 2,000 most variable genes → z-score → PCA-50 → t-SNE (perplexity 30, PCA init) + UMAP (n_neighbors 15, min_dist 0.3) on the 50 PCs. Stats: PC1 23.0%, PC1–10 49.8%, silhouettes PCA2 0.587 / t-SNE 0.794 / UMAP 0.803; PC1–purity r = 0.26, PC1–age r = −0.17. Top-40 box genes = full-data DEGs with largest |log2FC|/SE.</li>
<li><strong>100 figures by section.</strong> Volcano (4) from authors' <code>logFC/pvalue/FDR</code>; MA (4) from <code>logCPM/logFC</code>; Heatmap (4) from <code>top_gene_expression.csv</code> (log2 CPM, raw + row z-score); PCA/t-SNE/UMAP scatter (10 = 1 overview + 9 projection×colouring) from <code>embeddings.csv</code>; Box+dots (44 = 4 distributions + 40 genes); Violin (4 distributions); Dot/recurrence (10, one per N) from <code>gene_recurrence.csv.gz</code>; Bar (7 = 1 published + 4 thresholds + 2 rankings); Line (1 metrics trend); QC/diagnostics (6: 2 p-hist + scree + 2 PC-covariate + library size); Curves (2 PR at N = 5/10); Histograms (4 Jaccard at N = 3/5/10/15).</li>
<li><strong>Units.</strong> PCA: standardised log2 CPM units. t-SNE/UMAP: arbitrary units (neighbourhoods only — gaps/sizes not quantitative). Expression: log2 CPM (zeros at panel bottom). Volcano y = −log10 p (capped at 300); MA x = logCPM (log axis).</li>
</ul>
<pre>pip install -r requirements.txt
python workflow/run_sweep.py config/config.yaml                  # ~1 min: 1,000 cohorts
python workflow/embeddings.py config/config.yaml                 # PCA / t-SNE / UMAP + top genes
python workflow/validate_engine.py config/config.yaml            # report/validation.json
python workflow/make_gallery_figures.py config/config.yaml       # 70 PNGs
python workflow/make_gallery_bioplots.py                          # +8 volcano/heatmap → 82
python workflow/make_gallery_more.py                              # +18 MA/QC/curves → 100
python workflow/build_figures_page.py                            # this page → docs/figures.html</pre>
<p><strong>Validation (see dashboard).</strong> Full-data fast vs authors' DESeq2: Spearman logFC 0.965, 1,546/2,039 shared (Jaccard 0.659). Authors' N = 5 edgeR cohort: Spearman 0.966; engine 48 vs edgeR 340 DEGs (48 shared — conservative at small N). Fast-vs-PyDESeq2 medians agree within 0.04 at N = 7/10/15. <strong>Limitations.</strong> One of 18 data sets; no enrichment analysis; GEO sets outside network allowlist; reference = same-sample full-data (precision flattering at large N); anchor = 10 cohorts/size (noisy medians); cohort-size list + BSLA label taken from authors' repo, not paper Table 1.</p>
</section>
<footer>Data: <a href="https://github.com/pdegen/BootstrapSeq" target="_blank" rel="noopener">pdegen/BootstrapSeq</a> (MIT), derived from TCGA. Paper: Degen &amp; Medo, PLOS Comput Biol 2025. Independent reproduction.</footer>
</main>
<div class="overlay" id="ov" role="dialog" aria-modal="true"><div class="viewer">
<img id="vimg" alt=""><h3 id="vtitle"></h3><p id="vcap"></p>
<div class="nav"><div style="display:flex;gap:8px"><button id="prev">&larr; Prev</button><button id="next">Next &rarr;</button></div>
<div style="display:flex;gap:8px"><a id="dl" download style="font-size:13px;border:1px solid var(--axis);border-radius:8px;padding:6px 10px;text-decoration:none">Download PNG</a><button id="close">Close</button></div></div>
</div></div>
<script>
const G = __DATA__;
const secs=document.querySelector("#sections"),secnav=document.querySelector("#secnav"),q=document.querySelector("#q"),cat=document.querySelector("#cat"),pt=document.querySelector("#pt"),nsel=document.querySelector("#n"),sort=document.querySelector("#sort"),count=document.querySelector("#count"),empty=document.querySelector("#empty");
const ORDER=[["volcano","Volcano","DEG calling in one view: log2FC vs significance. All volcanoes live here."],["ma","MA","Mean-expression vs effect. All MA plots live here."],["heatmap","Heatmap","Top-gene expression across 252 tumours. All heatmaps live here."],["diagnostic","QC / diagnostics","p-values, PCA scree, PC-covariate checks, library sizes."],["curve","Precision-recall curves","Per-cohort tradeoff at N=5 and N=10."],["histogram","Histograms","Full Jaccard distributions per N."],["scatter","PCA / t-SNE / UMAP (scatter)","Embeddings of all tumours; PCA scores in log2-CPM units, t-SNE/UMAP arbitrary. All in one place."],["box","Box + dots","Distributions and single-gene expression with every point drawn."],["violin","Violin","Same distributions as violins: shape + median (white) + IQR (black)."],["dot","Dot / recurrence","Gene-level detection rate vs effect size, one panel per N."],["bar","Bar","Published comparison + threshold sensitivity + rankings."],["line","Line","Metric trend by cohort size."]];
[...new Set(G.map(d=>d.category))].sort().forEach(c=>{const o=document.createElement("option");o.value=c;o.textContent=c+" ("+G.filter(d=>d.category===c).length+")";cat.appendChild(o);});
[...new Set(G.map(d=>d.plot_type))].sort().forEach(k=>{const lab=(G.find(d=>d.plot_type===k)||{}).plot_type_label||k;const o=document.createElement("option");o.value=k;o.textContent=lab+" ("+G.filter(d=>d.plot_type===k).length+")";pt.appendChild(o);});
[...new Set(G.filter(d=>d.n!=null).map(d=>d.n))].sort((a,b)=>a-b).forEach(n=>{const o=document.createElement("option");o.value=n;o.textContent="N = "+n;nsel.appendChild(o);});
ORDER.forEach(([k,t])=>{const a=document.createElement("a");a.href="#sec-"+k;a.textContent=t;secnav.appendChild(a);});
let cur=[];
function card(d,i){
  const card=document.createElement("div");card.className="card";
  const img=document.createElement("img");img.loading="lazy";img.src=d.file;img.alt=d.title;
  img.addEventListener("click",()=>openV(i));
  const body=document.createElement("div");body.className="body";
  const b=document.createElement("span");b.className="badge";b.textContent=d.category+" · "+(d.plot_type_label||d.plot_type)+(d.n!=null?" · N="+d.n:"");
  const h=document.createElement("h3");h.textContent=d.title;
  const p=document.createElement("p");p.textContent=d.caption;
  const det=document.createElement("details");det.className="data";
  const sum=document.createElement("summary");sum.textContent="Data";det.appendChild(sum);
  const tbl=document.createElement("table");const rows=[["File",d.file],["Category",d.category],["Plot type",(d.plot_type_label||d.plot_type)+" ("+d.plot_type+")"],["N",d.n!=null?String(d.n):"—"],["Tags",(d.tags||[]).join(", ")||"—"],["Source","results/ + data/ (see Methodology)"]];
  rows.forEach(([k,v])=>{const tr=document.createElement("tr");const a=document.createElement("td");a.textContent=k;const c=document.createElement("td");c.textContent=v;tr.appendChild(a);tr.appendChild(c);tbl.appendChild(tr);});
  det.appendChild(tbl);
  const row=document.createElement("div");row.className="row";
  const v=document.createElement("button");v.textContent="View";v.addEventListener("click",()=>openV(i));
  const a=document.createElement("a");a.textContent="Download";a.href=d.file;a.download="";
  row.appendChild(v);row.appendChild(a);
  body.appendChild(b);body.appendChild(h);body.appendChild(p);body.appendChild(det);body.appendChild(row);
  card.appendChild(img);card.appendChild(body);return card;
}
function render(){
 const s=q.value.trim().toLowerCase(),c=cat.value,t=pt.value,n=nsel.value;
 cur=G.filter(d=>{
  if(c&&d.category!==c)return false;
  if(t&&d.plot_type!==t)return false;
  if(n&&String(d.n)!==String(n))return false;
  if(s&&!(d.title+" "+d.caption+" "+d.file+" "+(d.tags||[]).join(" ")+" "+(d.plot_type_label||"")).toLowerCase().includes(s))return false;
  return true;});
 if(sort.value==="cat")cur.sort((a,b)=>a.category.localeCompare(b.category)||a.file.localeCompare(b.file));
 else if(sort.value==="plot")cur.sort((a,b)=>(a.plot_type||"").localeCompare(b.plot_type||"")||a.file.localeCompare(b.file));
 else cur.sort((a,b)=>a.file.localeCompare(b.file));
 secs.replaceChildren();
 ORDER.forEach(([k,titl,desc])=>{
  const list=cur.filter(d=>d.plot_type===k);
  if(!list.length)return;
  const sec=document.createElement("section");sec.className="gsec";sec.id="sec-"+k;
  const h=document.createElement("h2");h.textContent=titl+" ("+list.length+")";
  const p=document.createElement("p");p.className="desc";p.textContent=desc;
  const gr=document.createElement("div");gr.className="grid";
  list.forEach(d=>gr.appendChild(card(d,cur.indexOf(d))));
  sec.appendChild(h);sec.appendChild(p);sec.appendChild(gr);secs.appendChild(sec);});
 count.textContent="Showing "+cur.length+" of "+G.length+" figures in "+secs.children.length+" sections";
 empty.hidden=cur.length>0;
}
let vi=0;const ov=document.querySelector("#ov"),vimg=document.querySelector("#vimg"),vtitle=document.querySelector("#vtitle"),vcap=document.querySelector("#vcap"),dl=document.querySelector("#dl");
function openV(i){vi=i;const d=cur[i];if(!d)return;vimg.src=d.file;vimg.alt=d.title;vtitle.textContent=d.title+" — "+d.category+(d.n!=null?" · N="+d.n:"");vcap.textContent=d.caption+"  [ "+d.file+" | tags: "+(d.tags||[]).join(", ")+" ]";dl.href=d.file;ov.classList.add("open");}
function step(k){if(!cur.length)return;vi=(vi+k+cur.length)%cur.length;openV(vi);}
document.querySelector("#prev").addEventListener("click",()=>step(-1));
document.querySelector("#next").addEventListener("click",()=>step(1));
document.querySelector("#close").addEventListener("click",()=>ov.classList.remove("open"));
ov.addEventListener("click",e=>{if(e.target===ov)ov.classList.remove("open");});
document.addEventListener("keydown",e=>{if(!ov.classList.contains("open"))return;if(e.key==="Escape")ov.classList.remove("open");if(e.key==="ArrowRight")step(1);if(e.key==="ArrowLeft")step(-1);});
[q,cat,pt,nsel,sort].forEach(el=>el.addEventListener("input",render));
document.querySelector("#theme").addEventListener("click",()=>{const r=document.documentElement;r.setAttribute("data-theme",r.getAttribute("data-theme")==="dark"?"light":"dark");});
render();
</script>
</body>
</html>
"""
html = html.replace("__DATA__", data_js)
open("docs/figures.html", "w").write(html)
print(f"docs/figures.html written ({len(html)/1024:.0f} KB) for {len(items)} figures")
