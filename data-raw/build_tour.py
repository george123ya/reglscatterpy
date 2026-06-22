"""Builds notebooks/reglscatterpy_tour.ipynb — the full feature tour.

Run:  python data-raw/build_tour.py   (then nbconvert --execute)
Kept in-repo so the tour can be regenerated when the API changes.
"""
import nbformat as nbf

nb = nbf.v4.new_notebook()
C = []
def md(s): C.append(nbf.v4.new_markdown_cell(s))
def code(s): C.append(nbf.v4.new_code_cell(s))

md("# reglscatterpy — feature tour\n\n"
   "Interactive WebGL scatterplots for single-cell data, with a **scanpy-style** API.\n\n"
   "Plots render **static by default** (a self-contained snapshot that reopens with no kernel, "
   "like a plotly figure). Pass **`interactive=True`** for the live, kernel-linked widget that "
   "round-trips a lasso selection back to Python (`w.selection`).\n\n"
   "Interactions: drag to pan, **Ctrl/Cmd + scroll to zoom** (plain scroll scrolls the page), "
   "**double-click to reset** the view.")

md("## Setup — a synthetic single-cell dataset\n"
   "Clusters with marker genes, a log-normalised `.raw` and a scaled `.X` (so the `use_raw` "
   "demo is meaningful), plus UMAP/PCA embeddings and numeric QC columns.")
code(
"import numpy as np, pandas as pd, anndata as ad\n"
"import reglscatterpy as rs\n"
"\n"
"rng = np.random.default_rng(0)\n"
"n, g, k = 3000, 50, 6\n"
"labels = rng.integers(0, k, n)\n"
"centers = rng.normal(0, 6, (k, 2))\n"
"umap = centers[labels] + rng.normal(0, 1.2, (n, 2))\n"
"counts = rng.poisson(0.5, (n, g)).astype(float)\n"
"for c in range(k):\n"
"    m = slice(c*6, c*6+3)            # 3 marker genes per cluster\n"
"    counts[labels == c, m] += rng.poisson(10, (int((labels == c).sum()), 3))\n"
"lib = counts.sum(1, keepdims=True); lib[lib == 0] = 1\n"
"lognorm = np.log1p(counts / lib * 1e4)             # .raw  (0..~5)\n"
"X = (lognorm - lognorm.mean(0)) / (lognorm.std(0) + 1e-9)   # .X (z-scored)\n"
"var = pd.DataFrame(index=[f'Gene{i}' for i in range(g)])\n"
"obs = pd.DataFrame({\n"
"    'cluster': pd.Categorical([f'C{l}' for l in labels]),\n"
"    'n_counts': counts.sum(1), 'n_genes': (counts > 0).sum(1),\n"
"    'pct_mito': rng.beta(1.5, 20, n) * 100, 'score': rng.random(n),\n"
"}, index=[f'cell{i}' for i in range(n)])\n"
"adata = ad.AnnData(X=X, obs=obs, var=var)\n"
"adata.raw = ad.AnnData(X=lognorm, var=var)\n"
"adata.obsm['X_umap'] = umap\n"
"adata.obsm['X_pca'] = rng.normal(0, 1, (n, 10))\n"
"adata")

md("## 1. Basics\n"
   "`basis=` picks the embedding (short names like `'umap'`/`'pca'` resolve to the `obsm` key). "
   "`color=` takes an `obs` column or a gene.")
code("rs.scatterplot(adata, basis='umap', color='cluster')   # an obs column")
code("rs.scatterplot(adata, basis='umap', color='Gene0')     # a gene")

md("## 2. Colour scales\n"
   "`cmap` (continuous), `palette` (categorical), `vmin`/`vmax` (numbers or `'p1'`/`'p99'` "
   "percentiles), `center_zero`.")
code("rs.scatterplot(adata, basis='umap', color='Gene0', cmap='magma', vmax='p99')")
code("rs.scatterplot(adata, basis='umap', color='cluster', palette='Dark2')")

md("### `use_raw` — matching scanpy's colour scale\n"
   "Like `sc.pl.umap`, a gene defaults to `.raw` (log-normalised, e.g. 0–5). Pass "
   "`use_raw=False` to colour by the scaled `.X` (z-scored, negative). Watch the colour bar.")
code("rs.scatterplot(adata, basis='umap', color='Gene0')                 # .raw (default)")
code("rs.scatterplot(adata, basis='umap', color='Gene0', use_raw=False)  # scaled .X")

md("## 3. Size & opacity\n"
   "`size=` is a scalar **or** a column/gene name (per-point). `opacity_by=` likewise.")
code("rs.scatterplot(adata, basis='umap', color='cluster', size=8)")
code("rs.scatterplot(adata, basis='umap', color='cluster', size='n_counts')   # size by a column")

md("## 4. Draw order (z-depth)\n"
   "`sort_order=True` (default) draws higher continuous values on top. `random_state=N` shuffles "
   "the draw order (seeded) so no category is systematically hidden by overplotting.")
code("rs.scatterplot(adata, basis='umap', color='Gene0', sort_order=True)   # high values on top")
code("rs.scatterplot(adata, basis='umap', color='cluster', random_state=0)  # seeded shuffle")

md("## 5. `na_color` & `groups`\n"
   "`groups=[...]` keeps only those categories coloured; the rest grey out to `na_color`.")
code("rs.scatterplot(adata, basis='umap', color='cluster', groups=['C0', 'C3'], na_color='#dddddd')")

md("## 6. Other components / embeddings\n"
   "`components=(i, j)` is 1-based (scanpy). Plot PC2 vs PC3:")
code("rs.scatterplot(adata, basis='pca', color='cluster', components=(2, 3))")

md("## Big data — atlas-scale rendering\n"
   "`scatterplot()` stays interactive on huge data **without silently hiding cells**. "
   "Let's make a larger synthetic atlas to demo the three modes.")
code(
"bn, bk = 1_500_000, 12\n"
"blab = rng.integers(0, bk, bn)\n"
"bcent = rng.normal(0, 8, (bk, 2))\n"
"bumap = bcent[blab] + rng.normal(0, 1.0, (bn, 2))\n"
"big = ad.AnnData(\n"
"    X=np.zeros((bn, 1), dtype='float32'),\n"
"    obs=pd.DataFrame({'cell_type': pd.Categorical([f'T{l}' for l in blab]),\n"
"                      'score': rng.random(bn).astype('float32')}),\n"
")\n"
"big.obsm['X_umap'] = bumap.astype('float32')\n"
"big")

md("### Auto density subsample (default)\n"
   "`max_points='auto'` (the default) caps at 500k via a **density-preserving** sketch that thins "
   "dense blobs but keeps rare types. It's honest about it: an on-plot `'X of Y shown'` caption + a "
   "one-time warning. `subsample='random'` is the uniform fallback.")
code("rs.scatterplot(big, basis='umap', color='cell_type')   # caption + downsample warning")

md("### All points resident (ABC-Atlas style)\n"
   "`max_points=None` puts **every** cell on the GPU; pan/zoom is camera-only. Smooth up to ~4M "
   "points on a decent GPU.")
code("rs.scatterplot(big, basis='umap', color='cell_type', max_points=None)")

md("### Progressive detail-on-zoom (>4M)\n"
   "`progressive=True` shows a light density overview, then re-renders **all cells inside the "
   "viewport** as you zoom in (no preprocessing; always the live widget). Tune with "
   "`progressive_opts={'detail_max_points': ..., 'overscan': ...}` (points per viewport / margin fetched).")
code("rs.scatterplot(big, basis='umap', color='cell_type', progressive=True,\n"
     "               progressive_opts={'detail_max_points': 300_000, 'overscan': 0.6})")

md("## 7. Multi-panel grid\n"
   "A **list** of names → one linked panel per value (genes and/or obs), camera + lasso synced. "
   "`ncols` sets the columns. (A single-element list is just one plot.)")
code("rs.scatterplot(adata, basis='umap', color=['cluster', 'Gene0', 'Gene6'], ncols=3)")

md("## 8. Filtering\n"
   "`filter_by=` shows distribution sliders: drag the **black grabbers**, or drag the **band** to "
   "pan the range; filtering is **live** while you drag (auto-deferred above 150k points).")
code("rs.scatterplot(adata, basis='umap', color='cluster',\n"
     "               filter_by=adata.obs[['n_counts', 'n_genes', 'pct_mito']])")

md("## 9. Richer tooltips\n"
   "`tooltip_by=` adds fields on hover (obs columns or genes).")
code("rs.scatterplot(adata, basis='umap', color='cluster',\n"
     "               tooltip_by=['n_genes', 'pct_mito', 'Gene0'])")

md("## 10. Toolbar & view\n"
   "`toolbar='left'`/`'top'`/`'none'`: pan, lasso, zoom-to-selection, reset, screenshot. "
   "`zoom_on_selection=True` auto-frames a lasso. Double-click resets; Ctrl/Cmd+scroll zooms.")
code("rs.scatterplot(adata, basis='umap', color='cluster', toolbar='left', zoom_on_selection=True)")

md("## 11. Interactive selection round-trip 🔴 *(needs `interactive=True` + a live kernel)*\n"
   "This is the live widget. **Lasso** some cells in the plot (toolbar → lasso), then read them "
   "back in the next cell. You can also set the selection from Python.")
code("w = rs.scatterplot(adata, basis='umap', color='cluster', interactive=True, toolbar='left')\n"
     "w   # lasso a cluster, then run the cells below")
code("# the lassoed cells, as positional indices into adata (empty until you lasso):\n"
     "w.selection")
code("# drive it from Python too — this highlights those points in the plot above:\n"
     "w.selection = list(np.where(adata.obs['cluster'] == 'C0')[0])\n"
     "len(w.selection)")

md("### Analyse the selection\n"
   "`subset`, `diff_expression` (uses scanpy `rank_genes_groups` when installed, else a "
   "Welch t-test), `composition`, and `annotate` (writes back to `adata.obs`).")
code("sub = w.subset()          # adata[w.selection]\n"
     "sub")
code("w.diff_expression(n=5)    # top markers of the selection vs the rest")
code("w.composition('cluster')  # what clusters the selection is made of")
code("w.annotate('my_label', 'group A')   # writes adata.obs['my_label'] for the selected cells\n"
     "rs.scatterplot(adata, basis='umap', color='my_label')")

md("## 12. Linked grid (`compose`)\n"
   "Pass plots to `rs.compose(...)` — pan/zoom and lasso stay in sync across panels. "
   "`compose` makes the panels interactive for you (no need for `interactive=True` on each).")
code("a = rs.scatterplot(adata, basis='umap', color='cluster')\n"
     "b = rs.scatterplot(adata, basis='pca',  color='cluster')\n"
     "rs.compose([a, b])")

md("## 13. Width\n"
   "Plots are **700 px** by default. Pass `width=` (px), or `width=None` to fill the cell.")
code("rs.scatterplot(adata, basis='umap', color='cluster', width=400)")

md("## 14. Other inputs & export\n"
   "DataFrames (give `x`/`y` columns) and numpy arrays work too. `save='plot.html'` writes a "
   "self-contained offline file; `w.to_html(...)` does the same for one plot.")
code("df = pd.DataFrame({'x': umap[:, 0], 'y': umap[:, 1], 'cluster': obs['cluster'].values})\n"
     "rs.scatterplot(df, x='x', y='y', color='cluster')")
code("import tempfile, os\n"
     "out = os.path.join(tempfile.gettempdir(), 'umap.html')\n"
     "rs.scatterplot(adata, basis='umap', color='cluster', save=out)\n"
     "print('wrote', out)")

md("---\n"
   "**Static vs interactive:** the default static plots reopen with no kernel (great for sharing / "
   "reports); `interactive=True` adds the Python round-trip but, like any Jupyter widget, needs a "
   "running kernel to render. Use the default for figures you keep, `interactive=True` while "
   "actively selecting.")

nb.cells = C
nb.metadata["kernelspec"] = {"display_name": "reglpy", "language": "python", "name": "python3"}
nbf.write(nb, "notebooks/reglscatterpy_tour.ipynb")
print("wrote notebooks/reglscatterpy_tour.ipynb with", len(C), "cells")
