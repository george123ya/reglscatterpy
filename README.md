# reglscatterpy

Interactive WebGL scatterplots for single-cell / spatial data in Python —
**AnnData, MuData, SpatialData**, pandas, numpy. Renders millions of points in
the browser via [`regl-scatterplot`](https://github.com/flekschas/regl-scatterplot),
in **Jupyter, JupyterLab, VS Code and Colab**.

This is the Python companion to the R package
[**reglScatterplotR**](https://github.com/george123ya/reglScatterplotR). Both
drive the *same* compiled widget, so a plot looks and behaves identically across
R and Python — the draggable legend, `filter_by` distribution sliders, lasso,
tooltips and PNG/SVG/PDF export all come from one shared codebase. (Equivalence
is locked down by `tests/test_payload_parity.py`, which checks the Python
payload byte-for-byte against R fixtures.)

## Install

```bash
pip install reglscatterpy            # numpy, pandas, anywidget
pip install anndata                  # for AnnData; mudata / spatialdata as needed
```

## Quick start

```python
import scanpy as sc
import reglscatterpy as rs

adata = sc.datasets.pbmc3k_processed()
rs.scatterplot(adata, x="X_umap", color_by="louvain")   # an obs column
rs.scatterplot(adata, x="X_umap", color_by="CST3")      # a gene
```

```python
import numpy as np, pandas as pd
df = pd.DataFrame({"x": np.random.rand(10_000), "y": np.random.rand(10_000),
                   "ct": np.random.choice(list("ABC"), 10_000)})
rs.scatterplot(df, x="x", y="y", color_by="ct")
```

Plots fill the notebook cell width by default; pass `width=` (pixels) for a
fixed size.

> **Note:** like other Jupyter widgets, a plot's large state isn't reliably
> saved into the `.ipynb`, so after **reopening** a notebook the cell may show
> blank (or `Could not render … widget-view`) until you **re-run** it. To keep
> an interactive copy that survives reopening — and to share a plot with someone
> who has no kernel — export it to a standalone HTML file (see below).

## Save a standalone HTML (offline, kernel-free)

The Python equivalent of R's `htmlwidgets::saveWidget`: write a single
self-contained `.html` that **inlines the widget and the plot's data**, so it
opens in any browser with no kernel and no internet:

```python
w = rs.scatterplot(adata, x="X_umap", color_by="leiden")
rs.save_html(w, "umap.html")      # or:  w.to_html("umap.html")
```

The saved file is fully interactive (pan/zoom, legend, lasso, tooltips,
PNG/SVG/PDF export) but it's a **snapshot** — it has no kernel, so the Python
round-trips (`w.selection`, `w.annotate`, …) only work in the live notebook. The
widget bundle is inlined gzip-compressed (~0.5 MB, decompressed in-browser), so
a one-plot file is well under 1 MB. No R is involved — it's pure Python.

### A whole notebook → one HTML report

Plain `jupyter nbconvert --to html` leaves the plots blank (the same widget-state
limitation). `save_notebook_html` re-executes the notebook and bakes every plot
in as an interactive, kernel-free figure, sharing **one** copy of the bundle
across all plots:

```python
rs.save_notebook_html("analysis.ipynb", "analysis_report.html")
```

Needs `nbconvert` + `ipykernel` (`pip install reglscatterpy[report]`). The plots
are fully offline; note that nbconvert's own page chrome (MathJax/RequireJS) is
still referenced from a CDN — use [`nb_offline_convert`](https://github.com/trungleduc/nb_offline_convert)
if you need the surrounding report shell to be 100% offline too.

## Selection round-trip

Lasso points in the plot, then read them back in another cell — or drive the
selection from Python:

```python
w = rs.scatterplot(adata, x="X_umap", color_by="leiden")
w                          # show it, lasso some cells in the widget

w.selection                # -> [12, 87, 134, ...]  positional indices
adata[w.selection]         # subset the AnnData directly
sub = w.subset()           # same thing, as a convenience

w.selection = list(range(100))   # or set it from Python to highlight points
```

## Annotate cells by lassoing

Lasso a population, label it, and the label is written straight back into
`adata.obs` (or a DataFrame column) — curate cell types interactively:

```python
w = rs.scatterplot(adata, x="X_umap", color_by="leiden")
w                                  # lasso a cluster
w.annotate("cell_type", "T cells") # -> writes adata.obs["cell_type"] for those cells
# lasso another, w.annotate("cell_type", "B cells"), ... then:
rs.scatterplot(adata, x="X_umap", color_by="cell_type")
```

## Differential expression of a selection

Lasso a population and get its top markers vs the rest (or vs another lasso):

```python
w = rs.scatterplot(adata, x="X_umap", color_by="leiden")
w                          # lasso a cluster
w.diff_expression(n=10)    # top genes for the selection vs all other cells
# or two saved selections:
a = w.selection            # after lassoing group A
# (lasso group B)
w.diff_expression(a, w.selection)
```

## Richer tooltips

Show extra fields on hover:

```python
rs.scatterplot(adata, x="X_umap", color_by="leiden",
               tooltip_by=["n_genes", "sample", "CST3"])   # obs cols or genes
```

## Composition of a selection

Lasso a region and see what it's made of:

```python
w = rs.scatterplot(adata, x="X_umap", color_by="leiden")
w                                  # lasso a region
w.composition("leiden")            # -> count + fraction per cluster in the selection
```

## Linked grid

Compare embeddings side by side — pan/zoom and lasso selection stay in sync:

```python
from reglscatterpy import scatterplot, compose

a = scatterplot(adata, x="X_umap", color_by="leiden")
b = scatterplot(adata, x="X_pca",  color_by="leiden")
compose([a, b])            # 2-up grid, linked camera + selection
```

## Toolbar & selection extras

`scatterplot(..., toolbar="left")` (or `"top"`, `"none"`) shows an in-plot
toolbar: pan, lasso, zoom-to-selection, reset, screenshot. Pass
`zoom_on_selection=True` to auto-frame a lasso selection.

Encode a numeric column on point **size** or **opacity** (in addition to
colour): `scatterplot(adata, x="X_umap", color_by="leiden", size_by="n_genes")`
or `opacity_by="total_counts"`.

## Supported objects

| Input | `x` (embedding) | `color_by` / `group_by` |
|-------|-----------------|-------------------------|
| `AnnData` | `obsm` key (`"X_umap"`, `"umap"`, `"spatial"`, …) | `obs` column or `var_names` feature |
| `MuData` | global `obsm` or `"modality:embedding"` | `obs` column or `"modality:feature"` |
| `SpatialData` | table's `obsm` (defaults to `"spatial"`) | table's `obs` / features |
| `pandas.DataFrame` | column name | column name or vector |
| `numpy.ndarray` | column index | vector |

## API parity with R

`rs.scatterplot(...)` mirrors R's `reglScatterplot(...)`: `color_by` / `group_by`,
`point_size`, `opacity`, `point_color`, `pixel_ratio`, `continuous_palette` /
`categorical_palette`, `custom_colors`, `vmin` / `vmax`, `center_zero`,
`filter_by`, legend styling, `enable_download`, and more.

> A `backend="jscatter"` option also exists if you'd rather render with
> [jupyter-scatter](https://github.com/flekschas/jupyter-scatter)
> (`pip install reglscatterpy[render]`); the default native widget is
> recommended.

## The widget bundle

`src/reglscatterpy/static/widget.js` is a **built artifact** (an anywidget ESM
bundle). Its source — the shared rendering widget plus the anywidget adapter —
lives in the **reglScatterplotR** repo under `js/`. To refresh it after a JS
change, build there and copy the result here:

```bash
# from a sibling checkout of reglScatterplotR
cd reglScatterplotR/js && npm install && npm run build
cp dist/widget.js ../../reglscatterpy/src/reglscatterpy/static/widget.js
```

## Develop / test

```bash
pip install -e .[dev]
pytest          # extraction tests skip cleanly without anndata/scipy
```
