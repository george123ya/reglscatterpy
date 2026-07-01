# reglscatterpy

[![PyPI](https://img.shields.io/pypi/v/reglscatterpy.svg)](https://pypi.org/project/reglscatterpy/)
[![Python versions](https://img.shields.io/pypi/pyversions/reglscatterpy.svg)](https://pypi.org/project/reglscatterpy/)
[![Docs](https://img.shields.io/badge/docs-github.io-blue.svg)](https://george123ya.github.io/reglscatterpy/)
[![Live demo](https://img.shields.io/badge/live%20demo-pbmc3k-brightgreen.svg)](https://george123ya.github.io/reglscatterpy/demo/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Fast, interactive scatterplots for single-cell and spatial data — right in your notebook.**

Plot millions of cells, then *use your mouse*: pan, zoom, **lasso a population**,
toggle cell types in the legend — and read your selection straight back into
Python. It speaks **AnnData, MuData, SpatialData**, pandas and numpy, and works
in **Jupyter, JupyterLab, VS Code and Colab**.

<p align="center">
  <img src="https://raw.githubusercontent.com/george123ya/reglscatterpy/main/assets/demo.gif"
       alt="Panning, lassoing and legend-filtering an interactive UMAP" width="760">
</p>

▶️ **[Try the live demo (real pbmc3k, no install) →](https://george123ya.github.io/reglscatterpy/demo/)**
 ·  📖 **[Full docs & API reference →](https://george123ya.github.io/reglscatterpy/)**

Under the hood it renders with [`regl-scatterplot`](https://github.com/flekschas/regl-scatterplot)
(WebGL), and it's the Python twin of the R package
[**reglScatterplotR**](https://github.com/george123ya/reglScatterplotR) — both
drive the *same* compiled widget, so a plot looks and behaves identically in R
and Python.

## Install

```bash
pip install reglscatterpy            # numpy, pandas, anywidget
pip install anndata                  # for AnnData; mudata / spatialdata as needed
```

## A 60-second tour

Give it an AnnData, say which embedding to show (`basis=`) and what to colour by
— an `obs` column **or a gene name**. That's it:

```python
import scanpy as sc
import reglscatterpy as rs

adata = sc.datasets.pbmc3k_processed()
rs.scatterplot(adata, basis="umap", color_by="louvain")   # colour by cell type
rs.scatterplot(adata, basis="umap", color_by="CST3")      # …or by a gene
```

`basis=` accepts short names like `"umap"` / `"pca"` (resolved to the `obsm` key
`X_umap`, … — case-insensitive). For a plain **DataFrame**, give the coordinate
columns with `x=` / `y=` instead:

```python
import numpy as np, pandas as pd
df = pd.DataFrame({"x": np.random.rand(10_000), "y": np.random.rand(10_000),
                   "ct": np.random.choice(list("ABC"), 10_000)})
rs.scatterplot(df, x="x", y="y", color_by="ct")
```

> Plots are **700 px wide by default** (not the full cell width). Pass `width=`
> (pixels) for another size, or `width=None` to fill the cell.

## Gallery

| Categorical colouring | Continuous (gene) colouring |
|---|---|
| ![Categorical UMAP with frosted legend](https://raw.githubusercontent.com/george123ya/reglscatterpy/main/assets/umap-categorical.png) | ![Gene-expression UMAP with colour bar](https://raw.githubusercontent.com/george123ya/reglscatterpy/main/assets/umap-continuous.png) |
| **`filter_by` distribution sliders** | **Linked grid (`compose`)** |
| ![Range-filter sliders with histograms](https://raw.githubusercontent.com/george123ya/reglscatterpy/main/assets/filter-sliders.png) | ![Two embeddings with synced camera and selection](https://raw.githubusercontent.com/george123ya/reglscatterpy/main/assets/linked-grid.png) |

The UMAP panels are rendered from **real pbmc3k** (2,638 PBMCs). You can
regenerate the figures and the
[live demo](https://george123ya.github.io/reglscatterpy/demo/) with the
maintainer scripts — see [`scripts/README.md`](https://github.com/george123ya/reglscatterpy/blob/main/scripts/README.md).

## The fun part: lasso, then ask questions

The point of an interactive plot is the **round-trip** — circle some cells and
keep working with them in Python. Pass `interactive=True` to get the live,
kernel-linked widget (the default plot is a portable snapshot; see
[Static vs. interactive](#static-by-default-interactive-on-request)):

```python
w = rs.scatterplot(adata, basis="umap", color_by="louvain", interactive=True)
w                          # show it, then lasso a population in the browser
```

**Read the selection back** — or set it from Python:

```python
w.selection                # -> [12, 87, 134, ...]  positional indices
adata[w.selection]         # subset the AnnData directly
sub = w.subset()           # same thing, convenience

w.selection = list(range(100))   # drive it from Python to highlight points
```

**Label cells you lassoed** (curate cell types by hand) — the label is written
straight into `adata.obs`:

```python
w.annotate("cell_type", "T cells")    # writes adata.obs["cell_type"] for those cells
# lasso another cluster, w.annotate("cell_type", "B cells"), …then re-plot:
rs.scatterplot(adata, basis="umap", color_by="cell_type")
```

**Get the top marker genes** of a selection vs the rest (or vs another lasso).
Both DE calls **return scanpy's native result** and auto-save it to
`adata.uns`, so the rest of the scanpy world just works:

```python
import scanpy as sc
w.diff_expression(n=10)                        # selection "A" vs rest -> adata.uns["rank_genes_groups"]
sc.get.rank_genes_groups_df(adata, group="A")  # tidy table when you want one
sc.pl.rank_genes_groups(adata)                 # scanpy's own plot, straight after

a = w.selection                                # save group A
# (lasso group B)
w.diff_expression(a, w.selection)              # A vs B
```

Or split a single lasso **by an `obs` column** (e.g. `condition` / `time`) and
compare its levels — one clean `adata.uns` entry, one column per level:

```python
res = w.diff_expression_by("condition")        # each level vs the rest
res["names"].dtype.names                        # ('D30', 'Y1', 'Y2')
w.diff_expression_by("condition", group_a="D30", group_b="Y1")   # one pair
w.diff_expression_by("condition", key_added="cond_de")           # choose the uns key
```

> **Which engine?** By default DE uses scanpy when it's installed and *warns* if
> it has to fall back. Pass `engine="scanpy"` to **require** it (raises instead
> of falling back), or run it on the **GPU**: `engine="gpu"` needs only `cupy`
> (~19× faster on 120k cells — a Welch t-test with the reductions on the GPU),
> and `engine="rapids"` uses [rapids-singlecell](https://rapids-singlecell.readthedocs.io/)
> for GPU logreg. See the [Widget API](https://george123ya.github.io/reglscatterpy/widget-api/#which-engine-ran-scanpy-gpu-fallback)
> for the full table.

**See what a region is made of:**

```python
w.composition("louvain")           # count + fraction per cluster in the selection
```

## Spatial data

Anything with coordinates in `obsm["spatial"]` plots the same way — just point
`basis` at it. Great for **Visium, Xenium, MERFISH, CosMx**, … (see the
[squidpy Xenium tutorial](https://squidpy.readthedocs.io/en/latest/notebooks/tutorials/tutorial_xenium.html)
for getting data in):

```python
import squidpy as sq
import reglscatterpy as rs

adata = sq.datasets.visium_hne_adata()         # public mouse-brain Visium
rs.scatterplot(adata, basis="spatial", color_by="cluster")    # tissue map
rs.scatterplot(adata, basis="spatial", color_by="Olfm1")      # a gene, in situ
```

A Xenium / MERFISH run read with `sq.read.xenium(...)` (or any AnnData with a
`spatial` embedding) works identically — bump `point_size` for sparse tissue:

```python
rs.scatterplot(adata, basis="spatial", color_by="leiden", point_size=3)
```

And because both layouts live in the same object, you can **morph between them**
— watch clusters snap from a UMAP into their tissue positions:

```python
w = rs.scatterplot(adata, basis="umap", color="leiden", interactive=True)
w                                  # the UMAP
w.morph_to("spatial")              # animate UMAP -> tissue coordinates
w.morph_to("umap", duration=800)   # …and back (ms)
```

The active filter, legend selection, colours and lasso are **kept through the
transition**. (Atlas-scale `progressive=True` plots don't morph — the point set
isn't fully resident.)

## Small-multiples and linked views

Colour one embedding by **several genes / columns at once** — pass a list and
you get a linked grid, one panel per value, with camera + lasso kept in sync:

```python
rs.scatterplot(adata, basis="umap", color_by=["louvain", "CST3", "NKG7"])
```

> A *list of names* means "one panel per name". A raw per-point colour vector
> must be a numpy array / pandas Series (not a list of strings).

Or compose pre-built plots — e.g. compare embeddings side by side. `compose()`
auto-upgrades plain plots to live widgets, so you don't need `interactive=True`
on each:

```python
from reglscatterpy import scatterplot, compose

a = scatterplot(adata, basis="umap", color_by="leiden")
b = scatterplot(adata, basis="pca",  color_by="leiden")
compose([a, b])            # 2-up grid, linked camera + selection
```

A lasso or legend filter on one panel propagates to the others **by original
cell** — even across panels coloured by different variables, and even when the
panels are `progressive=True`. A view reset on one panel resets the group. A
**list of bases** works the same way, and `basis` × `color` can both be lists
(the basis × colour grid, capped at 16 panels):

```python
rs.scatterplot(adata, basis=["umap", "tsne", "pca"], color="leiden")
```

## Atlas-scale rendering

By default `scatterplot()` keeps huge datasets interactive **without silently
hiding cells**, via `max_points` (default `"auto"`):

```python
# AUTO (default): caps at 500k with a density-preserving subsample that KEEPS
# rare cell types. The plot stays honest: an on-figure "500,000 of 3,900,000
# shown" caption, a one-time warning, and w.selection still indexes ALL rows.
rs.scatterplot(adata, basis="umap", color_by="cell_type")

# ALL POINTS RESIDENT (the Allen ABC-Atlas method): every cell on the GPU,
# camera-only pan/zoom — smooth up to ~4M cells on a decent GPU.
rs.scatterplot(adata, basis="umap", color_by="cell_type", max_points=None)
```

For datasets **larger than ~4M**, use `progressive=True` — a light density
overview that re-renders **all cells inside the viewport** as you zoom in (no
preprocessing, lasso stays complete):

```python
rs.scatterplot(adata, basis="umap", color_by="cell_type", progressive=True,
               progressive_opts={"detail_max_points": 300_000, "overscan": 0.6})
```

- `detail_max_points` — max points per zoomed-in viewport (lower = smoother pan).
- `overscan` — margin fetched around the view so panning has no hard cuts.

Rule of thumb: `max_points=None` for ~2–4M real atlases; `progressive=True`
beyond that. The two subsample modes are `"density"` (default, keeps rare types)
and `"random"` (uniform fallback).

## Static by default, interactive on request

By default a plot renders as a **self-contained snapshot** (a sandboxed
`<iframe>` with the WebGL bundle and data baked in) — like a plotly figure, it
shows in JupyterLab, Notebook 7, VS Code and Colab, and **survives reopening the
notebook with no kernel**. It stays fully interactive *visually* — pan, zoom,
lasso, legend, tooltips, PNG/SVG/PDF export — but, having no kernel link, it
can't send a selection back to Python.

For the **Python round-trip** (`w.selection`, `annotate`, `diff_expression`,
linked `compose` grids) pass `interactive=True`:

```python
w = rs.scatterplot(adata, basis="umap", color_by="leiden", interactive=True)
w                          # lasso some cells…
adata[w.selection]         # …read them back in Python
```

> Use the **default** for figures you want to keep/share, and `interactive=True`
> while you're actively selecting (it needs a running kernel and, like any
> Jupyter widget, may show blank on reopen).

## Save a standalone HTML (offline, kernel-free)

The Python equivalent of R's `htmlwidgets::saveWidget`: one self-contained
`.html` that **inlines the widget and the plot's data**, so it opens in any
browser with no kernel and no internet:

```python
w = rs.scatterplot(adata, x="X_umap", color_by="leiden")
rs.save_html(w, "umap.html")      # or:  w.to_html("umap.html")
```

It's fully interactive (pan/zoom, legend, lasso, tooltips, PNG/SVG/PDF export)
but a **snapshot** — no kernel, so the Python round-trips only work in the live
notebook. The bundle is inlined gzip-compressed (~0.5 MB), so a one-plot file is
well under 1 MB.

### A whole notebook → one HTML report (no re-running)

Plain `jupyter nbconvert --to html` leaves the plots blank. The fix that
**avoids re-executing a heavy notebook** is *record mode*: call
`rs.record_html()` once at the top, run your notebook normally, and each plot
bakes a static interactive copy into its own cell output. After that:

```python
import reglscatterpy as rs
rs.record_html()                 # run once near the top, then work as usual
# ... rs.scatterplot(...) cells ...
```

```bash
# now either of these makes a report WITHOUT re-running anything:
jupyter nbconvert --to html analysis.ipynb
reglscatterpy-report analysis.ipynb -o analysis_report.html
```

`reglscatterpy-report` (and `rs.save_notebook_html(...)`) default to **not**
re-executing — they reuse recorded outputs and share **one** copy of the bundle
across all plots. Pass `--execute` / `execute=True` to re-run a notebook that
wasn't recorded. Needs `nbconvert` + `ipykernel`
(`pip install 'reglscatterpy[report]'`).

> Recorded plots are a one-way snapshot — `w.selection` / `w.annotate` no longer
> round-trip. Call `rs.record_html(False)` to go back to the live widget.

## Theme (light / dark / auto)

Plots are a **white "figure card" by default** — portable, and matching the
exported HTML. To make the live widget follow your notebook, pass `theme=`:

```python
rs.scatterplot(adata, basis="umap", color="leiden", theme="auto")   # dark in a dark theme
rs.scatterplot(adata, basis="umap", color="leiden", theme="dark")   # always dark
```

- `"light"` (default) — white card.
- `"dark"` — dark card with light axes/legend.
- `"auto"` (alias `"system"`) — dark **only** when the host (VS Code /
  JupyterLab) is in a dark theme.

Set it once per session with `rs.set_theme("auto")` (`rs.get_theme()` reads it
back). A per-call `theme=` overrides the global; an explicit `background_color=`
/ `axis_color=` always wins. The theme only affects the **live** widget; an
exported `.html` stays portably light.

## More knobs worth knowing

**Richer tooltips** — show extra fields (obs columns or genes) on hover:

```python
rs.scatterplot(adata, x="X_umap", color_by="leiden",
               tooltip_by=["n_genes", "sample", "CST3"])
```

**Outlines & highlighting** — `add_outline=True` rings *every* point
(scanpy-style, small/medium plots); `w.highlight([...])` persistently marks a
**chosen subset** that survives double-click and new lassoes:

```python
rs.scatterplot(adata, basis="umap", color_by="cluster", add_outline=True)

w = rs.scatterplot(adata, basis="umap", color_by="cluster", interactive=True)
w.highlight([12, 87, 134], color="red")   # ring + size bump on these cells
w.highlight([])                            # clear
```

**Encode more variables** — put a numeric column or a gene on point **size** or
**opacity** too:

```python
rs.scatterplot(adata, basis="umap", color_by="leiden", size_by="n_genes")  # or size_by="CST3"
rs.scatterplot(adata, basis="umap", color_by="leiden", opacity_by="total_counts")
```

**Toolbar & legend** — `toolbar="left"` (or `"top"`, `"none"`) shows an in-plot
toolbar (pan, lasso, zoom-to-selection, reset, screenshot); `zoom_on_selection=True`
auto-frames a lasso; the categorical legend shows a live per-category count
(`legend_counts=False` to hide); `alpha=` sets a global transparency.

## Supported objects

| Input | `x` (embedding) | `color_by` / `group_by` |
|-------|-----------------|-------------------------|
| `AnnData` | `obsm` key (`"X_umap"`, `"umap"`, `"spatial"`, …) | `obs` column or `var_names` feature |
| `MuData` | global `obsm` or `"modality:embedding"` | `obs` column or `"modality:feature"` |
| `SpatialData` | table's `obsm` (defaults to `"spatial"`) | table's `obs` / features |
| `pandas.DataFrame` | column name | column name or vector |
| `numpy.ndarray` | column index | vector |

## Same plot in R and Python

`rs.scatterplot(...)` mirrors R's `reglScatterplot(...)`: `color_by` / `group_by`,
`point_size`, `opacity`, `point_color`, `pixel_ratio`, `continuous_palette` /
`categorical_palette`, `custom_colors`, `vmin` / `vmax`, `center_zero`,
`filter_by`, legend styling, `enable_download`, and more. Equivalence is locked
down by `tests/test_payload_parity.py`, which checks the Python payload
byte-for-byte against R fixtures.

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
