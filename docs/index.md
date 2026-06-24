# reglscatterpy

Interactive **WebGL scatterplots** for single-cell / spatial data in Python —
**AnnData, MuData, SpatialData**, pandas and numpy. Renders millions of points
in the browser via [`regl-scatterplot`](https://github.com/flekschas/regl-scatterplot),
in **Jupyter, JupyterLab, VS Code and Colab**.

<p align="center">
  <img src="https://raw.githubusercontent.com/george123ya/reglscatterpy/main/assets/demo.gif"
       alt="Panning, lassoing and legend-filtering an interactive UMAP" width="760">
</p>

This is the Python companion to the R package
[**reglScatterplotR**](https://github.com/george123ya/reglScatterplotR) — both
drive the *same* compiled widget, so a plot looks and behaves identically across
R and Python.

## Install

```bash
pip install reglscatterpy
pip install anndata          # for AnnData; mudata / spatialdata as needed
```

## 30-second start

```python
import scanpy as sc
import reglscatterpy as rs

adata = sc.datasets.pbmc3k_processed()
rs.scatterplot(adata, basis="umap", color="louvain")   # an obs column
rs.scatterplot(adata, basis="umap", color="CST3")      # a gene
```

## Where to next

<div class="grid cards" markdown>

- :material-book-open-variant: **[User guide](guide.md)** — every feature with
  runnable examples: big-data rendering, selection round-trip, annotate, DE,
  linked grids, HTML export, embedding animations.

- :material-function: **[Functions API](api.md)** — `scatterplot()` and its 80+
  arguments, plus `compose`, `save_html`, `record_html`, `extract`.

- :material-cursor-default-click: **[Widget object](widget-api.md)** — the live
  widget's methods and properties (`selection`, `annotate`, `diff_expression`,
  `composition`, `highlight`, `morph_to`, `to_html`).

</div>

## Feature highlights

- **Atlas-scale**: `max_points="auto"` density-subsamples while keeping rare
  cell types; `max_points=None` puts every cell on the GPU; `progressive=True`
  does detail-on-zoom beyond ~4M cells.
- **Selection round-trip**: lasso in the browser, read `w.selection` back in
  Python; `annotate`, `diff_expression`, `composition` on the selection.
- **Linked grids**: `color=[...]` or `basis=[...]` builds a synced multi-panel
  grid; `compose([...])` links arbitrary plots (camera + selection + filters).
- **Embedding animation**: `w.morph_to("spatial")` tweens points (and axes)
  from one embedding to another, preserving filter/selection/colour.
- **Offline export**: `save_html` writes a self-contained interactive `.html`;
  `record_html` turns a whole notebook into a static report with no re-run.
