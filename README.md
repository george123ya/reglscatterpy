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

Plots default to a fixed **700 px** width (like matplotlib / plotly); pass
`width=` for another size or `width=None` to fill the cell.

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

### Backends

By default `scatterplot()` renders reglscatterpy's own widget (the shared one).
`backend="jscatter"` is an optional, lighter alternative that renders through
[jupyter-scatter](https://github.com/flekschas/jupyter-scatter) instead, without
this package's legend / filter / export UI (`pip install reglscatterpy[render]`).

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
