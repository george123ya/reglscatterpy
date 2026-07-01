# Widget object

The live, kernel-linked widget returned by [`scatterplot(..., interactive=True)`](api.md#reglscatterpy.scatterplot) (and by every panel inside [`compose`](api.md#reglscatterpy.compose)). These methods need a running kernel; a static/exported plot keeps the *visual* interactions (pan/zoom/lasso/legend/export) but cannot round-trip to Python.

Most methods accept a `selection=` override; when omitted they act on the **current lasso selection** (`w.selection`). Indices are always given and returned in **original data order**.

## Properties

### `w.selection`

Indices of the lasso-selected points (read or assign), always in
**data order** — translated through the draw-order permutation when
the plot was z-ordered (sort_order / random_state).

Live (`interactive=True`) only — on a static plot this stays empty
because there is no kernel link.

### `w.filtered`

Original indices of the cells currently passing the in-plot filters
(range sliders + legend categories). When no filter is active this is
**all shown cells** (everything passes). Live (`interactive=True`)
only — like `selection`, but for the filter instead of the lasso.

### `w.colors`

The categorical colour map as `{category: '#rrggbb'}` (the rendered
palette, including any in-plot recolours via the legend colorpicker).
`None` for a continuous / single-colour plot. Save it scanpy-style with
e.g. `adata.uns['louvain_colors'] = list(w.colors.values())`.

## Methods

### `w.subset(selection=None)`

The source object subset to the selected cells (`adata[w.selection]`).

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `selection` | `list[int]` \| `list[str]` \| `bool mask` \| `None` | Cells to subset. `None` uses the current lasso selection. | `None` |

**Returns** the same type as the source — `adata[sel]` for AnnData/MuData, `df.iloc[sel]` for a DataFrame.

### `w.annotate(key, label, selection=None)`

Write `label` onto the lasso-selected cells in `obs[key]` / column `key`
of the source object. Useful for curating cell types interactively.

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `key` | `str` | `obs` / DataFrame column to write into (created if missing). | *required* |
| `label` | `str` | Value to assign to the selected cells. | *required* |
| `selection` | `list[int]` \| `list[str]` \| `bool mask` \| `None` | Cells to label. `None` uses the current lasso selection. | `None` |

**Returns** the annotated source object (the column becomes a `pandas.Categorical`).

### `w.composition(by, selection=None, normalize=True)`

Count + fraction of the selected cells in each category of `by` — "what is
this region made of?".

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `by` | `str` | Categorical `obs` / DataFrame column to break the selection down by. | *required* |
| `selection` | `list[int]` \| `list[str]` \| `bool mask` \| `None` | Cells to summarise. `None` uses the current lasso selection. | `None` |
| `normalize` | `bool` | Also add a `fraction` column (counts ÷ total). | `True` |

**Returns** a `pandas.DataFrame` indexed by category with a `count` (and, if
`normalize`, `fraction`) column.

### `w.diff_expression(group_a=None, group_b=None, n=10, layer=None, method='wilcoxon', key_added=None, use_raw=None, engine='auto')`

Top differential genes between two cell groups (markers of a lasso vs the
rest, or one lasso vs another).

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `group_a` | `list[int]` \| `list[str]` \| `bool mask` \| `None` | The "test" group. `None` uses the current lasso selection. | `None` |
| `group_b` | `list[int]` \| `list[str]` \| `bool mask` \| `None` | The reference group. `None` means **all other cells** ("rest"). | `None` |
| `n` | `int` | Number of top genes to keep per group. | `10` |
| `layer` | `str` \| `None` | `adata.layers` key to test on. `None` auto-picks a matrix with finite logFC (see note). | `None` |
| `method` | `str` | scanpy test: `'wilcoxon'`, `'t-test'`, `'t-test_overestim_var'`, `'logreg'`. | `'wilcoxon'` |
| `key_added` | `str` \| `False` \| `None` | `adata.uns` key to save under. `None` → `'rank_genes_groups'`; `False` → don't save. | `None` |
| `use_raw` | `bool` \| `None` | Force testing on `adata.raw`. `None` lets the matrix-picker decide. | `None` |
| `engine` | `str` | Compute backend / scanpy requirement — see the [engine table](#which-engine-ran-scanpy-gpu-fallback). | `'auto'` |

When **scanpy** is installed (and the source is an AnnData) this runs
`sc.tl.rank_genes_groups` on a copy; otherwise it falls back to a Welch
t-test. AnnData/MuData only.

**Returns the scanpy-native result** — the `params` + rec.array dict
(`names` / `scores` / `logfoldchanges` / `pvals` / `pvals_adj`), identical to
what scanpy stores in `adata.uns`. The return value and the saved entry are
the same object, and a tidy table is one
`sc.get.rank_genes_groups_df(adata, group="A")` away.

!!! note "Finite logFC"
    A scaled / z-scored `.X` has negative values, which make scanpy emit
    `NaN` logfoldchanges. When you don't pin a `layer` / `use_raw`, the
    matrix-picker routes to log-normalised expression for the **same genes**
    (`adata.raw` restricted to `var_names`), keeping `params` clean
    (`layer: None`, `use_raw: False`).

### `w.diff_expression_by(by, group_a=None, group_b=None, selection=None, n=10, layer=None, method='wilcoxon', key_added=None, use_raw=None, min_cells=2)`

Differential expression **between the levels of an obs column**, restricted to
the lasso selection — e.g. lasso a region, then compare `condition` or `time`
*within* it.

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `by` | `str` | The `obs` column whose levels you compare (e.g. `"condition"`, `"time"`). | *required* |
| `group_a` | `str` \| `None` | A single level of `by`. With `group_b` → A-vs-B; alone → A-vs-rest. | `None` |
| `group_b` | `str` \| `None` | A second level of `by` to use as the reference. | `None` |
| `selection` | `list[int]` \| `list[str]` \| `bool mask` \| `None` | Cells to restrict to. `None` uses the lasso; empty selection → **all** cells. | `None` |
| `n` | `int` | Number of top genes to keep per group. | `10` |
| `layer` | `str` \| `None` | `adata.layers` key to test on. `None` auto-picks a finite-logFC matrix. | `None` |
| `method` | `str` | scanpy test: `'wilcoxon'`, `'t-test'`, `'t-test_overestim_var'`, `'logreg'`. | `'wilcoxon'` |
| `key_added` | `str` \| `False` \| `None` | `adata.uns` key to save under. `None` → `'rank_genes_groups'`; `False` → don't save. | `None` |
| `use_raw` | `bool` \| `None` | Force testing on `adata.raw`. `None` lets the matrix-picker decide. | `None` |
| `min_cells` | `int` | In the all-levels mode, skip (with a warning) any level with fewer cells than this in the selection. | `2` |
| `engine` | `str` | Compute backend / scanpy requirement — see the [engine table](#which-engine-ran-scanpy-gpu-fallback). | `'auto'` |

The three modes:

- `group_a` **and** `group_b` → a single A-vs-B comparison (e.g. `group_a="D30", group_b="Y1"`).
- `group_a` only → that level vs the pooled rest of the selection.
- **neither** → every present level vs the rest, run as a **single** scanpy `rank_genes_groups` call (one-vs-rest, the scanpy idiom).

**Returns the scanpy-native result** in every mode — same shape as
`diff_expression`. The *real* `by` name and *real* level names appear in
`params` / the rec.array fields (no synthetic `_rs_by` / `_rs_lognorm`), so
`sc.pl.rank_genes_groups(adata)` and
`sc.get.rank_genes_groups_df(adata, group=...)` work straight after. The same
object is **auto-saved** to `adata.uns` (`key_added` controls the key). Also
exposed on a [`compose`](api.md#reglscatterpy.compose) grid
(`grid.diff_expression_by(...)`, proxied to a shared panel).

```python
import scanpy as sc
w.selection = some_cells                     # or lasso in the UI

# every condition in the lasso vs the rest -> scanpy-native result (== adata.uns):
res = w.diff_expression_by("condition")       # params + rec.arrays, one col per level
res["names"].dtype.names                      # ('D30', 'Y1', 'Y2')
sc.get.rank_genes_groups_df(adata, group="D30")   # tidy table for any level

# one specific pair, saved under a custom uns key:
res = w.diff_expression_by("condition", group_a="D30", group_b="Y1",
                           key_added="D30_vs_Y1")
# one level vs the rest of the lasso:
res = w.diff_expression_by("condition", group_a="Y1")
# also available straight off the selection:
res = w.selection.diff_expression_by("time")
```

#### Which engine ran (scanpy / GPU / fallback)

Both DE methods take `engine=` to pick the backend **and** to control whether
scanpy is *required* — so a Welch-t-test fallback is never silent:

| `engine` | Backend | If unavailable |
|----------|---------|----------------|
| `'auto'` *(default)* | scanpy (`sc.tl.rank_genes_groups`) if importable | **warns**, then CPU Welch t-test |
| `'scanpy'` | scanpy — always | raises `ImportError` (never falls back) |
| `'gpu'` (alias `'cupy'`) | **GPU** Welch t-test via `cupy` (per-cell reductions on the GPU) | raises `ImportError` |
| `'rapids'` | **GPU** logistic-regression DE via [rapids-singlecell](https://rapids-singlecell.readthedocs.io/) | raises `ImportError` |
| `'ttest'` | the built-in CPU Welch t-test | — |

Two ways to be sure scanpy ran:

1. **Require it** — `engine="scanpy"` raises instead of falling back silently.
2. **Check the result** — scanpy records its correction in `params`:

```python
res = w.diff_expression(engine="scanpy")
res["params"]["method"]        # -> 'wilcoxon'  (your real test)
res["params"]["corr_method"]   # -> 'benjamini-hochberg'  (scanpy's default)
# both Welch t-test paths (CPU + GPU) instead stamp method='t-test', corr_method='none'
```

!!! note "GPU acceleration"
    `sc.tl.rank_genes_groups` is **CPU-only**, so GPU DE uses a different test:

    * **`engine="gpu"`** (needs only `cupy`, `pip install cupy-cuda12x`) runs a
      Welch **t-test** with the heavy `O(cells × genes)` mean/variance reductions
      on the GPU. The formula matches scipy's `ttest_ind(..., equal_var=False)`,
      so GPU and CPU results agree to floating point. On an RTX 2060, DE over
      **120k cells × 2000 genes** ran in **~1.9 s** vs **~36 s** on the CPU
      (≈19×). It pays off when a group has many cells; for a few-thousand-cell
      lasso the CPU path is already fast.
    * **`engine="rapids"`** uses [rapids-singlecell](https://rapids-singlecell.readthedocs.io/)
      (full RAPIDS / `cuml`) for GPU **logreg** DE — `method=` is ignored. Heavier
      install, but returns the same scanpy-native `adata.uns` dict.

### `w.highlight(indices, color=None)`

Persistently mark points with a crisp ring + size bump (the engine's
selection look) — but this is **not** the selection, so it survives a
double-click and a new lasso.

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `indices` | `list[int]` \| `[]` \| `None` | Original data indices to mark. `[]` / `None` clears the highlight. | *required* |
| `color` | `str` \| `None` | Ring colour (hex / CSS colour). `None` keeps the default. | `None` |

Live (`interactive=True`) only — needs the kernel link.

!!! note
    With `progressive=True` the highlight marks the currently-shown cells; it
    doesn't yet follow new cells streamed in on zoom.

### `w.morph_to(basis, duration=1200)`

Animate the points from their current layout to another embedding — e.g.
`w.morph_to('spatial')` morphs a UMAP into the spatial layout (and back with
`w.morph_to('umap')`). Positions tween; colours and sizes stay.

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `basis` | `str` | Target `obsm` embedding to morph into (`'umap'`, `'spatial'`, …). | *required* |
| `duration` | `int` | Animation length in milliseconds. | `1200` |

Live (`interactive=True`) only; not for `progressive=True` (only a subset is
resident). Returns `self` so calls can chain.

### `w.to_html(path, title='reglscatterpy plot')`

Save this plot as a standalone, offline HTML file (like R's
`htmlwidgets::saveWidget`) — inlines the bundle and data, stays interactive
with no kernel.

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `path` | `str` \| `Path` | Where to write the `.html` file. | *required* |
| `title` | `str` | `<title>` of the saved page. | `'reglscatterpy plot'` |

On a live (`interactive=True`) widget it captures the **current** view /
selection / filter, so the file opens exactly as you left it.

### `w.update(spec)`

Push a partial spec update to a live widget (e.g. recolour, change point size)
without rebuilding it.

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `spec` | `dict` | Plot properties to merge into the current spec. | *required* |

Returns `self`. Live (`interactive=True`) only.
