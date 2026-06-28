# Widget object
The live, kernel-linked widget returned by [`scatterplot(..., interactive=True)`](api.md#reglscatterpy.scatterplot) (and by every panel inside [`compose`](api.md#reglscatterpy.compose)). These methods need a running kernel; a static/exported plot keeps the *visual* interactions (pan/zoom/lasso/legend/export) but cannot round-trip to Python.

## Properties

### `w.selection`

Indices of the lasso-selected points (read or assign), always in
**data order** — translated through the draw-order permutation when
the plot was z-ordered (sort_order / random_state).

Live (``interactive=True``) only — on a static plot this stays empty
because there is no kernel link.

### `w.filtered`

Original indices of the cells currently passing the in-plot filters
(range sliders + legend categories). When no filter is active this is
**all shown cells** (everything passes). Live (``interactive=True``)
only — like :attr:`selection`, but for the filter instead of the lasso.

### `w.colors`

The categorical colour map as ``{category: '#rrggbb'}`` (the rendered
palette, including any in-plot recolours via the legend colorpicker).
``None`` for a continuous / single-colour plot. Save it scanpy-style with
e.g. ``adata.uns['louvain_colors'] = list(w.colors.values())``.

## Methods

### `w.subset(selection=None)`

The source object subset to the selected cells (``adata[w.selection]``).

### `w.annotate(key, label, selection=None)`

Write ``label`` onto the lasso-selected cells in ``obs[key]`` /
column ``key`` of the source object. Returns the annotated object.

### `w.composition(by, selection=None, normalize=True)`

Count + fraction of the selected cells in each category of ``by``.

### `w.diff_expression(group_a=None, group_b=None, n=10, layer=None, method='wilcoxon', key_added=None, use_raw=None)`

Top differential genes between two cell groups.

``group_a`` defaults to the lasso selection; ``group_b`` to the rest.
Groups accept integer positions, obs_names, or a boolean mask. When
**scanpy** is installed (and the source is an AnnData) this runs
``sc.tl.rank_genes_groups`` on a copy and returns its result frame
(names / scores / logfoldchanges / pvals / pvals_adj). Otherwise it
falls back to a Welch t-test. AnnData/MuData only.

When the source is an **AnnData** the result is **auto-saved** to
``adata.uns`` (scanpy-style) — default key ``"rank_genes_groups"`` (the
scanpy convention), or ``key_added`` if you pass one. Pass
``key_added=False`` to skip saving.

### `w.diff_expression_by(by, group_a=None, group_b=None, selection=None, n=10, layer=None, method='wilcoxon', key_added=None, use_raw=None, min_cells=2)`

Differential expression **between the levels of an obs column**, restricted to the
lasso selection.

Lasso a group of cells, then split them by ``by`` (an ``obs`` column such as
``"time"`` or ``"condition"``) and compare its levels:

- ``group_a`` **and** ``group_b`` given → a single A-vs-B comparison
  (e.g. ``group_a="D30", group_b="Y1"``); returns one DataFrame.
- ``group_a`` only → that level vs the pooled rest of the selection; one DataFrame.
- **neither** → every present level vs the rest, run as a **single** scanpy
  ``rank_genes_groups`` call (one-vs-rest, the scanpy idiom); returns a ``dict``
  ``{level: df}``.

Whichever the mode, the result is saved to ``adata.uns`` as **one clean
scanpy-native entry** — the *real* ``by`` name and *real* level names in
``params`` / columns (no synthetic ``_rs_by`` / ``_rs_lognorm``), so
``sc.pl.rank_genes_groups(adata)`` works straight after. Default key
``"rank_genes_groups"``; pass ``key_added`` to choose the key, or
``key_added=False`` to skip saving.

Cells default to the current lasso ``selection`` (pass ``selection=`` to override:
integer positions / obs_names / a boolean mask); if nothing is selected it falls
back to **all** cells. In the all-levels mode, levels with fewer than ``min_cells``
cells in the selection are skipped (with a warning); an explicit ``group_a`` /
``group_b`` is always honoured. Uses ``sc.tl.rank_genes_groups`` when scanpy is
installed (same finite-logFC matrix routing as ``diff_expression``), else a Welch
t-test. AnnData/MuData only. Also exposed on a [`compose`](api.md#reglscatterpy.compose)
grid (``grid.diff_expression_by(...)``, proxied to a shared panel).

```python
w.selection = some_cells                     # or lasso in the UI
# every condition in the lasso vs the rest (one clean adata.uns entry):
res = w.diff_expression_by("condition")       # {"D30": df, "Y1": df, "Y2": df}
# one specific pair, saved under a custom uns key:
df  = w.diff_expression_by("condition", group_a="D30", group_b="Y1",
                           key_added="D30_vs_Y1")
# one level vs the rest of the lasso:
df  = w.diff_expression_by("condition", group_a="Y1")
# also available straight off the selection:
res = w.selection.diff_expression_by("time")
```

### `w.highlight(indices, color=None)`

Persistently mark points with a crisp ring + size bump (the engine's
selection look) — but this is **not** the selection, so it survives a
double-click and a new lasso. ``indices`` are original data indices;
``color`` sets the ring colour (a hex / CSS colour). Pass ``[]`` / ``None``
to clear. Live (``interactive=True``) only — needs the kernel link.

Note: with ``progressive=True`` the highlight marks the currently-shown
cells; it doesn't yet follow new cells streamed in on zoom.

### `w.morph_to(basis, duration=1200)`

Animate the points from their current layout to another embedding —
e.g. ``w.morph_to('spatial')`` morphs a UMAP into the spatial layout
(and back with ``w.morph_to('umap')``). Positions tween; colours and
sizes stay. Live (``interactive=True``) only; not for ``progressive=True``
(only a subset is resident). Returns ``self`` so calls can chain.

### `w.to_html(path, title='reglscatterpy plot')`

Save this plot as a standalone, offline HTML file (like R's
``htmlwidgets::saveWidget``) — inlines the bundle and data, stays
interactive with no kernel. On a live (``interactive=True``) widget it
captures the CURRENT view / selection / filter so the file opens exactly
as you left it.

### `w.update(spec: 'dict') -> "'ReglScatter'"`

_(undocumented)_
