"""The public ``scatterplot()`` entry point.

By default this renders through reglscatterpy's **own** anywidget, which drives
the very same compiled widget the R package ships (legend, lasso, tooltips,
sync, PNG/SVG/PDF export) - so a plot is pixel-identical across R and Python.
It works in Jupyter, JupyterLab, VS Code, Colab and (via ``shinywidgets``)
Shiny for Python, and embeds into self-contained HTML.

Pass ``backend="jscatter"`` to render through `jupyter-scatter
<https://github.com/flekschas/jupyter-scatter>`_ instead (the JS author's
binding) - lighter, but without the package's own legend / sync / export UI.

The signature mirrors the R ``reglScatterplot()`` so the two feel like one tool.
"""

from __future__ import annotations

from typing import Any, Optional, Sequence, Union

import numpy as np
import pandas as pd

from ._extract import (
    ColorSpec,
    PlotData,
    _is_anndata,
    _is_mudata,
    _is_spatialdata,
    _not_found,
    _resolve_anndata_vec,
    extract,
)
from ._palettes import CONTINUOUS, QUALITATIVE
from ._payload import build_payload

__all__ = ["scatterplot"]

_TOOLBAR_CHOICES = (None, "left", "top", "none")
_DEFAULT_MAX_POINTS = 500_000   # auto-cap for smooth exploration of huge data


class _Unset:
    def __repr__(self):
        return "<unset>"


_UNSET = _Unset()   # sentinel: tells "scanpy alias not passed" from "passed None"


def _is_name_list(spec) -> bool:
    """True for a list/tuple of strings (= panel names for a color-by grid),
    as opposed to a raw per-point vector (which must be an ndarray/Series)."""
    return (
        isinstance(spec, (list, tuple))
        and len(spec) > 0
        and all(isinstance(e, str) for e in spec)
    )


def _viewport_payload(widget, bounds, pad=0.6):
    """Re-render the cells inside the current view (detail-on-zoom). Zoomed in,
    the viewport holds few cells so we draw them ALL (full detail + complete
    lasso); zoomed out it falls back to the density budget. Maps the rendered
    indices back to ORIGINAL rows so w.selection stays correct.

    ``pad`` overscans the fetched region by that fraction of the view span on
    every side, so a pan within the margin doesn't expose empty 'hard cuts'
    before the next refresh arrives."""
    vp = getattr(widget, "_vp", None)
    if not vp:
        return
    x0, y0, x1, y1 = (float(b) for b in bounds)
    xr, yr = vp.get("xrange"), vp.get("yrange")
    # Zoomed back out to ~the full extent? Push the cached overview channels
    # instantly (no recompute, no rebuild) — and via plot.draw, not a re-render.
    if (xr and xr[0] is not None
            and (x1 - x0) >= 0.8 * (xr[1] - xr[0])
            and (y1 - y0) >= 0.8 * (yr[1] - yr[0])):
        if vp.get("showing_overview"):
            return                              # already at overview -> no redundant redraw
        vp["showing_overview"] = True
        widget._inv_draw_order = None
        widget._draw_order = vp["overview_draw_order"]
        widget._source = vp["data"]
        widget.send({"type": "vp_overview",     # JS redraws its cached overview buffers
                     "select": _sel_positions(vp.get("sel"), vp["overview_draw_order"])})
        return
    vp["showing_overview"] = False
    dx, dy = (x1 - x0) * pad, (y1 - y0) * pad      # overscan margin
    x0, x1, y0, y1 = x0 - dx, x1 + dx, y0 - dy, y1 + dy
    fx, fy = vp["x"], vp["y"]
    cand = _grid_query(vp["grid"], x0, y0, x1, y1) if vp.get("grid") is not None else None
    if cand is not None:                           # deep zoom: grid candidates + exact filter
        m = ((fx[cand] >= x0) & (fx[cand] <= x1)
             & (fy[cand] >= y0) & (fy[cand] <= y1))
        in_idx = cand[m]
    else:                                          # wide view / no grid: vectorised full scan
        m = (fx >= x0) & (fx <= x1) & (fy >= y0) & (fy <= y1)
        in_idx = np.where(m)[0]
    if in_idx.size == 0:
        return
    data = vp["data"]
    sub = data[in_idx] if hasattr(data, "obs") else data.iloc[in_idx]
    extra = {"max_points": vp["budget"]}
    if xr and xr[0] is not None:                   # keep the overview domain
        extra["xrange"] = vp["xrange"]
        extra["yrange"] = vp["yrange"]
    if vp.get("point_size") is not None:           # keep size constant
        extra["point_size"] = vp["point_size"]
    if vp.get("vmin") is not None:                 # keep the colour scale constant
        extra["vmin"] = vp["vmin"]
        extra["vmax"] = vp["vmax"]
    if vp.get("categories") is not None:           # keep categorical colours/codes constant
        extra["_color_categories"] = vp["categories"]
    w2 = scatterplot(sub, **{**vp["base"], **extra}, **vp["bk"])
    do = w2._draw_order
    orig = in_idx if do is None else in_idx[np.asarray(do)]
    widget._inv_draw_order = None
    widget._draw_order = orig
    widget._source = data                          # keep selection/DE on the full object
    # Swap points on the EXISTING plot (custom message -> plot.draw); never touches
    # _spec, so no re-mount / spinner / plot recreation. Re-map the persisted lasso
    # to the new in-view positions so it follows the same cells.
    _msg = _vp_channels(w2._spec)
    _msg["select"] = _sel_positions(vp.get("sel"), orig)
    widget.send(_msg)


def _vp_channels(spec):
    """The point channels to ship for an in-place plot.draw update."""
    return {"type": "vp_update",
            "x": spec.get("x"), "y": spec.get("y"), "z": spec.get("z"),
            "w": spec.get("w"), "group_data": spec.get("group_data"),
            "n_points": spec.get("n_points")}


def _sel_positions(sel, draw_order):
    """Map a logical (original-index) selection to positions in a draw order, so a
    persisted lasso re-highlights the SAME cells after a viewport swap."""
    if not sel or draw_order is None:
        return []
    inv = {int(o): p for p, o in enumerate(np.asarray(draw_order))}
    return [inv[o] for o in sel if o in inv]


def _build_grid_index(fx, fy, g=512):
    """One-time 2D grid index over the embedding so each detail-on-zoom viewport
    query is ~O(cells-in-view + hits) instead of an O(n) scan of every cell — the
    key speedup for 10M+ atlases. Points are bucketed into a g×g grid; indices are
    sorted by cell so a viewport's cell-rectangle maps to contiguous slices."""
    fx = np.asarray(fx, "float64"); fy = np.asarray(fy, "float64")
    x_min, x_max = float(np.nanmin(fx)), float(np.nanmax(fx))
    y_min, y_max = float(np.nanmin(fy)), float(np.nanmax(fy))
    xs = (x_max - x_min) or 1.0
    ys = (y_max - y_min) or 1.0
    gx = np.clip(((fx - x_min) / xs * g).astype(np.int64), 0, g - 1)
    gy = np.clip(((fy - y_min) / ys * g).astype(np.int64), 0, g - 1)
    cell = gx * g + gy
    order = np.argsort(cell, kind="stable")
    return {"g": g, "x_min": x_min, "y_min": y_min, "xs": xs, "ys": ys,
            "order": order, "cell_sorted": cell[order]}


def _grid_query(idx, x0, y0, x1, y1):
    """Candidate point indices whose grid cells overlap [x0,x1]×[y0,y1] (a superset
    of the exact box — caller does the precise bounds filter on this small set).
    Returns None for a WIDE box (spans most of the grid) so the caller uses a plain
    vectorised full scan, which is faster than looping ~every grid row."""
    g = idx["g"]; cs = idx["cell_sorted"]; order = idx["order"]
    def _gi(v, lo, span):
        return int(np.clip((v - lo) / span * g, 0, g - 1))
    gx0 = _gi(x0, idx["x_min"], idx["xs"]); gx1 = _gi(x1, idx["x_min"], idx["xs"])
    gy0 = _gi(y0, idx["y_min"], idx["ys"]); gy1 = _gi(y1, idx["y_min"], idx["ys"])
    if (gx1 - gx0 + 1) > g // 8 or (gy1 - gy0 + 1) > g // 8:
        return None                              # wide -> full scan is cheaper
    parts = []
    for gx in range(gx0, gx1 + 1):           # each grid row is a contiguous slice
        a = np.searchsorted(cs, gx * g + gy0, "left")
        b = np.searchsorted(cs, gx * g + gy1 + 1, "left")
        if b > a:
            parts.append(order[a:b])
    if not parts:
        return np.empty(0, dtype=np.int64)
    return np.concatenate(parts)


def _viewport_handler(widget, content, buffers):
    try:
        if isinstance(content, dict) and content.get("type") == "viewport":
            _viewport_payload(widget, content["bounds"])
    except Exception:
        pass


def _maybe_save(obj, save):
    """scanpy-style ``save=``: write the plot to a file. ``.html`` is supported
    programmatically (self-contained); image formats need the in-plot button."""
    if not save:
        return
    import pathlib

    ext = pathlib.Path(str(save)).suffix.lower()
    if ext in (".html", ".htm"):
        from ._export import save_html      # handles single plots AND compose grids

        save_html(obj, save)
    else:
        raise ValueError(
            f"save={save!r}: only '.html' is supported programmatically right now "
            "(a self-contained, offline file). For PNG/SVG/PDF, show the plot with "
            "enable_download=True and use its download button."
        )


def _density_sketch(x, y, target, seed):
    """Density-preserving subsample of indices: bin the 2D embedding into a grid
    and keep an even number per occupied cell, so dense blobs are thinned but
    sparse / rare regions survive (atlas-safe — unlike uniform random sampling,
    which would drop rare cell types). Vectorised, ~O(n log n). None if target>=n.
    """
    n = int(x.shape[0])
    if target >= n:
        return None
    rng = np.random.RandomState(0 if seed is None else int(seed))
    g = max(1, int(np.sqrt(target)))   # grid side scaled to the target budget

    def _bin(a):
        a = np.asarray(a, "float64")
        lo, hi = np.nanmin(a), np.nanmax(a)
        if not np.isfinite(lo) or hi <= lo:
            return np.zeros(n, dtype=np.int64)
        return np.clip(((a - lo) / (hi - lo) * g).astype(np.int64), 0, g - 1)

    cell = _bin(x) * (g + 1) + _bin(y)
    # random order within each cell: sort by cell, break ties randomly
    order = np.argsort(cell + rng.random(n), kind="stable")
    cs = cell[order]
    change = np.empty(n, dtype=bool); change[0] = True; change[1:] = cs[1:] != cs[:-1]
    run_start = np.maximum.accumulate(np.where(change, np.arange(n), 0))
    rank = np.arange(n) - run_start                     # 0-based position within its cell
    per = max(1, target // max(1, int(change.sum())))   # even budget per occupied cell
    idx = order[rank < per]
    if idx.shape[0] > target:                            # trim overshoot
        idx = rng.choice(idx, target, replace=False)
    elif idx.shape[0] < target:                          # top up to fill the budget
        mask = np.ones(n, dtype=bool); mask[idx] = False
        rest = np.where(mask)[0]
        if rest.shape[0]:
            add = rng.choice(rest, min(target - idx.shape[0], rest.shape[0]), replace=False)
            idx = np.concatenate([idx, add])
    return np.sort(idx)


def _render_index(color, n, sort_order, random_state, max_points,
                  x=None, y=None, subsample="density"):
    """Indices (into the original n points) to actually render, in draw order.

    When ``max_points`` is set and ``n`` is larger, **subsample** to that many
    points so huge datasets stay interactive (every pan re-renders only the
    subset, and points can be round again instead of perf-mode squares). Then
    apply z-depth ordering: ``random_state`` shuffles (seeded), else
    ``sort_order`` draws higher *continuous* values on top (scanpy default).
    Returns ``None`` to render all points in natural order. Indices are in
    ORIGINAL coordinates, so ``w.selection`` still maps to rows of the source.
    """
    idx = None
    if max_points is not None and n > int(max_points):
        if subsample == "density" and x is not None and y is not None:
            idx = _density_sketch(x, y, int(max_points), random_state)  # atlas-safe
        if idx is None:                                                  # random fallback
            rng = np.random.RandomState(0 if random_state is None else int(random_state))
            idx = np.sort(rng.choice(n, int(max_points), replace=False))
    m = n if idx is None else int(idx.shape[0])
    arr = np.asarray(color) if color is not None else None
    order = None
    if random_state is not None:
        order = np.random.RandomState(int(random_state)).permutation(m)
    elif sort_order and arr is not None and np.issubdtype(arr.dtype, np.number):
        sub = arr if idx is None else arr[idx]
        order = np.argsort(sub, kind="stable")   # ascending -> high on top
    if order is not None:
        idx = order if idx is None else idx[order]
    return idx


def _resolve_numeric(spec, data, layer=None, param="size_by", use_raw=None):
    """Resolve a size/opacity encoding to a numeric array: a column name (in a
    DataFrame or AnnData ``.obs``), a ``var_names`` feature, or a raw vector.
    Returns None if not set."""
    if spec is None:
        return None
    if isinstance(spec, str):
        if hasattr(data, "columns") and spec in getattr(data, "columns", []):
            return pd.to_numeric(data[spec]).to_numpy()
        obs = getattr(data, "obs", None)
        if obs is not None and spec in getattr(obs, "columns", []):
            return pd.to_numeric(obs[spec]).to_numpy()
        # feature (gene) fallback, reusing the color resolver (re-label its
        # error with the actual parameter name)
        try:
            if _is_anndata(data):
                vec, _ = _resolve_anndata_vec(data, spec, layer, use_raw)
                return np.asarray(vec, dtype="float64")
            if _is_mudata(data) and ":" in spec:
                mod, key = spec.split(":", 1)
                if mod in data.mod:
                    vec, _ = _resolve_anndata_vec(data.mod[mod], key, layer, use_raw)
                    return np.asarray(vec, dtype="float64")
        except KeyError:
            pass
        avail = list(getattr(data, "columns", []) or [])
        if obs is not None:
            avail += list(getattr(obs, "columns", []))
        if hasattr(data, "var_names"):
            avail += list(getattr(data, "var_names", []))
        raise _not_found(param, spec, avail, searched="columns / .obs / .var_names")
    return np.asarray(spec, dtype="float64")


def _resolve_cols(spec, data):
    """Resolve tooltip fields to a dict {name: array}: column name(s) (DataFrame
    or AnnData obs / var feature), a dict, or a DataFrame. None if not set."""
    if spec is None:
        return None
    import numpy as np
    if isinstance(spec, dict):
        return spec
    if hasattr(spec, "columns") and hasattr(spec, "iloc"):  # DataFrame
        return {c: spec[c].to_numpy() for c in spec.columns}
    names = [spec] if isinstance(spec, str) else list(spec)
    out = {}
    for nm in names:
        if hasattr(data, "columns") and nm in getattr(data, "columns", []):
            out[nm] = data[nm].to_numpy()
        elif hasattr(data, "obs") and nm in getattr(data.obs, "columns", []):
            out[nm] = data.obs[nm].to_numpy()
        elif hasattr(data, "var_names") and nm in list(getattr(data, "var_names", [])):
            col = data[:, nm].X
            out[nm] = col.toarray().ravel() if hasattr(col, "toarray") else np.asarray(col).ravel()
        else:
            avail = list(getattr(data, "columns", []) or [])
            if hasattr(data, "obs"):
                avail += list(getattr(data.obs, "columns", []))
            if hasattr(data, "var_names"):
                avail += list(getattr(data, "var_names", []))
            raise _not_found("tooltip_by", nm, avail,
                             searched="columns / .obs / .var_names")
    return out


def scatterplot(
    data: Any = None,
    *,
    # --- scanpy-style names (preferred) -------------------------------------
    color: Any = _UNSET,            # alias of color_by (name, list of names, or vector)
    size: Any = _UNSET,             # scalar -> point_size; name/array -> size_by
    cmap: Any = _UNSET,             # alias of continuous_palette
    palette: Any = _UNSET,          # alias of categorical_palette
    components: Any = _UNSET,       # 1-based embedding dims, e.g. (1, 2); alias of dims
    ncols: Optional[int] = None,    # grid columns when color is a list
    save: Optional[str] = None,     # write to a file (.html now; see docs)
    na_color: str = "lightgray",    # colour for NaN / un-selected categories
    groups: Any = None,             # show only these categories; grey the rest
    sort_order: bool = True,        # draw higher continuous values on top (scanpy)
    random_state: Optional[int] = None,  # seed -> random draw order (reduce overplotting)
    max_points: Any = "auto",            # subsample huge data for smooth exploration; "auto" caps at 500k, None = all points
    subsample: str = "density",           # "density" (atlas-safe: keeps rare cells) or "random"
    # --- original names (still supported as aliases) ------------------------
    basis: Optional[Union[str, int]] = None,
    x: Optional[Union[str, int]] = None,
    y: Optional[Union[str, int]] = None,
    color_by: ColorSpec = None,
    group_by: ColorSpec = None,
    layer: Optional[str] = None,
    use_raw: Optional[bool] = None,
    dims: Optional[tuple] = None,
    table: Optional[str] = None,
    point_size: Optional[float] = None,
    opacity: Optional[float] = None,
    point_color: Optional[str] = None,
    size_by: Any = None,
    opacity_by: Any = None,
    tooltip_by: Any = None,
    pixel_ratio: Optional[float] = None,
    categorical_palette: str = "Set1",
    continuous_palette: str = "viridis",
    custom_palette: Optional[Sequence[str]] = None,
    custom_colors: Optional[dict] = None,
    vmin: Any = None,
    vmax: Any = None,
    center_zero: bool = False,
    title: Optional[str] = None,
    xlab: Optional[str] = None,
    ylab: Optional[str] = None,
    legend_title: Optional[str] = None,
    show_axes: bool = True,
    show_tooltip: bool = True,
    background_color: Optional[str] = None,
    axis_color: str = "#333333",
    legend_bg: str = "#ffffff",
    legend_text: str = "#000000",
    legend_opacity: Optional[float] = None,
    legend_blur: Optional[float] = None,
    toolbar: str = "left",
    zoom_on_selection: bool = False,
    legend_position: str = "top-right",
    draggable_legend: bool = True,
    enable_download: bool = False,
    font_size: int = 12,
    legend_font_size: int = 12,
    auto_fit: bool = False,
    range_padding: float = 0.15,
    xrange: Optional[Sequence[float]] = None,
    yrange: Optional[Sequence[float]] = None,
    filter_by: Any = None,
    point_labels: Optional[Sequence] = None,
    plot_id: Optional[str] = None,
    width: Optional[int] = 700,
    height: int = 500,
    backend: str = "regl",
    interactive: bool = False,
    fast: bool = False,             # experimental: binary transfer (implies interactive)
    progressive: bool = False,      # experimental: density-sketch overview + full detail as you zoom in
    detail_on_zoom: bool = False,   # internal: emit viewport messages so the kernel re-renders in-view cells
    performance_mode: Optional[bool] = None,  # squares+no-blend (faster) vs round circles; None=auto (n>500k), forced on for progressive unless set False
    _color_categories=None,         # internal: pin the full categorical level set (detail-on-zoom colour consistency)
    show: bool = True,
    **backend_kwargs: Any,
):
    """Interactive WebGL scatterplot from single-cell / tabular data.

    Parameters
    ----------
    data
        ``AnnData``, ``MuData``, ``SpatialData``, pandas ``DataFrame`` or numpy
        array. See :func:`reglscatterpy._extract.extract`.
    basis
        Embedding selector for single-cell objects (an ``obsm`` key, e.g.
        ``"umap"``/``"X_umap"``/``"pca"``; short names are auto-prefixed and
        matched case-insensitively). Preferred over ``x`` for embeddings.
    x, y
        For single-cell objects, ``x`` is an alias for ``basis`` (same meaning).
        For DataFrame / array inputs, ``x``/``y`` are the column names / indices
        of the coordinates. ``basis`` is not valid for tables/arrays.
    color_by
        Pass a single name to colour by one ``obs`` column / feature, or a
        **list of names** (genes and/or ``obs`` columns) to render a linked
        multi-panel grid — one panel per value (a raw per-point colour vector
        must be a numpy array / pandas Series, not a list of strings).
    color_by, group_by
        ``obs`` column, a feature in ``var_names``, or a raw vector
        (``"modality:feature"`` for ``MuData``).
    point_color
        A fixed hex colour overriding ``color_by``.
    pixel_ratio
        Device-pixel-ratio for the WebGL backing store (defaults to a crisp
        ``max(devicePixelRatio, 2)`` in-widget).
    categorical_palette
        A ColorBrewer qualitative palette name (``"Set1"``, ``"Dark2"`` ...).
    continuous_palette
        A viridis-family name (``"viridis"``, ``"magma"`` ...).
    custom_colors, custom_palette
        Explicit per-level mapping (dict) or ordered list of colours.
    filter_by
        A dict / DataFrame of numeric columns shown as interactive range
        filters.
    backend
        ``"regl"`` (default) renders the package's own widget; ``"jscatter"``
        uses jupyter-scatter.
    show
        When ``True`` (default) returns the displayable widget.

    Returns
    -------
    The widget object (displays inline in notebooks).

    Notes
    -----
    Argument names follow scanpy where possible — ``color``, ``size``, ``cmap``,
    ``palette``, ``components`` (1-based), ``ncols``, ``layer``, ``vmin``/``vmax``,
    ``save`` — and the original names (``color_by``, ``point_size``,
    ``continuous_palette``, ``categorical_palette``, ``dims`` …) keep working.
    """
    # snapshot the call args up front so progressive mode can re-render the FULL
    # data through this exact same pipeline (keeps lasso/filter/tooltip identical)
    _params = {k: v for k, v in locals().items() if k != "data"}

    # --- scanpy-style aliases win over the original names when given --------
    if color is not _UNSET:
        color_by = color
    if size is not _UNSET:
        if isinstance(size, bool):
            raise TypeError("size must be a number or a column/feature name.")
        if isinstance(size, (int, float)):
            point_size = size           # scalar -> global point size
        else:
            size_by = size              # name / array -> per-point size
    if cmap is not _UNSET:
        continuous_palette = cmap
    if palette is not _UNSET:
        categorical_palette = palette
    if components is not _UNSET and components is not None:
        dims = tuple(int(c) - 1 for c in components)   # 1-based -> 0-based
    _auto_max = (max_points == "auto")
    if _auto_max:
        max_points = _DEFAULT_MAX_POINTS   # smooth by default; pass max_points=None for all points
    if fast:
        interactive = True   # binary transfer rides the live comm; needs a widget
    if pixel_ratio is not None and pixel_ratio < 1:
        import warnings
        warnings.warn(
            f"pixel_ratio={pixel_ratio} < 1 downsamples (blurry) and can break "
            "rendering; clamped to 1. pixel_ratio supersamples for crispness — use "
            ">= 1 (e.g. 1.5 or 2).",
            stacklevel=2,
        )
        pixel_ratio = 1.0

    # --- progressive: subset now (instant), full streamed in the background ---
    if progressive and not _is_name_list(color_by) and backend == "regl":
        import uuid
        interactive = True
        # detail-on-zoom ("in-memory tiling"): show a density-sketch overview, then
        # re-render ALL cells inside the viewport as you zoom in (the view holds few
        # cells when zoomed, so we can draw them all -> full detail + complete lasso,
        # with NO preprocessing). Shared plotId so each re-render replaces the plot.
        budget = max_points if (isinstance(max_points, int) and max_points) else _DEFAULT_MAX_POINTS
        # A lighter overview re-renders fast on zoom-out (the heaviest refresh);
        # zoom-in detail still uses the full budget.
        overview_budget = min(budget, 250_000)
        pid = plot_id or ("rs_" + uuid.uuid4().hex[:10])
        bk = dict(_params.pop("backend_kwargs", {}) or {})
        # base64 channels (fast=False) so in-view updates can ship over a custom
        # message and be applied via plot.draw (no _spec change -> no re-render).
        base = {**_params, "progressive": False, "interactive": True, "show": False,
                "fast": False, "plot_id": pid, "detail_on_zoom": True}
        w = scatterplot(data, **{**base, "max_points": overview_budget}, **bk)   # overview
        try:
            _io = _is_anndata(data) or _is_mudata(data) or _is_spatialdata(data)
            _eff = basis if (_io and basis is not None) else x
            _full = extract(data, x=_eff, y=y, layer=layer, use_raw=use_raw,
                            dims=dims, table=table)
            _lg = w._spec.get("legend") or {}
            _fx = np.asarray(_full.x, "float64")
            _fy = np.asarray(_full.y, "float64")
            w._vp = {"data": data, "base": base, "bk": bk, "budget": budget,
                     "x": _fx, "y": _fy,
                     "grid": None,   # built in a background thread (below) so it
                                     # doesn't delay first paint; full-scan until ready

                     # keep the overview's data domain so the camera stays put and
                     # the in-view points land in the right place.
                     "xrange": (w._spec.get("x_min"), w._spec.get("x_max")),
                     "yrange": (w._spec.get("y_min"), w._spec.get("y_max")),
                     # keep point size + colour scale fixed across refreshes.
                     "point_size": (w._spec.get("options") or {}).get("size"),
                     "vmin": _lg.get("minVal"), "vmax": _lg.get("maxVal"),
                     # full categorical level set -> detail views keep identical
                     # colours + codes (else zoomed-in colours shift & legend mis-filters).
                     "categories": (_lg.get("names")
                                    if _lg.get("var_type") == "categorical" else None),
                     # the JS caches the overview buffers; zoom-out just asks it to
                     # redraw them. Track state to skip redundant redraws on pan.
                     "showing_overview": True,
                     # logical lasso (ORIGINAL indices) so it persists across zoom;
                     # updated only by real user lassos (re-applies use msg.select,
                     # which never touches the _selection trait -> no observer fire).
                     "sel": set(),
                     "overview_draw_order": w._draw_order}
            w.on_msg(_viewport_handler)

            def _on_user_sel(change, _w=w):
                vp = getattr(_w, "_vp", None)
                if vp is None or vp.get("_setting"):
                    return
                try:
                    do = _w._draw_order
                    pos = change.get("new") or []
                    if do is not None:
                        od = np.asarray(do)
                        vp["sel"] = {int(od[p]) for p in pos if 0 <= p < od.size}
                    else:
                        vp["sel"] = {int(p) for p in pos}
                except Exception:
                    pass
            w.observe(_on_user_sel, names="_selection")

            def _mk_grid(_w=w, _fx=_fx, _fy=_fy):
                try:
                    _w._vp["grid"] = _build_grid_index(_fx, _fy)
                except Exception:
                    pass
            import threading
            threading.Thread(target=_mk_grid, daemon=True).start()
        except Exception:
            pass
        return w

    # --- validate enum-ish arguments up front (fail before doing work) ------
    if backend not in ("regl", "jscatter"):
        raise ValueError(
            f"backend={backend!r} is invalid; choose 'regl' or 'jscatter'."
        )
    if toolbar not in _TOOLBAR_CHOICES:
        raise ValueError(
            f"toolbar={toolbar!r} is invalid; choose one of "
            f"{[c for c in _TOOLBAR_CHOICES if c is not None]} or None."
        )

    # A single-element list is just one plot (scanpy's color=['gene']) — unwrap
    # so it renders as a normal single plot (700px), not a full-width 1-up grid.
    if _is_name_list(color_by) and len(color_by) == 1:
        color_by = color_by[0]

    # --- color_by as a list of names -> one linked panel per value ----------
    if _is_name_list(color_by):
        if backend != "regl":
            raise ValueError(
                "color_by as a list of names (multi-panel grid) requires "
                "backend='regl'."
            )
        from ._compose import compose

        panels = [
            scatterplot(
                data, basis=basis, x=x, y=y, color_by=name, group_by=group_by,
                layer=layer, use_raw=use_raw, dims=dims, table=table, point_size=point_size,
                opacity=opacity, point_color=point_color, size_by=size_by,
                opacity_by=opacity_by, tooltip_by=tooltip_by,
                pixel_ratio=pixel_ratio, categorical_palette=categorical_palette,
                continuous_palette=continuous_palette, custom_palette=custom_palette,
                custom_colors=custom_colors, vmin=vmin, vmax=vmax,
                center_zero=center_zero, na_color=na_color, groups=groups,
                sort_order=sort_order, random_state=random_state, max_points=max_points,
                subsample=subsample,
                title=(title or name), xlab=xlab, ylab=ylab,
                legend_title=legend_title, show_axes=show_axes,
                show_tooltip=show_tooltip, background_color=background_color,
                axis_color=axis_color, legend_bg=legend_bg, legend_text=legend_text,
                legend_opacity=legend_opacity, legend_blur=legend_blur,
                toolbar=toolbar, zoom_on_selection=zoom_on_selection,
                legend_position=legend_position, draggable_legend=draggable_legend,
                enable_download=enable_download, font_size=font_size,
                legend_font_size=legend_font_size, auto_fit=auto_fit,
                range_padding=range_padding, xrange=xrange, yrange=yrange,
                filter_by=filter_by, point_labels=point_labels, plot_id=None,
                width=None, height=height, backend=backend,  # responsive in the grid
                interactive=True, show=False,  # linked sync needs live widgets
                **backend_kwargs,
            )
            for name in color_by
        ]
        grid = compose(panels, cols=ncols)
        _maybe_save(grid, save)
        return grid

    # --- resolve the effective embedding (basis is preferred; x is alias) ---
    is_object = _is_anndata(data) or _is_mudata(data) or _is_spatialdata(data)
    if basis is not None and not is_object:
        raise ValueError(
            "basis= applies to AnnData/MuData/SpatialData inputs; for a "
            "DataFrame/array use x=/y= column selectors."
        )
    eff_x = basis if (is_object and basis is not None) else x

    pd_data: PlotData = extract(
        data, x=eff_x, y=y, color_by=color_by, group_by=group_by,
        layer=layer, use_raw=use_raw, dims=dims, table=table,
    )

    # adaptive defaults, matching the R heuristics
    n = pd_data.n
    if point_size is None:
        point_size = 1 if n > 500_000 else 5 if n < 5_000 else 4 if n < 50_000 else 3
    if opacity is None:
        opacity = 1.0 if n > 500_000 else 0.8

    if backend == "jscatter":
        return _scatterplot_jscatter(
            pd_data, point_size=point_size, opacity=opacity,
            categorical_palette=custom_colors or custom_palette or categorical_palette,
            continuous_palette=continuous_palette, title=title, xlab=xlab, ylab=ylab,
            show_axes=show_axes, show_tooltip=show_tooltip,
            background_color=background_color, width=width, height=height,
            show=show, **backend_kwargs,
        )
    # palette names are validated on the regl path (jscatter uses its own)
    if categorical_palette not in QUALITATIVE:
        raise _not_found("categorical_palette", categorical_palette,
                         list(QUALITATIVE), searched="built-in qualitative palettes")
    if continuous_palette not in CONTINUOUS:
        raise _not_found("continuous_palette", continuous_palette,
                         list(CONTINUOUS), searched="built-in continuous palettes")

    size_values = _resolve_numeric(size_by, data, layer, "size_by", use_raw)
    opacity_values = _resolve_numeric(opacity_by, data, layer, "opacity_by", use_raw)
    tooltip_cols = _resolve_cols(tooltip_by, data)
    for pname, vec in (("size_by", size_values), ("opacity_by", opacity_values)):
        if vec is not None and len(vec) != n:
            raise ValueError(
                f"{pname} has length {len(vec)} but the data has {n} points."
            )

    # Render index = which points to draw (subsample for huge data) + draw order
    # (z-depth). Subset/reorder ALL per-point arrays consistently; the index is
    # stored on the widget and w.selection is translated back to ORIGINAL data
    # indices at the Python boundary, so the JS never needs to know.
    draw_order = _render_index(pd_data.color, n, sort_order, random_state, max_points,
                               x=pd_data.x, y=pd_data.y, subsample=subsample)
    if draw_order is not None:
        import dataclasses

        pd_data = dataclasses.replace(
            pd_data,
            x=np.asarray(pd_data.x)[draw_order], y=np.asarray(pd_data.y)[draw_order],
            color=(np.asarray(pd_data.color)[draw_order] if pd_data.color is not None else None),
            group=(np.asarray(pd_data.group)[draw_order] if pd_data.group is not None else None),
        )
        if size_values is not None:
            size_values = np.asarray(size_values)[draw_order]
        if opacity_values is not None:
            opacity_values = np.asarray(opacity_values)[draw_order]
        if tooltip_cols:
            tooltip_cols = {k: np.asarray(v)[draw_order] for k, v in tooltip_cols.items()}
        if filter_by is not None:
            filter_by = pd.DataFrame(filter_by).iloc[draw_order].reset_index(drop=True)
        if point_labels is not None:
            point_labels = [point_labels[int(i)] for i in draw_order]

    spec = build_payload(
        pd_data,
        point_size=point_size, opacity=opacity, point_color=point_color,
        size_values=size_values, opacity_values=opacity_values,
        tooltip_cols=tooltip_cols,
        pixel_ratio=pixel_ratio,
        categorical_palette=categorical_palette, continuous_palette=continuous_palette,
        custom_palette=custom_palette, custom_colors=custom_colors,
        vmin=vmin, vmax=vmax, center_zero=center_zero,
        na_color=na_color, groups=groups, categories=_color_categories,
        xrange=xrange, yrange=yrange, range_padding=range_padding,
        xlab=xlab, ylab=ylab, title=title, legend_title=legend_title,
        show_axes=show_axes, show_tooltip=show_tooltip,
        background_color=background_color, axis_color=axis_color,
        legend_bg=legend_bg, legend_text=legend_text,
        legend_opacity=legend_opacity, legend_blur=legend_blur,
        toolbar_position=toolbar, zoom_on_selection=zoom_on_selection,
        legend_position=legend_position, draggable_legend=draggable_legend,
        enable_download=enable_download, font_size=font_size,
        legend_font_size=legend_font_size, auto_fit=auto_fit,
        point_labels=point_labels, plot_id=plot_id, filter_by=filter_by,
        binary=fast,
    )
    if detail_on_zoom:
        spec["detailOnZoom"] = True   # client emits viewport msgs -> kernel re-renders in-view cells
        # squares + no alpha-blend are much cheaper per draw; default on for the
        # zoom loop, but honour an explicit performance_mode=False (round circles).
        spec["performanceMode"] = True if performance_mode is None else bool(performance_mode)
    elif performance_mode is not None:
        spec["performanceMode"] = bool(performance_mode)

    # Be honest about subsampling: caption the plot ("X of Y shown") and, when the
    # downsample was automatic (not user-requested), warn — so a subsampled plot is
    # never mistaken for the full data.
    _rendered = n if draw_order is None else int(np.asarray(draw_order).size)
    if _rendered < n:
        spec["caption"] = f"{_rendered:,} of {n:,} shown"
        if _auto_max and not detail_on_zoom:
            import warnings
            warnings.warn(
                f"Showing a {_rendered:,}-point density-preserving subsample of "
                f"{n:,} points for smooth rendering (rare groups are kept). "
                f"Pass max_points=None for all points, or set max_points=N.",
                stacklevel=2,
            )

    w = int(width) if width else 0   # 0 => responsive (100%)
    if interactive:
        # Live, kernel-linked widget: w.selection round-trips (needs a kernel).
        from ._widget import ReglScatter

        widget = ReglScatter()
        widget._height = int(height)
        widget._width = w
        widget._source = data   # so w.annotate(...) can write back to obs/colData
        widget._draw_order = draw_order   # selection translates through this
        widget._spec = spec
        _maybe_save(widget, save)
        return widget

    # Default: a static, self-contained plot — renders with no comm and reopens
    # with no kernel (like a plotly figure), and is NOT an ipywidget, so nothing
    # is written to the notebook's widget-state.
    from ._widget import StaticPlot

    plot = StaticPlot(spec=spec, source=data, height=int(height), width=w)
    plot._draw_order = draw_order
    _maybe_save(plot, save)
    return plot


def _scatterplot_jscatter(
    pd_data: PlotData, *, point_size, opacity, categorical_palette,
    continuous_palette, title, xlab, ylab, show_axes, show_tooltip,
    background_color, width, height, show, **jscatter_kwargs,
):
    try:
        import jscatter
    except ModuleNotFoundError as exc:  # pragma: no cover - import guard
        raise ModuleNotFoundError(
            "backend='jscatter' needs 'jupyter-scatter'. "
            "Install with: pip install reglscatterpy[render]"
        ) from exc

    df = pd.DataFrame({"x": pd_data.x, "y": pd_data.y})
    color_col = None
    if pd_data.color is not None:
        color_col = pd_data.color_name or "color"
        series = pd.Series(pd_data.color)
        if not pd.api.types.is_numeric_dtype(series):
            series = series.astype("category")
        df[color_col] = series.to_numpy()

    opts: dict[str, Any] = {
        "data": df, "x": "x", "y": "y", "size": point_size, "opacity": opacity,
        "height": height, "axes": show_axes,
    }
    if color_col is not None:
        opts["color_by"] = color_col
        if not pd.api.types.is_numeric_dtype(df[color_col]):
            if categorical_palette is not None:
                opts["color_map"] = categorical_palette
        else:
            opts["color_map"] = continuous_palette
    if background_color is not None:
        opts["background_color"] = background_color
    if width is not None:
        opts["width"] = width
    if title is not None:
        opts["title"] = title
    opts.update(jscatter_kwargs)

    scatter = jscatter.plot(**opts)
    try:
        scatter.label(x=xlab or pd_data.xlab, y=ylab or pd_data.ylab)
    except Exception:  # pragma: no cover
        pass
    if not show_tooltip:
        try:
            scatter.tooltip(False)
        except Exception:  # pragma: no cover
            pass
    return scatter.show() if show else scatter
