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

import logging
from typing import Any, Optional, Sequence, Union

import numpy as np
import pandas as pd

_log = logging.getLogger("reglscatterpy")

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


def _build_detail(vp, idx, max_points):
    """Build a detail widget for the given ORIGINAL indices by numpy-indexing the
    cached full channels — no AnnData re-slice / re-extract / re-resolution per
    viewport (the dominant cost on big sparse atlases). build_payload still does
    the channel packing + colour encoding, so the result is byte-identical to a
    full rebuild (same colours/codes — no shift on zoom)."""
    import dataclasses
    fpd = vp["full_pd"]
    idx = np.asarray(idx)
    sub_pd = dataclasses.replace(
        fpd,
        x=np.asarray(fpd.x)[idx], y=np.asarray(fpd.y)[idx],
        color=(np.asarray(fpd.color)[idx] if fpd.color is not None else None),
        group=(np.asarray(fpd.group)[idx] if fpd.group is not None else None),
    )
    fs, fo, ft, ff = (vp["full_size"], vp["full_opacity"],
                      vp["full_tooltip"], vp["full_filter"])
    sub_sz = None if fs is None else np.asarray(fs)[idx]
    sub_op = None if fo is None else np.asarray(fo)[idx]
    sub_tt = None if not ft else {k: np.asarray(v)[idx] for k, v in ft.items()}
    sub_ft = None if ff is None else ff.iloc[idx].reset_index(drop=True)
    extra = {"max_points": max_points}
    xr = vp.get("xrange")
    if xr and xr[0] is not None:                   # keep the overview domain
        extra["xrange"] = vp["xrange"]; extra["yrange"] = vp["yrange"]
    if vp.get("point_size") is not None:           # keep size constant
        extra["point_size"] = vp["point_size"]
    if vp.get("vmin") is not None:                 # keep the colour scale constant
        extra["vmin"] = vp["vmin"]; extra["vmax"] = vp["vmax"]
    if vp.get("categories") is not None:           # keep categorical colours/codes constant
        extra["_color_categories"] = vp["categories"]
    return scatterplot(None, **{**vp["base"], **extra,
                                "_pd_data": sub_pd, "_size_values": sub_sz,
                                "_opacity_values": sub_op, "_tooltip_cols": sub_tt,
                                "_filter_df": sub_ft}, **vp["bk"])


def _viewport_payload(widget, bounds, pad=0.6, seq=None):
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
    pad = vp.get("pad", pad)
    x0, y0, x1, y1 = (float(b) for b in bounds)
    xr, yr = vp.get("xrange"), vp.get("yrange")
    # Zoomed back out to ~the full extent? Push the cached overview channels
    # instantly (no recompute, no rebuild) — and via plot.draw, not a re-render.
    if (xr and xr[0] is not None
            and (x1 - x0) >= 0.8 * (xr[1] - xr[0])
            and (y1 - y0) >= 0.8 * (yr[1] - yr[0])):
        if vp.get("showing_overview"):
            widget.send({"type": "vp_noop", "seq": seq})   # already at overview; clear loader
            return
        vp["showing_overview"] = True
        vp["last_fetch"] = None                 # so the next zoom-in re-fetches detail
        sel = vp.get("sel")
        if not sel:
            widget._inv_draw_order = None
            widget._draw_order = vp["overview_draw_order"]
            widget._source = vp["data"]
            widget.send({"type": "vp_overview", "seq": seq,   # cached overview (fast)
                         "select": [], "keep": _keep_positions(vp, vp["overview_draw_order"])})
            return
        # selection present -> render the overview density-sketch UNION the selected
        # cells so they're all drawn + highlighted (else only the ~1% in the sketch show).
        ov = np.asarray(vp["overview_draw_order"])
        aug = np.unique(np.concatenate([ov, np.fromiter(sel, dtype=np.int64)]))
        data = vp["data"]
        w2 = _build_detail(vp, aug, None)       # keep ALL selected (don't re-subsample)
        do = w2._draw_order
        orig = aug if do is None else aug[np.asarray(do)]
        widget._inv_draw_order = None
        widget._draw_order = orig
        widget._source = data
        _vp_send(widget, w2._spec, _sel_positions(sel, orig), seq, _keep_positions(vp, orig))
        return
    vp["showing_overview"] = False
    # Skip the round-trip if the new view is still inside the region we last fetched
    # (a pan within the overscan margin) and we haven't zoomed in enough to need denser
    # detail — the points are already on screen, so re-fetching just causes flicker/lag.
    _vspan = x1 - x0
    lf = vp.get("last_fetch")
    if lf is not None:
        px0, py0, px1, py1, lspan = lf
        if (x0 >= px0 and x1 <= px1 and y0 >= py0 and y1 <= py1
                and _vspan >= 0.85 * lspan):
            widget.send({"type": "vp_noop", "seq": seq})   # covered; clear the loader
            return
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
        widget.send({"type": "vp_noop", "seq": seq})       # nothing in view; clear loader
        return
    # Dense region: a tight zoom into a packed cluster can still hold millions of
    # cells. Pre-cap the candidates (cheap random) BEFORE building so subset+density
    # -sketch stay fast; the density-sketch below still refines to the budget.
    budget = int(vp["budget"])
    if in_idx.size > budget * 3:
        _r = np.random.RandomState(0)
        in_idx = np.sort(_r.choice(in_idx, budget * 3, replace=False))
    data = vp["data"]
    w2 = _build_detail(vp, in_idx, vp["budget"])
    do = w2._draw_order
    orig = in_idx if do is None else in_idx[np.asarray(do)]
    widget._inv_draw_order = None
    widget._draw_order = orig
    widget._source = data                          # keep selection/DE on the full object
    # Swap points on the EXISTING plot (custom message -> plot.draw); never touches
    # _spec, so no re-mount / spinner / plot recreation. Re-map the persisted lasso
    # to the new in-view positions so it follows the same cells.
    vp["last_fetch"] = (x0, y0, x1, y1, _vspan)    # padded region + the view span we fetched at
    _vp_send(widget, w2._spec, _sel_positions(vp.get("sel"), orig), seq, _keep_positions(vp, orig))


def _vp_send(widget, spec, select, seq, keep=None):
    """Ship an in-place plot.draw update. The big per-point channels (x/y/z/w) go
    as raw binary comm BUFFERS (no base64 inflate/decode); the small bits
    (group/tooltip/filter) stay in the JSON content."""
    order, buffers = [], []
    for k in ("x", "y", "z", "w"):
        ch = spec.get(k)
        if ch is not None:
            order.append(k)
            buffers.append(ch if isinstance(ch, memoryview) else memoryview(np.ascontiguousarray(ch)))
    content = {"type": "vp_update", "n_points": spec.get("n_points"),
               "buf_order": order,
               "group_data": spec.get("group_data"),
               "tooltip_data": spec.get("tooltip_data"),
               "filter_data": spec.get("filter_data"),
               "select": select, "keep": keep, "seq": seq}
    widget.send(content, buffers=buffers)


def _keep_positions(vp, draw_order):
    """Legend-filter keep positions for the current draw order (None = no filter).

    ``filter_keep`` is a numpy array of ORIGINAL cells to keep (set once from the
    filtered panel's codes — works cross-variable). Each panel maps it to its own
    displayed positions with a vectorised ``np.isin`` over the ~cells-in-view draw
    order — no Python set of millions, no per-element dict lookups (the old version
    of both was why legend filtering felt slow on 20M-cell atlases)."""
    ko = vp.get("filter_keep")
    if ko is None:
        return None
    do = vp.get("overview_draw_order") if draw_order is None else draw_order
    if do is None:
        return None
    return np.where(np.isin(np.asarray(do), ko))[0].tolist()


_INV_CACHE = {}  # id(draw_order) -> (draw_order_ref, {original_index: position})


def _inv_map(draw_order):
    """Cached original-index -> draw-position map. Rebuilding this dict on every
    viewport swap (per panel) is O(n) and was dominating selection/filter syncing
    for 10M-cell atlases; draw_order is stable per panel, so cache by identity."""
    key = id(draw_order)
    cached = _INV_CACHE.get(key)
    if cached is not None and cached[0] is draw_order:
        return cached[1]
    inv = {int(o): p for p, o in enumerate(np.asarray(draw_order))}
    if len(_INV_CACHE) > 32:        # bound the cache (a session rarely needs many)
        _INV_CACHE.clear()
    _INV_CACHE[key] = (draw_order, inv)
    return inv


def _sel_positions(sel, draw_order):
    """Map a logical (original-index) selection to positions in a draw order, so a
    persisted lasso re-highlights the SAME cells after a viewport swap."""
    if not sel or draw_order is None:
        return []
    inv = _inv_map(draw_order)
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


def _points_in_polygon(fx, fy, poly):
    """Boolean mask of points inside the polygon (data coords). matplotlib when
    available, else a vectorised ray-casting fallback."""
    poly = np.asarray(poly, "float64")
    x = np.asarray(fx, "float64"); y = np.asarray(fy, "float64")
    try:
        from matplotlib.path import Path
        return Path(poly).contains_points(np.column_stack([x, y]))
    except Exception:
        n = poly.shape[0]
        inside = np.zeros(x.shape[0], dtype=bool)
        j = n - 1
        for i in range(n):
            xi, yi = poly[i]; xj, yj = poly[j]
            cond = ((yi > y) != (yj > y)) & (
                x < (xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi)
            inside ^= cond
            j = i
        return inside


def _ack_work(widget):
    """Mark the front-end's latest HEAVY-work request (_work_req) as DONE now that
    this async handler has finished, so a synchronous w.selection/w.filtered read
    blocks for completion rather than just delivery (see _PlotAPI._pump). Light
    interactions (clicks / non-progressive trait sets) use the _sel_gen barrier,
    which the kernel acks automatically."""
    try:
        widget._work_done = int(getattr(widget, "_work_req", 0) or 0)
    except Exception:
        pass


def _lasso_payload(widget, polygon):
    """Full-region lasso: select EVERY cell inside the polygon on the full dataset
    (not just the rendered subset), then re-highlight whatever's currently in view."""
    vp = getattr(widget, "_vp", None)
    if not vp or not polygon or len(polygon) < 3:
        return
    mask = _points_in_polygon(vp["x"], vp["y"], polygon)
    sel = {int(i) for i in np.where(mask)[0]}        # original cells in the lasso
    # propagate to every linked panel (same embedding/cells) — each highlights its
    # own in-view subset of the SAME original cells.
    for panel in (vp.get("group") or [widget]):
        pvp = getattr(panel, "_vp", None)
        if pvp is None:
            continue
        pvp["sel"] = sel
        pvp["showing_overview"] = False              # rebuild overview∪selection on zoom-out
        panel._source = pvp["data"]
        panel.send({"type": "vp_select",
                    "select": _sel_positions(sel, panel._draw_order)})


def _compute_filter_keep(vp):
    """Intersect the active legend-category AND range-slider filters over the FULL
    dataset. Returns kept ORIGINAL indices (np.int64 array) or None when nothing is
    filtering (so w.filtered / the plot show every cell)."""
    mask = None
    # legend categories (clicked panel's categorical colour)
    cats = vp.get("filter_cats")
    codes, levels = vp.get("color_codes"), vp.get("color_levels")
    if cats and codes is not None and levels is not None:
        want = {str(c) for c in cats}
        sel_codes = np.array([i for i, lv in enumerate(levels) if lv in want],
                             dtype=np.int64)
        m = np.isin(np.asarray(codes), sel_codes)
        mask = m if mask is None else (mask & m)
    # range sliders over the pre-extracted FULL filter columns (vp["full_filter"])
    ranges = vp.get("filter_ranges") or {}
    ff = vp.get("full_filter")
    if ranges and ff is not None:
        cols = getattr(ff, "columns", [])
        for col, rng in ranges.items():
            if rng is None or col not in cols:
                continue
            vals = np.asarray(ff[col], dtype="float64")
            m = (vals >= float(rng[0])) & (vals <= float(rng[1]))
            mask = m if mask is None else (mask & m)
    return None if mask is None else np.where(mask)[0]


def _apply_progressive_filters(widget):
    """Recompute the combined keep on the acting panel and share it (by ORIGINAL
    cell) to every linked panel, pushing each a vp_filter to redraw its in-view
    subset. Sets vp["filter_keep"] — the source w.filtered reads."""
    vp = getattr(widget, "_vp", None)
    if not vp:
        return
    keep = _compute_filter_keep(vp)
    for panel in (vp.get("group") or [widget]):
        pvp = getattr(panel, "_vp", None)
        if pvp is None:
            continue
        pvp["filter_keep"] = keep
        try:
            panel.send({"type": "vp_filter",
                        "keep": _keep_positions(pvp, panel._draw_order)})
        except Exception:
            pass


def _legend_filter_payload(widget, cats):
    """Legend category filter (progressive): store the kept categories, then re-apply
    the combined (category ∩ range) filter over the FULL dataset, synced by ORIGINAL
    cell across linked panels (cross-variable — others just hide the same cells)."""
    vp = getattr(widget, "_vp", None)
    if not vp:
        return
    vp["filter_cats"] = list(cats) if cats else None
    _apply_progressive_filters(widget)


def _range_filter_payload(widget, ranges):
    """Range-slider filter (progressive): store the per-column [lo, hi] ranges, then
    re-apply the combined filter over the FULL dataset — so w.filtered and the plot
    reflect every cell in range, not just the drawn subset."""
    vp = getattr(widget, "_vp", None)
    if not vp:
        return
    vp["filter_ranges"] = {str(k): [float(v[0]), float(v[1])]
                           for k, v in (ranges or {}).items()
                           if v is not None and len(v) == 2}
    _apply_progressive_filters(widget)


def _viewport_handler(widget, content, buffers):
    try:
        if not isinstance(content, dict):
            return
        t = content.get("type")
        if t == "viewport":
            if content.get("reset"):
                # double-click reset clears the selection too (matches the client
                # deselect); drop the persisted lasso so the snap-back overview
                # doesn't re-draw/re-highlight stale cells. Propagate to ALL linked
                # panels (a reset on one panel resets the group).
                vp = getattr(widget, "_vp", None)
                if vp is not None:
                    for _panel in (vp.get("group") or [widget]):
                        _pvp = getattr(_panel, "_vp", None)
                        if _pvp is not None:
                            _pvp["sel"] = set()
                        if _panel is not widget:
                            try:
                                _panel.send({"type": "vp_select", "select": []})
                            except Exception:
                                pass
            _viewport_payload(widget, content["bounds"], seq=content.get("seq"))
        elif t == "lasso":
            _lasso_payload(widget, content.get("polygon"))
        elif t == "deselect":
            # Authoritative, in-band selection clear (same channel as viewport msgs):
            # zero vp["sel"] group-wide AND push a clearing vp_select [] so a stale
            # pan/zoom redraw can't leave the old region (or the overview-union "select
            # all") highlighted after the user deselected.
            vp = getattr(widget, "_vp", None)
            if vp is not None:
                for _panel in (vp.get("group") or [widget]):
                    _pvp = getattr(_panel, "_vp", None)
                    if _pvp is not None:
                        _pvp["sel"] = set()
                        _pvp["showing_overview"] = False
                    try:
                        _panel.send({"type": "vp_select", "select": []})
                    except Exception:
                        pass
        elif t == "legend_filter":
            _legend_filter_payload(widget, content.get("cats"))
        elif t == "range_filter":
            _range_filter_payload(widget, content.get("ranges"))
    except Exception:
        # Detail-on-zoom must never wedge the UI: log for diagnosis (silent by
        # default; surfaces under logging.basicConfig(level=DEBUG)) and clear the
        # client's loading spinner so it doesn't spin forever on a failed fetch.
        _log.debug("viewport handler failed", exc_info=True)
        try:
            widget.send({"type": "vp_noop"})
        except Exception:
            pass
    finally:
        # ANY processed message (incl. early-returns / errors) marks the heavy-work
        # request done, so a synchronous read can never hang to the cap waiting on a
        # handler that skipped its ack. A plain pan didn't bump _work_req, so this is
        # a no-op there; the front-end keeps _work_req ahead only until its message
        # is processed (FIFO: the bump rides immediately before the work message).
        _ack_work(widget)


def _resolve_filter_by(filter_by, data):
    """Allow ``filter_by`` to be obs/column NAMES (a str or list of str), or ``True``
    for every numeric obs column — not only a pre-built dict / DataFrame."""
    if filter_by is None or filter_by is False or isinstance(filter_by, (pd.DataFrame, dict)):
        return None if filter_by is False else filter_by
    obs = getattr(data, "obs", None)
    table = obs if obs is not None else (data if isinstance(data, pd.DataFrame) else None)
    if table is None:
        return filter_by                                  # no table to resolve names against
    if filter_by is True:                                 # every numeric obs column
        num = table.select_dtypes(include="number")
        return num if num.shape[1] else None
    if isinstance(filter_by, str):
        filter_by = [filter_by]
    if (isinstance(filter_by, (list, tuple)) and filter_by
            and all(isinstance(c, str) for c in filter_by)):
        missing = [c for c in filter_by if c not in table.columns]
        if missing:
            raise _not_found("filter_by", missing[0], list(table.columns),
                             searched="obs / DataFrame columns")
        return table[list(filter_by)]
    return filter_by                                      # arrays etc. -> validated downstream


def _rgb01(c, fallback=(0.0, 0.0, 0.0)):
    """A colour name / hex string -> normalised (r, g, b) in 0..1 for a GL uniform."""
    try:
        import matplotlib.colors as mcolors
        return list(mcolors.to_rgb(c))
    except Exception:
        pass
    s = str(c).strip().lower()
    _named = {"black": (0, 0, 0), "white": (1, 1, 1), "gray": (.5, .5, .5),
              "grey": (.5, .5, .5), "lightgray": (.83, .83, .83),
              "lightgrey": (.83, .83, .83), "none": tuple(fallback)}
    if s in _named:
        return list(_named[s])
    if s.startswith("#") and len(s) in (4, 7):
        if len(s) == 4:
            s = "#" + "".join(ch * 2 for ch in s[1:])
        return [int(s[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    return list(fallback)


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
    add_outline: bool = False,      # scanpy-style outline (dark halo) around each point
    outline_width: Sequence[float] = (0.3, 0.05),   # (outline, gap) as fractions of the point radius
    outline_color: Sequence = ("black", "white"),   # (outline colour, gap colour); gap defaults to the background
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
    alpha: Optional[float] = None,   # scanpy alias for global opacity (0-1)
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
    progressive: bool = False,      # density-sketch overview + full detail as you zoom in (for >~4M points)
    progressive_opts: Optional[dict] = None,  # tuning: {"detail_max_points": int, "overscan": float}
    _fast: bool = False,            # internal: binary channel transfer (auto for live regl widgets)
    _detail_on_zoom: bool = False,  # internal: emit viewport messages so the kernel re-renders in-view cells
    _performance_mode: Optional[bool] = None,  # internal: squares+no-blend; None=auto (n>500k / progressive)
    _color_categories=None,         # internal: pin the full categorical level set (detail-on-zoom colour consistency)
    _pd_data: Any = None,           # internal: pre-extracted PlotData (skip extract on each detail-on-zoom viewport)
    _size_values: Any = _UNSET,     # internal: pre-resolved size channel (skip _resolve_numeric)
    _opacity_values: Any = _UNSET,  # internal: pre-resolved opacity channel
    _tooltip_cols: Any = _UNSET,    # internal: pre-resolved tooltip dict
    _filter_df: Any = _UNSET,       # internal: pre-sliced filter DataFrame
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
        Pass a **list** (e.g. ``["umap", "tsne"]``) for a linked multi-panel grid,
        one panel per embedding. A list for ``basis`` AND ``color_by`` renders their
        cross-product (one panel per basis×colour), capped at 16 panels.
    alpha
        Global point opacity (0–1), a scanpy alias for ``opacity`` applied to every
        panel. For PER-POINT opacity use ``opacity_by``; ``alpha`` is a single
        number and does not create panels.
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
        Numeric columns shown as interactive range-filter sliders. Give a
        dict / DataFrame, **or just obs / column name(s)** (a str or list of str)
        to pull them from ``adata.obs`` automatically, or ``True`` for every
        numeric ``obs`` column.
    max_points
        Cap on points actually drawn (huge data stays interactive). ``"auto"``
        (default) caps at 500k via a density-preserving subsample; ``None`` draws
        every point (ABC-Atlas style, smooth to ~4M); an int sets a custom cap.
        The plot is captioned ``"X of Y shown"`` and an automatic cap warns once.
    subsample
        ``"density"`` (default, atlas-safe: keeps rare cells) or ``"random"``.
    add_outline
        Draw a scanpy-style outline (a dark ring + a background-coloured gap)
        around every point, to make clusters pop. ``outline_width=(outline, gap)``
        are fractions of the point radius and ``outline_color=(outline, gap)`` set
        the two colours (the gap defaults to the background). Rendered in a single
        shader pass — no extra points, no performance cost — but the ring is only
        visible with reasonably large points (use it for small/medium plots, not
        20M-cell atlases where points are ~1px).
    interactive
        ``True`` returns the live, kernel-linked widget (needed for
        ``w.selection`` / ``subset`` / ``annotate`` / linked ``compose``).
        The default is a self-contained static snapshot.
    progressive
        For datasets beyond ~4M: show a density-sketch overview and re-render all
        cells inside the viewport as you zoom in (detail-on-zoom, no
        preprocessing). Always uses the live widget.
    progressive_opts
        Tuning dict for ``progressive``: ``detail_max_points`` (max points per
        zoomed-in viewport; default ``max_points``/500k) and ``overscan``
        (fraction of margin fetched around the view; default ``0.6``).
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
    # Back-compat: a few args were made internal / folded into progressive_opts.
    # The old top-level names now land in **backend_kwargs; accept them (with a
    # deprecation warning) instead of letting them vanish silently into the backend.
    _DEPRECATED = ("fast", "detail_on_zoom", "performance_mode",
                   "detail_max_points", "overscan")
    _dep = {k: backend_kwargs.pop(k) for k in _DEPRECATED if k in backend_kwargs}
    if _dep:
        import warnings
        warnings.warn(
            "These scatterplot() arguments changed: 'fast'/'performance_mode'/"
            "'detail_on_zoom' are now internal, and 'detail_max_points'/'overscan' "
            "moved into progressive_opts={'detail_max_points': ..., 'overscan': ...}. "
            "The old names still work for now but are deprecated.",
            DeprecationWarning, stacklevel=2,
        )
        if "fast" in _dep:
            _fast = bool(_dep["fast"])
        if "detail_on_zoom" in _dep:
            _detail_on_zoom = bool(_dep["detail_on_zoom"])
        if "performance_mode" in _dep:
            _performance_mode = _dep["performance_mode"]
        if "detail_max_points" in _dep:
            progressive_opts = {**(progressive_opts or {}),
                                "detail_max_points": _dep["detail_max_points"]}
        if "overscan" in _dep:
            progressive_opts = {**(progressive_opts or {}), "overscan": _dep["overscan"]}

    # Fail fast on a misspelled / unknown argument (e.g. interactive_mode, color_bys)
    # instead of silently swallowing it and plotting the wrong thing. Only the
    # jscatter backend legitimately forwards extra kwargs.
    if backend != "jscatter" and backend_kwargs:
        import difflib
        import inspect
        _valid = [p for p, par in inspect.signature(scatterplot).parameters.items()
                  if par.kind is not par.VAR_KEYWORD and not p.startswith("_")]
        _bad = next(iter(backend_kwargs))
        _sugg = difflib.get_close_matches(_bad, _valid, n=1)
        raise TypeError(
            f"scatterplot() got an unexpected keyword argument {_bad!r}"
            + (f". Did you mean {_sugg[0]!r}?" if _sugg else "")
            + " (extra keyword arguments are only forwarded when backend='jscatter')."
        )

    # `_detail_on_zoom` is the INTERNAL engine behind `progressive=` (it wires up
    # the viewport round-trip). A user reaching it directly — without progressive —
    # would get viewport messages with no handler, i.e. an endless spinner. Promote
    # it to `progressive=True`. `plot_id` is set only by the internal overview
    # recursion, which must NOT re-enter this guard (would recurse forever).
    if _detail_on_zoom and not progressive and plot_id is None:
        progressive = True
        _detail_on_zoom = False

    del _DEPRECATED, _dep
    # snapshot the call args up front so progressive mode can re-render the FULL
    # data through this exact same pipeline (keeps lasso/filter/tooltip identical)
    _params = {k: v for k, v in locals().items() if k != "data"}

    # progressive tuning knobs (folded into one dict to keep the signature lean)
    _popts = progressive_opts or {}
    _detail_max_points = _popts.get("detail_max_points")
    _overscan = float(_popts.get("overscan", 0.6))

    # filter_by may be obs/column NAMES (str / list) or True (all numeric obs) — pull
    # the actual columns now so you don't have to hand-build a DataFrame.
    filter_by = _resolve_filter_by(filter_by, data)

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
    if alpha is not None:               # scanpy alias: global point opacity
        if isinstance(alpha, (list, tuple)):
            raise TypeError(
                "alpha must be a single number (global opacity 0-1); for per-point "
                "opacity use opacity_by=, and it does not create panels."
            )
        opacity = float(alpha)
    if components is not _UNSET and components is not None:
        dims = tuple(int(c) - 1 for c in components)   # 1-based -> 0-based
    _auto_max = (max_points == "auto")
    if _auto_max:
        max_points = _DEFAULT_MAX_POINTS   # smooth by default; pass max_points=None for all points
    if _fast:
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

    # A single-element colour/basis list is just one plot — collapse it BEFORE the
    # progressive/fan-out checks so color=["x"] / basis=["umap"] take the normal path.
    if _is_name_list(color_by) and len(color_by) == 1:
        color_by = color_by[0]
    if isinstance(basis, (list, tuple)) and len(basis) == 1:
        basis = basis[0]

    # --- progressive: subset now (instant), full streamed in the background ---
    if progressive and not _is_name_list(color_by) and not isinstance(basis, (list, tuple)) and backend == "regl":
        import uuid
        interactive = True
        # detail-on-zoom ("in-memory tiling"): show a density-sketch overview, then
        # re-render ALL cells inside the viewport as you zoom in (the view holds few
        # cells when zoomed, so we can draw them all -> full detail + complete lasso,
        # with NO preprocessing). Shared plotId so each re-render replaces the plot.
        budget = (_detail_max_points if (isinstance(_detail_max_points, int) and _detail_max_points)
                  else (max_points if (isinstance(max_points, int) and max_points)
                        else _DEFAULT_MAX_POINTS))
        # A lighter overview re-renders fast on zoom-out (the heaviest refresh);
        # zoom-in detail still uses the full budget.
        overview_budget = min(budget, 250_000)
        pid = plot_id or ("rs_" + uuid.uuid4().hex[:10])
        bk = dict(_params.pop("backend_kwargs", {}) or {})
        # base64 channels (fast=False) so in-view updates can ship over a custom
        # message and be applied via plot.draw (no _spec change -> no re-render).
        base = {**_params, "progressive": False, "interactive": True, "show": False,
                "save": None,   # save once on the final widget, not the overview recursion
                "_fast": True, "plot_id": pid, "_detail_on_zoom": True}   # binary channels
        w = scatterplot(data, **{**base, "max_points": overview_budget}, **bk)   # overview
        try:
            _io = _is_anndata(data) or _is_mudata(data) or _is_spatialdata(data)
            _eff = basis if (_io and basis is not None) else x
            _full = extract(data, x=_eff, y=y, color_by=color_by, group_by=group_by,
                            layer=layer, use_raw=use_raw, dims=dims, table=table)
            # pre-resolve the other per-point channels ONCE so each viewport just
            # numpy-indexes them (no AnnData re-slice / re-resolution per zoom).
            _full_size = _resolve_numeric(size_by, data, layer, "size_by", use_raw)
            _full_opacity = _resolve_numeric(opacity_by, data, layer, "opacity_by", use_raw)
            _full_tooltip = _resolve_cols(tooltip_by, data)
            _full_filter = None if filter_by is None else pd.DataFrame(filter_by)
            _lg = w._spec.get("legend") or {}
            _fx = np.asarray(_full.x, "float64")
            _fy = np.asarray(_full.y, "float64")
            # full categorical colour (codes+levels) so a legend filter can resolve to
            # ORIGINAL cells and sync across linked panels.
            _ccodes = _clevels = None
            if _full.color is not None and _lg.get("var_type") == "categorical":
                # The overview spec's legend was derived from the density-sketch
                # SUBSAMPLE, so its LEVELS can miss categories absent from the sketch
                # and its COUNTS reflect only the subsample. Rebuild the legend from
                # the FULL dataset (same NA/groups/palette rules as build_payload) so
                # names, colours and counts are all real. This also gives the legend
                # filter the full level/code set. The JS reads the legend once from
                # the spec (never recomputed per-viewport), so this sticks.
                from ._payload import _resolve_categorical_palette
                _s = pd.Series(_full.color).astype("object")
                if _s.isna().any():
                    _s = _s.where(_s.notna(), "NA")
                _ccat = pd.Categorical(_s.astype(str))
                _ccodes = _ccat.codes
                _clevels = [str(c) for c in _ccat.categories]
                _cols = _resolve_categorical_palette(
                    _clevels, custom_colors, custom_palette, categorical_palette)
                _keep = None if groups is None else {str(g) for g in groups}
                for _i, _lv in enumerate(_clevels):
                    if _lv == "NA" or (_keep is not None and _lv not in _keep):
                        _cols[_i] = na_color
                _fc = (pd.Series(_ccat).value_counts()
                       .reindex(_ccat.categories).fillna(0).astype("int64"))
                _lg = {**_lg, "names": _clevels, "colors": _cols,
                       "counts": [int(c) for c in _fc.to_numpy()]}
                w._spec = {**w._spec, "legend": _lg}
            w._vp = {"data": data, "base": base, "bk": bk, "budget": budget,
                     "x": _fx, "y": _fy,
                     # full pre-extracted channels -> each viewport numpy-indexes
                     # these instead of re-slicing the AnnData + re-resolving color.
                     "full_pd": _full, "full_size": _full_size,
                     "full_opacity": _full_opacity, "full_tooltip": _full_tooltip,
                     "full_filter": _full_filter,
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
                     "pad": float(_overscan),   # overscan margin (tunable)
                     # legend-filter sync across linked panels (original-cell based)
                     "color_codes": _ccodes, "color_levels": _clevels,
                     "filter_keep": None, "group": None,
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
                    vp["showing_overview"] = False   # rebuild overview∪selection next zoom-out
                except Exception:
                    _log.debug("selection observer failed", exc_info=True)
            w.observe(_on_user_sel, names="_selection")

            def _mk_grid(_w=w, _fx=_fx, _fy=_fy):
                try:
                    _w._vp["grid"] = _build_grid_index(_fx, _fy)
                except Exception:
                    _log.debug("grid index build failed", exc_info=True)
            import threading
            threading.Thread(target=_mk_grid, daemon=True).start()
        except Exception:
            _log.debug("progressive setup failed; falling back to plain widget",
                       exc_info=True)
        if save:
            # save= was silently ignored for progressive plots. A static HTML has
            # no kernel, so detail-on-zoom can't work there — write an honest
            # static overview snapshot (non-binary so it's JSON-serializable).
            import warnings
            warnings.warn(
                "save= with progressive=True writes a static overview snapshot; "
                "the saved HTML has no kernel, so detail-on-zoom is unavailable.",
                stacklevel=2,
            )
            try:
                snap = scatterplot(data, **{**base, "_fast": False,
                                            "interactive": False,
                                            "_detail_on_zoom": False,
                                            "max_points": overview_budget}, **bk)
                _maybe_save(snap, save)
            except Exception:
                _log.debug("progressive save snapshot failed", exc_info=True)
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

    # --- basis and/or color_by as lists -> a linked grid over their cross-product --
    basis_multi = isinstance(basis, (list, tuple))
    color_multi = _is_name_list(color_by)
    if basis_multi or color_multi:
        if backend != "regl":
            raise ValueError(
                "A multi-panel grid (a list for basis= and/or color_by=) requires "
                "backend='regl'."
            )
        from ._compose import compose

        basis_list = list(basis) if basis_multi else [basis]
        color_list = list(color_by) if color_multi else [color_by]
        n_combo = len(basis_list) * len(color_list)
        _CAP = 16
        if n_combo > _CAP:
            raise ValueError(
                f"{len(basis_list)} basis × {len(color_list)} color_by = {n_combo} "
                f"panels exceeds the {_CAP}-panel cap. Shorten one of the lists "
                "(or call scatterplot separately per embedding)."
            )

        def _panel_title(b, c):
            parts = []
            if basis_multi:
                parts.append(str(b))
            if color_multi:
                parts.append(str(c))
            return title or " · ".join(parts) or None

        panels = [
            scatterplot(
                data, basis=b, x=x, y=y, color_by=c, group_by=group_by,
                layer=layer, use_raw=use_raw, dims=dims, table=table, point_size=point_size,
                opacity=opacity, point_color=point_color, size_by=size_by,
                opacity_by=opacity_by, tooltip_by=tooltip_by,
                pixel_ratio=pixel_ratio, categorical_palette=categorical_palette,
                continuous_palette=continuous_palette, custom_palette=custom_palette,
                custom_colors=custom_colors, vmin=vmin, vmax=vmax,
                center_zero=center_zero, na_color=na_color, groups=groups,
                add_outline=add_outline, outline_width=outline_width,
                outline_color=outline_color,
                # Linked panels are the SAME cells and the widget syncs selection +
                # filters POSITIONALLY across them, so every panel MUST draw cells in
                # the same order. A per-colour value-sort (sort_order=True) would give
                # each panel a different order -> a lasso/filter in one panel would hit
                # the wrong cells in the others. Force a colour-independent shared order
                # (subsample/random_state are identical across panels; only the value-
                # sort differs, so we drop it here). NB: different bases have different
                # coordinates but the SAME cell rows, so positional sync still holds.
                sort_order=False, random_state=random_state, max_points=max_points,
                subsample=subsample,
                progressive=progressive, progressive_opts=progressive_opts,
                title=_panel_title(b, c), xlab=xlab, ylab=ylab,
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
                # interactive=True -> linked live GridBox (needs the widgets frontend);
                # default static -> an HTML iframe-grid that renders anywhere (e.g. an
                # HPC JupyterLab without the widget manager), just not camera-linked.
                interactive=interactive, show=False,
                **backend_kwargs,
            )
            # cross-product: rows = basis, cols = color (panels fill row-major)
            for b in basis_list
            for c in color_list
        ]
        # when BOTH are lists, width = #colours so each basis is its own row;
        # otherwise let compose pick a near-square layout.
        _cols = ncols or (len(color_list) if (basis_multi and color_multi) else None)
        grid = compose(panels, cols=_cols)
        _maybe_save(grid, save)
        return grid

    # --- resolve the effective embedding (basis is preferred; x is alias) ---
    if _pd_data is not None:
        # detail-on-zoom fast path: channels were pre-extracted once and are
        # indexed per viewport (skips re-slicing the AnnData + re-resolving color,
        # the dominant per-viewport cost on big sparse atlases).
        pd_data = _pd_data
    else:
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
    if add_outline:
        # The outline ring is a dark disc drawn BEHIND each point; a semi-transparent
        # body lets it bleed through and muddies/darkens the fill (the gap pass uses
        # the canvas background, which is transparent here, so it can't mask it).
        # Render bodies opaque so the colour stays true and only the ring shows.
        opacity = 1.0

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

    size_values = (_resolve_numeric(size_by, data, layer, "size_by", use_raw)
                   if _size_values is _UNSET else _size_values)
    opacity_values = (_resolve_numeric(opacity_by, data, layer, "opacity_by", use_raw)
                      if _opacity_values is _UNSET else _opacity_values)
    tooltip_cols = (_resolve_cols(tooltip_by, data)
                    if _tooltip_cols is _UNSET else _tooltip_cols)
    if _filter_df is not _UNSET:
        filter_by = _filter_df   # pre-sliced to the in-view cells
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
        binary=_fast,
    )
    # Draw order (original cell index per drawn position). Lets the client sync a
    # selection across LINKED panels by ORIGINAL cell even when panels are drawn in
    # different orders (e.g. compose([a, b]) coloured by different variables).
    if draw_order is not None:
        spec["drawOrder"] = np.asarray(draw_order).astype("int32").tolist()
    if _detail_on_zoom:
        spec["detailOnZoom"] = True   # client emits viewport msgs -> kernel re-renders in-view cells
        # squares + no alpha-blend are much cheaper per draw; default on for the
        # zoom loop, but honour an explicit _performance_mode=False (round circles).
        spec["performanceMode"] = True if _performance_mode is None else bool(_performance_mode)
    elif _performance_mode is not None:
        spec["performanceMode"] = bool(_performance_mode)

    if add_outline:
        # scanpy-style add_outline: a crisp ring + background gap behind EVERY point,
        # using the engine's antialiased outline passes (the look of a selection).
        # It's 2 extra full-cloud draw passes, so cap it at a size where the ring is
        # actually visible AND cheap — above that it's slow and invisible (~1px), so
        # skip with a warning (use w.highlight() to mark a subset instead).
        _rendered = n if draw_order is None else int(np.asarray(draw_order).size)
        if _rendered > 150_000:
            import warnings
            warnings.warn(
                f"add_outline skipped: {_rendered:,} drawn points is too many to "
                "outline without slowing rendering (and the ring is invisible at "
                "that scale). Use it on smaller plots, or w.highlight() to mark a "
                "subset of cells.",
                stacklevel=2,
            )
        else:
            ow = outline_width or (0.3, 0.05)
            oc = outline_color or ("black", "white")
            spec["addOutline"] = True
            spec["performanceMode"] = False   # need round points for a round ring
            spec["outlineAllWidth"] = max(1.0, float(ow[0]) * float(point_size or 4))
            _rgb = _rgb01(oc[0])
            spec["outlineColor"] = "#%02x%02x%02x" % tuple(
                int(round(max(0.0, min(1.0, c)) * 255)) for c in _rgb)

    # Be honest about subsampling: caption the plot ("X of Y shown") and, when the
    # downsample was automatic (not user-requested), warn — so a subsampled plot is
    # never mistaken for the full data.
    _rendered = n if draw_order is None else int(np.asarray(draw_order).size)
    if _rendered < n:
        spec["caption"] = f"{_rendered:,} of {n:,} shown"
        if _auto_max and not _detail_on_zoom:
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
