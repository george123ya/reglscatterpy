"""Arrange several scatterplots in a linked grid.

``compose([...])`` lays out multiple :func:`scatterplot` widgets in a grid and
links them so panning/zooming one moves the others and a lasso selection is
mirrored across the group (handy for comparing embeddings of the same cells).

Linking reuses the widget's own ``syncPlots`` mechanism: each plot is given a
shared sync-group id, so the (shared, page-global) widget registry wires their
cameras and selections together client-side - no extra plumbing.
"""

from __future__ import annotations

from math import ceil, sqrt
from typing import Optional, Sequence

__all__ = ["compose"]


def compose(plots: Sequence, cols: Optional[int] = None, sync="auto"):
    """Lay out scatterplot widgets in a linked grid.

    Parameters
    ----------
    plots
        A list of widgets returned by :func:`scatterplot`.
    cols
        Number of columns (defaults to a near-square grid).
    sync
        ``"auto"`` (default) links pan/zoom, lasso selection and filters across
        panels **only when they share the same data object** — panels from
        *different* AnnData files stay independent (linking by cell index would be
        meaningless). ``True`` / ``False`` force linking on / off.

    Returns
    -------
    A ``GridBox`` (subclass) of the plots. When linked it exposes the group's
    ``selection`` / ``filtered`` / ``subset`` / ``composition`` / ``diff_expression``;
    either way ``grid.panels`` gives the per-panel widgets (e.g.
    ``grid.panels[0].filtered``).
    """
    try:
        import ipywidgets
    except ModuleNotFoundError as exc:  # pragma: no cover - anywidget pulls ipywidgets
        raise ModuleNotFoundError("compose() needs ipywidgets.") from exc

    plots = list(plots)
    if not plots:
        raise ValueError("compose() needs at least one plot.")

    from ._widget import ReglScatter, is_live_widget

    # Static panels (the default) -> render as an HTML iframe-grid that needs NO
    # ipywidgets frontend, so multi-panel plots show up even in an HPC / SCC
    # JupyterLab where the widget manager isn't active. (Not camera/lasso-linked —
    # linking needs live widgets; pass interactive=True / sync=True for that.)
    if sync is not True and all(not is_live_widget(p) for p in plots):
        cols = cols or ceil(sqrt(len(plots)))
        return _HtmlGrid(plots, cols)

    # A linked grid syncs over the kernel, so each panel must be a live widget.
    # Auto-upgrade plain (default static) plots so you DON'T have to pass
    # interactive=True to every scatterplot() — just compose([a, b]).
    from ._widget import ReglScatter, is_live_widget

    def _as_live(p):
        if is_live_widget(p):
            return p
        spec = dict(getattr(p, "_spec", {}) or {})
        if not spec:
            raise TypeError("compose() takes reglscatterpy plots (from scatterplot()).")
        w = ReglScatter()
        w._height = int(getattr(p, "_height", 500) or 500)
        w._width = int(getattr(p, "_width", 0) or 0)
        w._source = getattr(p, "_source", None)
        w._draw_order = getattr(p, "_draw_order", None)
        w._spec = spec
        return w

    plots = [_as_live(p) for p in plots]
    if len(plots) > 16:
        import warnings
        warnings.warn(
            f"compose() with {len(plots)} linked plots may be slow to render; "
            "consider fewer or smaller plots.",
            stacklevel=2,
        )

    # Decide whether to LINK: by default only when every panel shares the same data
    # object (same cells) — different AnnData files can't be synced by cell index.
    _src0 = getattr(plots[0], "_source", None)
    same_source = _src0 is not None and all(
        getattr(p, "_source", None) is _src0 for p in plots)
    if sync == "auto":
        do_sync = same_source
        if not do_sync and len(plots) > 1:
            import warnings
            warnings.warn(
                "compose(): panels come from different data sources, so selection / "
                "filter are NOT linked (cells aren't comparable across them). Use "
                "grid.panels[i] for each; pass sync=True to force linking if the cells "
                "do correspond.",
                stacklevel=2,
            )
    else:
        do_sync = bool(sync)
        # Forcing sync across different sources links by CELL INDEX (position i in
        # one panel == position i in another) — valid only when the panels are the
        # same cells in the same order. Warn if the sizes don't even match.
        if do_sync and not same_source:
            def _nobs(p):
                src = getattr(p, "_source", None)
                return getattr(src, "n_obs", None) or getattr(src, "shape", [None])[0]
            sizes = [_nobs(p) for p in plots]
            if len({s for s in sizes if s is not None}) > 1:
                import warnings
                warnings.warn(
                    f"compose(sync=True): panels have different sizes {sizes}; "
                    "index-based sync only lines up when the panels are the SAME "
                    "cells in the SAME order. Selection/filter will misalign here.",
                    stacklevel=2,
                )

    ids = [f"sp_compose_{i}" for i in range(len(plots))]
    group = ids if do_sync else None
    for w, pid in zip(plots, ids):
        spec = dict(getattr(w, "_spec", {}) or {})
        spec["plotId"] = pid
        spec["syncPlots"] = group
        spec["syncState"] = bool(do_sync)
        w._width = 0   # panels are responsive so they fit their grid column
        w._spec = spec  # re-renders the widget with the sync wiring

    # progressive panels coordinate selection/filter across the group by ORIGINAL
    # cell index (each maps to its own draw order) — give each a ref to the siblings.
    if do_sync:
        for w in plots:
            vp = getattr(w, "_vp", None)
            if vp is not None:
                vp["group"] = plots

    cols = cols or ceil(sqrt(len(plots)))
    grid = _get_composed_class()(
        plots,
        layout=ipywidgets.Layout(
            grid_template_columns=f"repeat({cols}, 1fr)",
            grid_gap="8px",
        ),
    )
    grid._panels = plots
    grid._synced = do_sync
    return grid


_COMPOSED_CLS = None


def _get_composed_class():
    global _COMPOSED_CLS
    if _COMPOSED_CLS is None:
        _COMPOSED_CLS = _composed_plots_class()
    return _COMPOSED_CLS


def _composed_plots_class():
    """A GridBox that also exposes the group's synced selection / filter / analysis
    (the panels are linked, so reading any one gives the group's state)."""
    import ipywidgets

    class _ComposedPlots(ipywidgets.GridBox):
        @property
        def panels(self):
            """The individual plot widgets, e.g. ``grid.panels[0].filtered``."""
            return list(getattr(self, "_panels", []))

        @property
        def _primary(self):
            if not getattr(self, "_synced", True):
                raise AttributeError(
                    "These panels are independent (different data sources), so there "
                    "is no shared selection/filter. Use grid.panels[i].selection / "
                    ".filtered for each, or compose([...], sync=True) to link them."
                )
            for p in getattr(self, "_panels", []):
                if callable(getattr(p, "send", None)):
                    return p                      # first live panel
            return self._panels[0]

        @property
        def selection(self):
            """The group's lasso selection (original indices) — assignable to all panels."""
            return self._primary.selection

        @selection.setter
        def selection(self, value):
            for p in self._panels:
                try:
                    p.selection = value
                except Exception:
                    pass

        @property
        def filtered(self):
            """Original indices passing the group's (synced) filters."""
            return self._primary.filtered

        def subset(self, *a, **k):
            return self._primary.subset(*a, **k)

        def diff_expression(self, *a, **k):
            return self._primary.diff_expression(*a, **k)

        def composition(self, *a, **k):
            return self._primary.composition(*a, **k)

        def annotate(self, *a, **k):
            return self._primary.annotate(*a, **k)

        def highlight(self, *a, **k):
            for p in self._panels:
                try:
                    p.highlight(*a, **k)
                except Exception:
                    pass
            return self

    return _ComposedPlots


class _HtmlGrid:
    """A static multi-panel grid rendered as an HTML CSS-grid of ``<iframe>``s.

    Needs NO ipywidgets frontend, so it shows up even in an HPC / BU-SCC JupyterLab
    where the widget manager isn't active (a live ``GridBox`` would there fall back
    to its text ``repr``). Panels are independent (static, no kernel link) — use
    ``grid.panels[i]`` for each; pass ``interactive=True`` for a linked live grid.
    """

    def __init__(self, panels, cols):
        self._panels = list(panels)
        self._cols = cols

    @property
    def panels(self):
        return list(self._panels)

    def _repr_mimebundle_(self, **kwargs):
        from ._export import iframe_srcdoc

        cells = "\n".join(
            iframe_srcdoc(p) for p in self._panels if getattr(p, "_spec", None)
        )
        html = (
            f'<div style="display:grid; '
            f'grid-template-columns:repeat({self._cols},1fr); gap:8px; '
            f'align-items:start;">{cells}</div>'
        )
        return {
            "text/html": html,
            "text/plain": f"reglscatterpy grid ({len(self._panels)} panels)",
        }
