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


def compose(plots: Sequence, cols: Optional[int] = None, sync: bool = True):
    """Lay out scatterplot widgets in a linked grid.

    Parameters
    ----------
    plots
        A list of widgets returned by :func:`scatterplot`.
    cols
        Number of columns (defaults to a near-square grid).
    sync
        When ``True`` (default), pan/zoom and lasso selection are synchronised
        across all plots.

    Returns
    -------
    An ``ipywidgets.GridBox`` containing the linked plots.
    """
    try:
        import ipywidgets
    except ModuleNotFoundError as exc:  # pragma: no cover - anywidget pulls ipywidgets
        raise ModuleNotFoundError("compose() needs ipywidgets.") from exc

    plots = list(plots)
    if not plots:
        raise ValueError("compose() needs at least one plot.")

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

    ids = [f"sp_compose_{i}" for i in range(len(plots))]
    group = ids if sync else None
    for w, pid in zip(plots, ids):
        spec = dict(getattr(w, "_spec", {}) or {})
        spec["plotId"] = pid
        spec["syncPlots"] = group
        spec["syncState"] = bool(sync)
        w._width = 0   # panels are responsive so they fit their grid column
        w._spec = spec  # re-renders the widget with the sync wiring

    # progressive panels coordinate selection/filter across the group by ORIGINAL
    # cell index (each maps to its own draw order) — give each a ref to the siblings.
    if sync:
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
        def _primary(self):
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
