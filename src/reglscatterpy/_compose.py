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
    from ._widget import is_live_widget
    if not all(is_live_widget(p) for p in plots):
        raise ValueError(
            "compose() needs live (interactive) plots — a linked grid syncs over "
            "the kernel. Build each plot with scatterplot(..., interactive=True), "
            "or pass a list to color_by= which does this for you."
        )
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
        w._spec = spec  # re-renders the widget with the sync wiring

    cols = cols or ceil(sqrt(len(plots)))
    return ipywidgets.GridBox(
        plots,
        layout=ipywidgets.Layout(
            grid_template_columns=f"repeat({cols}, 1fr)",
            grid_gap="8px",
        ),
    )
