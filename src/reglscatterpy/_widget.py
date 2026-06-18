"""The anywidget that renders the shared reglScatterplot widget.

``static/widget.js`` is an ESM bundle built from ``js/src/anywidget.js`` (see
the repo's ``js/`` directory). It loads the *same* compiled widget the R package
uses via a tiny ``HTMLWidgets`` shim and drives it directly, so a plot looks and
behaves identically whether created from R (htmlwidgets) or Python (anywidget):
the draggable legend, lasso, tooltips, sync and PNG/SVG/PDF export all come from
one codebase.

The Python side hands the widget a ``_spec`` dict built by
:func:`reglscatterpy._payload.build_payload`, which is a byte-for-byte port of
the R payload (locked down by ``tests/test_payload_parity``).
"""

from __future__ import annotations

import pathlib

__all__ = ["ReglScatter"]

_STATIC = pathlib.Path(__file__).parent / "static" / "widget.js"


def _make_class():
    try:
        import anywidget
        import traitlets
    except ModuleNotFoundError as exc:  # pragma: no cover - import guard
        raise ModuleNotFoundError(
            "reglscatterpy's default renderer needs 'anywidget'. "
            "Install with: pip install reglscatterpy"
        ) from exc

    class ReglScatter(anywidget.AnyWidget):
        _esm = _STATIC
        _spec = traitlets.Dict().tag(sync=True)
        _height = traitlets.Int(500).tag(sync=True)
        # 0 => responsive (100% of the cell); a positive value => fixed px width.
        _width = traitlets.Int(0).tag(sync=True)
        # Selected point indices, kept in sync both ways with the lasso.
        _selection = traitlets.List(trait=traitlets.Int()).tag(sync=True)

        def update(self, spec: dict) -> "ReglScatter":
            """Swap in a new payload and re-render in place."""
            self._spec = spec
            return self

        @property
        def selection(self):
            """Indices of the lasso-selected points (read or assign)."""
            return list(self._selection)

        @selection.setter
        def selection(self, indices):
            self._selection = [int(i) for i in (indices or [])]

        def annotate(self, key, label, selection=None):
            """Write a label onto the lasso-selected cells.

            Lasso a population in the plot, then ``w.annotate("cell_type",
            "T cells")`` writes that label into ``obs[key]`` (AnnData / MuData)
            or the column ``key`` (DataFrame) of the object this plot was made
            from, for the currently selected rows. Call repeatedly with
            different labels to build up an annotation; re-plot ``color_by=key``
            to see it. Returns the annotated object.
            """
            import numpy as np
            import pandas as pd

            sel = self.selection if selection is None else [int(i) for i in selection]
            src = getattr(self, "_source", None)
            if src is None:
                raise ValueError(
                    "This plot has no source object to annotate "
                    "(it was built from raw arrays)."
                )
            has_obs = hasattr(src, "obs")
            frame = src.obs if has_obs else src
            if not hasattr(frame, "columns"):
                raise TypeError("annotate() supports AnnData, MuData and DataFrame.")
            n = frame.shape[0]
            if key in frame.columns:
                col = np.asarray(frame[key].astype("object")).copy()
            else:
                col = np.array([None] * n, dtype=object)
            for i in sel:
                if 0 <= i < n:
                    col[i] = label
            new = pd.Categorical(col)
            if has_obs:
                src.obs[key] = new
            else:
                src[key] = new
            return src

    return ReglScatter


# Lazily built so importing reglscatterpy (e.g. just for `extract`) does not hard
# require anywidget; the class is created on first use.
_CLASS = None


def ReglScatter(*args, **kwargs):  # noqa: N802 - factory mimics a class
    global _CLASS
    if _CLASS is None:
        _CLASS = _make_class()
    return _CLASS(*args, **kwargs)
