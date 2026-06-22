"""Plot objects: a static (no-widget) default and a live anywidget.

``scatterplot()`` returns one of two objects that share the same analysis API
(``selection`` / ``subset`` / ``annotate`` / ``composition`` /
``diff_expression`` / ``to_html``):

* :class:`StaticPlot` (default) — a plain Python object, **not** an ipywidget.
  Its ``_repr_mimebundle_`` emits a self-contained ``<iframe srcdoc>`` snapshot,
  so it renders in JupyterLab / Notebook 7 / VS Code and **survives reopening
  with no kernel** (like a plotly figure). Because it is not a widget, nothing
  is written to the notebook's widget-state, so the ``.ipynb`` stays small.
  Trade-off: no live Python round-trip (``w.selection`` is empty unless a future
  local bridge fills it).

* ``ReglScatter`` (``interactive=True``) — the anywidget that drives the shared
  reglScatterplot bundle over the kernel comm, so ``w.selection`` round-trips
  live (needed for the single-cell workflows and for linked ``compose`` grids).

Both render the *same* compiled ``static/widget.js`` (a byte-for-byte port of the
R payload), so plots look identical across R and Python.
"""

from __future__ import annotations

import pathlib

__all__ = ["ReglScatter", "StaticPlot", "is_live_widget"]

_STATIC = pathlib.Path(__file__).parent / "static" / "widget.js"


def _make_classes():
    try:
        import anywidget
        import traitlets
    except ModuleNotFoundError as exc:  # pragma: no cover - import guard
        raise ModuleNotFoundError(
            "reglscatterpy's renderer needs 'anywidget'. "
            "Install with: pip install reglscatterpy"
        ) from exc

    class _PlotAPI:
        """Analysis API shared by the static plot and the live widget.

        Operates on ``self._spec``, ``self._source`` and ``self._selection``
        (a synced trait on the live widget, a plain list on the static plot).
        """

        def to_html(self, path, title="reglscatterpy plot"):
            """Save this plot as a standalone, offline HTML file (like R's
            ``htmlwidgets::saveWidget``) — inlines the bundle and data, stays
            interactive with no kernel."""
            from ._export import save_html

            return save_html(self, path, title=title)

        @property
        def selection(self):
            """Indices of the lasso-selected points (read or assign), always in
            **data order** — translated through the draw-order permutation when
            the plot was z-ordered (sort_order / random_state).

            Live (``interactive=True``) only — on a static plot this stays empty
            because there is no kernel link.
            """
            sel = list(self._selection)
            perm = getattr(self, "_draw_order", None)
            if perm is not None:
                return [int(perm[p]) for p in sel if 0 <= p < len(perm)]
            return [int(p) for p in sel]

        @selection.setter
        def selection(self, indices):
            idx = [int(i) for i in (indices or [])]
            perm = getattr(self, "_draw_order", None)
            if perm is not None:
                inv = getattr(self, "_inv_draw_order", None)
                if inv is None:
                    # original data index -> rendered position; a subsample omits
                    # most originals, so use a dict and drop the not-rendered ones.
                    inv = {int(o): p for p, o in enumerate(perm)}
                    self._inv_draw_order = inv
                idx = [inv[d] for d in idx if d in inv]
            self._selection = idx

        def subset(self, selection=None):
            """The source object subset to the selected cells (``adata[w.selection]``)."""
            sel = self.selection if selection is None else [int(i) for i in selection]
            src = getattr(self, "_source", None)
            if src is None:
                raise ValueError("This plot has no source object to subset.")
            if hasattr(src, "obs"):       # AnnData / MuData
                return src[sel]
            return src.iloc[sel]          # DataFrame

        def annotate(self, key, label, selection=None):
            """Write ``label`` onto the lasso-selected cells in ``obs[key]`` /
            column ``key`` of the source object. Returns the annotated object."""
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

        def composition(self, by, selection=None, normalize=True):
            """Count + fraction of the selected cells in each category of ``by``."""
            import pandas as pd

            sel = self.selection if selection is None else [int(i) for i in selection]
            if not sel:
                raise ValueError("Nothing selected - lasso some points first.")
            src = getattr(self, "_source", None)
            if src is None:
                raise ValueError("This plot has no source object.")
            frame = src.obs if hasattr(src, "obs") else src
            sub = frame.iloc[sel]
            counts = sub[by].value_counts(dropna=False)
            out = pd.DataFrame({"count": counts})
            if normalize:
                out["fraction"] = counts / counts.sum()
            return out

        def diff_expression(self, group_a=None, group_b=None, n=10, layer=None,
                            method="wilcoxon"):
            """Top differential genes between two cell groups.

            ``group_a`` defaults to the lasso selection; ``group_b`` to the rest.
            When **scanpy** is installed (and the source is an AnnData) this runs
            ``sc.tl.rank_genes_groups`` on a copy and returns its result frame
            (names / scores / logfoldchanges / pvals / pvals_adj). Otherwise it
            falls back to a Welch t-test. AnnData/MuData only.
            """
            import numpy as np
            import pandas as pd

            src = getattr(self, "_source", None)
            if src is None or not hasattr(src, "X"):
                raise TypeError("diff_expression() needs an AnnData/MuData with .X.")
            n_obs = src.n_obs
            a_idx = self.selection if group_a is None else [int(i) for i in group_a]
            if not a_idx:
                raise ValueError("Group A is empty - lasso some cells first.")
            labels = np.array(["rest"] * n_obs, dtype=object)
            labels[a_idx] = "A"
            ref = "rest"
            if group_b is not None:
                labels[[int(i) for i in group_b]] = "B"
                ref = "B"

            # Preferred path: scanpy's rank_genes_groups (AnnData only).
            if type(src).__name__ == "AnnData":
                try:
                    import scanpy as sc

                    ad = src.copy()
                    ad.obs["_rs_grp"] = pd.Categorical(labels)
                    sc.tl.rank_genes_groups(
                        ad, "_rs_grp", groups=["A"], reference=ref,
                        method=method, layer=layer, n_genes=n,
                    )
                    return sc.get.rank_genes_groups_df(ad, group="A").head(n).reset_index(drop=True)
                except ImportError:
                    pass  # fall back to the built-in test below

            a_mask = labels == "A"
            b_mask = labels == ref
            X = src.layers[layer] if layer else src.X
            Xa, Xb = X[a_mask], X[b_mask]
            if hasattr(Xa, "toarray"):
                Xa, Xb = Xa.toarray(), Xb.toarray()
            Xa, Xb = np.asarray(Xa, dtype="float64"), np.asarray(Xb, dtype="float64")
            ma, mb = Xa.mean(0), Xb.mean(0)
            lfc = np.log2((ma + 1e-9) / (mb + 1e-9))
            try:
                from scipy import stats
                stat, pval = stats.ttest_ind(Xa, Xb, axis=0, equal_var=False)
            except Exception:  # pragma: no cover - scipy optional
                denom = Xa.std(0) + Xb.std(0) + 1e-9
                stat, pval = (ma - mb) / denom, np.full(ma.shape, np.nan)
            res = pd.DataFrame({
                "gene": np.asarray(src.var_names),
                "logFC": lfc, "stat": stat, "pval": pval,
                "mean_A": ma, "mean_B": mb,
            })
            res = res.reindex(res["stat"].abs().sort_values(ascending=False).index)
            return res.head(n).reset_index(drop=True)

        def __repr__(self):
            spec = getattr(self, "_spec", None) or {}
            n = spec.get("n_points")
            cap = spec.get("caption")
            by = spec.get("colorVar") or spec.get("groupVar")
            bits = ["reglscatterpy plot"]
            if cap:
                bits.append(cap)            # "X of Y shown" when subsampled
            elif n is not None:
                bits.append(f"{n:,} points")
            if by:
                bits.append(f"color_by={by!r}")
            return "<" + ", ".join(bits) + ">"

        def _export_mimebundle(self):
            """Static (no-comm) mimebundle, honouring report/record export modes."""
            from . import _export

            if _export._report_repr_enabled():
                html = _export.report_fragment(self)
            elif _export._record_enabled():
                html = _export.record_fragment(self)
            else:
                html = _export.iframe_srcdoc(self)
            return {"text/html": html, "text/plain": repr(self)}

    class StaticPlot(_PlotAPI):
        """Default plot object: a self-contained iframe snapshot, not a widget."""

        def __init__(self, spec=None, source=None, height=500, width=0):
            self._spec = dict(spec or {})
            self._source = source
            self._height = int(height)
            self._width = int(width)
            self._selection = []

        def update(self, spec):
            self._spec = dict(spec)
            return self

        def _repr_mimebundle_(self, **kwargs):
            return self._export_mimebundle()

    class ReglScatter(anywidget.AnyWidget, _PlotAPI):
        """Live anywidget (``interactive=True``): kernel-linked, w.selection round-trips."""

        _esm = _STATIC
        _spec = traitlets.Dict().tag(sync=True)
        _height = traitlets.Int(500).tag(sync=True)
        _width = traitlets.Int(0).tag(sync=True)
        _selection = traitlets.List(trait=traitlets.Int()).tag(sync=True)

        def update(self, spec: dict) -> "ReglScatter":
            self._spec = spec
            return self

        def _repr_mimebundle_(self, **kwargs):
            # In an export kernel (report/record mode) emit the static snapshot;
            # otherwise the live widget view for the kernel round-trip.
            from . import _export

            try:
                if _export._report_repr_enabled() or _export._record_enabled():
                    return self._export_mimebundle()
            except Exception:
                pass
            return super()._repr_mimebundle_(**kwargs)

    return StaticPlot, ReglScatter


# Lazily built so importing reglscatterpy (e.g. just for `extract`) doesn't hard
# require anywidget; the classes are created on first use.
_CLASSES = None


def _classes():
    global _CLASSES
    if _CLASSES is None:
        _CLASSES = _make_classes()
    return _CLASSES


def StaticPlot(*args, **kwargs):  # noqa: N802 - factory mimics a class
    return _classes()[0](*args, **kwargs)


def ReglScatter(*args, **kwargs):  # noqa: N802 - factory mimics a class
    return _classes()[1](*args, **kwargs)


def is_live_widget(obj) -> bool:
    """True if ``obj`` is the live anywidget (vs a static plot / other)."""
    return isinstance(obj, _classes()[1])
