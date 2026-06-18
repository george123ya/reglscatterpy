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

        def to_html(self, path, title="reglscatterpy plot"):
            """Save this plot as a standalone, offline HTML file.

            Like R's ``htmlwidgets::saveWidget`` - the file inlines the widget
            and the plot's data, so it stays interactive after a notebook is
            closed and reopened (a kernel-free snapshot; no Python
            round-trip). ``w.to_html("umap.html")``.
            """
            from ._export import save_html

            return save_html(self, path, title=title)

        def _repr_mimebundle_(self, **kwargs):
            # In a live notebook this returns the normal interactive widget
            # view. During a report export (save_notebook_html flips report
            # mode on, in that kernel only) it instead returns a self-contained
            # text/html snapshot — and drops the widget view so nbconvert bakes
            # the static, kernel-free plot in unambiguously.
            from . import _export

            if _export._report_repr_enabled():
                try:
                    return {
                        "text/html": _export.report_fragment(self),
                        "text/plain": repr(self),
                    }
                except Exception:
                    pass
            return super()._repr_mimebundle_(**kwargs)

        def __repr__(self):
            # Shown as the text/plain fallback when a host can't render the
            # live widget view - e.g. a reopened notebook whose widget state
            # wasn't saved. A clean summary beats the default object address;
            # the hint tells the user how to get the plot back.
            spec = self._spec or {}
            n = spec.get("n_points")
            by = spec.get("colorVar") or spec.get("groupVar")
            bits = ["reglscatterpy plot"]
            if n is not None:
                bits.append(f"{n:,} points")
            if by:
                bits.append(f"color_by={by!r}")
            head = ", ".join(bits)
            return f"<{head}> (re-run this cell to render the interactive plot)"

        @property
        def selection(self):
            """Indices of the lasso-selected points (read or assign)."""
            return list(self._selection)

        @selection.setter
        def selection(self, indices):
            self._selection = [int(i) for i in (indices or [])]

        def subset(self, selection=None):
            """The source object subset to the selected cells.

            Equivalent to ``adata[w.selection]`` (the indices are positional, in
            the plotted order), returned as an ``AnnData`` / ``MuData`` view or a
            DataFrame slice — so you can keep analysing the lassoed cells:
            ``sub = w.subset(); sc.tl.rank_genes_groups(sub, ...)``.
            """
            sel = self.selection if selection is None else [int(i) for i in selection]
            src = getattr(self, "_source", None)
            if src is None:
                raise ValueError("This plot has no source object to subset.")
            if hasattr(src, "obs"):       # AnnData / MuData
                return src[sel]
            return src.iloc[sel]          # DataFrame

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

        def composition(self, by, selection=None, normalize=True):
            """Composition of the lasso-selected cells by an ``obs`` column.

            ``w.composition("leiden")`` returns a DataFrame of the count and
            fraction of the selected cells in each category of ``by`` - e.g. to
            see which clusters a lassoed region is made of. ``normalize=False``
            drops the fraction column.
            """
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

        def diff_expression(self, group_a=None, group_b=None, n=10, layer=None):
            """Top differential genes between two cell groups.

            ``group_a`` defaults to the current lasso selection; ``group_b``
            defaults to all other cells ("rest"). Pass explicit index lists to
            compare two saved selections (e.g. ``a = w.selection`` after one
            lasso, then ``w.diff_expression(a, w.selection)`` after another).
            Ranks genes by a Welch t-statistic (falls back to a standardised
            mean difference if SciPy is absent) and returns the top ``n`` by
            absolute effect, with ``logFC`` and group means. AnnData/MuData only.
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
            a_mask = np.zeros(n_obs, dtype=bool); a_mask[a_idx] = True
            if group_b is None:
                b_mask = ~a_mask
            else:
                b_mask = np.zeros(n_obs, dtype=bool); b_mask[[int(i) for i in group_b]] = True

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

    return ReglScatter


# Lazily built so importing reglscatterpy (e.g. just for `extract`) does not hard
# require anywidget; the class is created on first use.
_CLASS = None


def ReglScatter(*args, **kwargs):  # noqa: N802 - factory mimics a class
    global _CLASS
    if _CLASS is None:
        _CLASS = _make_class()
    return _CLASS(*args, **kwargs)
