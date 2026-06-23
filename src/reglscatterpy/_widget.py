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


class _SelectionResult(list):
    """The selected ORIGINAL indices — a plain ``list`` that *also* carries the
    analysis helpers, so both ``w.composition('louvain')`` and the clearer
    ``w.selection.composition('louvain')`` work (likewise ``.subset()``,
    ``.diff_expression()``, ``.annotate()``)."""

    def __init__(self, indices, widget=None):
        super().__init__(indices)
        self._w = widget

    def composition(self, by, **kw):
        return self._w.composition(by, selection=list(self), **kw)

    def subset(self):
        return self._w.subset(selection=list(self))

    def diff_expression(self, group_b=None, **kw):
        return self._w.diff_expression(group_a=list(self), group_b=group_b, **kw)

    def annotate(self, key, label):
        return self._w.annotate(key, label, selection=list(self))



def _make_classes():
    try:
        import anywidget
        import traitlets
    except ModuleNotFoundError as exc:  # pragma: no cover - import guard
        raise ModuleNotFoundError(
            "reglscatterpy's renderer needs 'anywidget'. "
            "Install with: pip install reglscatterpy"
        ) from exc

    def _alias_trap(wrong, right):
        """A property that errors on read OR assignment of a common misspelling —
        Python otherwise silently lets ``w.<wrong> = ...`` create a dead attribute."""
        msg = (f"'{wrong}' is not a reglscatterpy plot attribute — did you mean "
               f"'.{right}'? (Assigning a misspelled attribute would otherwise "
               "silently do nothing.)")

        def _raise(self, *args):
            raise AttributeError(msg)

        return property(_raise, _raise)

    class _PlotAPI:
        """Analysis API shared by the static plot and the live widget.

        Operates on ``self._spec``, ``self._source`` and ``self._selection``
        (a synced trait on the live widget, a plain list on the static plot).
        """

        # Trap the typos that silently no-op as attribute assignments.
        selected = _alias_trap("selected", "selection")
        select = _alias_trap("select", "selection")
        selections = _alias_trap("selections", "selection")
        highlights = _alias_trap("highlights", "highlight")
        filter = _alias_trap("filter", "filtered")
        filters = _alias_trap("filters", "filtered")

        def to_html(self, path, title="reglscatterpy plot"):
            """Save this plot as a standalone, offline HTML file (like R's
            ``htmlwidgets::saveWidget``) — inlines the bundle and data, stays
            interactive with no kernel. On a live (``interactive=True``) widget it
            captures the CURRENT view / selection / filter so the file opens exactly
            as you left it."""
            from ._export import save_html

            return save_html(self, path, title=title)

        def _snapshot_spec(self):
            """A copy of the spec with the CURRENT live camera / selection / filter
            baked in, so an exported HTML opens as the interactive widget looks now.
            (No live state -> identical to the plain spec.)"""
            import copy
            spec = copy.deepcopy(dict(getattr(self, "_spec", {}) or {}))
            cam = list(getattr(self, "_camera", []) or [])
            if cam:
                spec["cameraView"] = [float(v) for v in cam]
            sel = list(getattr(self, "_selection", []) or [])      # rendered positions
            if sel:
                spec["init_selected_indices"] = [int(p) for p in sel]
            if getattr(self, "_filtered_on", False):               # filtered -> positions
                keep = list(getattr(self, "_filtered", []) or [])  # original indices
                perm = getattr(self, "_draw_order", None)
                if perm is not None:
                    inv = getattr(self, "_inv_draw_order", None)
                    if inv is None:
                        inv = {int(o): p for p, o in enumerate(perm)}
                    spec["init_server_indices"] = [inv[o] for o in keep if o in inv]
                else:
                    spec["init_server_indices"] = [int(o) for o in keep]
            return spec

        def _pump(self):
            """Drain pending front-end messages so a just-made lasso / selection /
            filter is reflected BEFORE we read it back.

            ``w.selection`` is a plain Python read, but the lasso travels to the
            kernel as an async widget message; if the kernel hasn't processed it
            yet (it falls behind on big data — 10M-point polygon tests + viewport
            redraws), the read returns the *previous* state (one step stale). This
            forces the queued messages through first. No-op on a static plot or if
            ``jupyter_ui_poll`` isn't installed (then reads may lag by one step)."""
            if not callable(getattr(self, "send", None)):
                return
            try:
                from jupyter_ui_poll import ui_events
            except Exception:
                return
            try:
                import time as _t
                with ui_events() as poll:
                    poll(512)            # process whatever the front-end already sent
                    _t.sleep(0.03)       # tiny nudge to catch an in-flight message...
                    poll(512)            # ...then drain that too
            except Exception:
                pass

        @property
        def selection(self):
            """Indices of the lasso-selected points (read or assign), always in
            **data order** — translated through the draw-order permutation when
            the plot was z-ordered (sort_order / random_state).

            Live (``interactive=True``) only — on a static plot this stays empty
            because there is no kernel link.
            """
            self._pump()    # make sure a just-drawn lasso has been applied
            # detail-on-zoom keeps the logical selection (original indices) so it
            # survives viewport swaps; prefer it when present.
            vp = getattr(self, "_vp", None)
            if vp is not None and "sel" in vp:
                return _SelectionResult(sorted(int(i) for i in vp["sel"]), self)
            sel = list(self._selection)
            perm = getattr(self, "_draw_order", None)
            if perm is not None:
                return _SelectionResult((int(perm[p]) for p in sel if 0 <= p < len(perm)), self)
            return _SelectionResult((int(p) for p in sel), self)

        @property
        def filtered(self):
            """Original indices of the cells currently passing the in-plot filters
            (range sliders + legend categories). When no filter is active this is
            **all shown cells** (everything passes). Live (``interactive=True``)
            only — like :attr:`selection`, but for the filter instead of the lasso.
            """
            if not callable(getattr(self, "send", None)):
                raise AttributeError(
                    "w.filtered needs a live widget — build it with "
                    "interactive=True. A static plot has no kernel link, so the "
                    "filter state can't be read back into Python."
                )
            self._pump()    # apply any just-made filter before reading it back
            if getattr(self, "_filtered_on", False):
                return sorted(int(i) for i in getattr(self, "_filtered", []))
            # no active filter -> every shown cell passes
            do = getattr(self, "_draw_order", None)
            if do is not None:
                return sorted(int(i) for i in do)
            n = (getattr(self, "_spec", None) or {}).get("n_points")
            return list(range(int(n))) if n else []

        def _resolve_orig_indices(self, indices):
            """Coerce a selection-like input to ORIGINAL integer indices. Accepts
            integer positions, obs_names / DataFrame index labels (strings), or a
            boolean mask (length == n)."""
            import numpy as np
            if indices is None:
                return []
            seq = list(indices)
            if not seq:
                return []
            first = seq[0]
            if isinstance(first, (bool, np.bool_)):
                return [int(i) for i, b in enumerate(seq) if b]
            if isinstance(first, (int, np.integer)):
                return [int(i) for i in seq]
            # names (e.g. adata[mask].obs_names) -> positions in the source
            src = getattr(self, "_source", None)
            names = None
            if src is not None:
                if hasattr(src, "obs_names"):
                    names = list(src.obs_names)
                elif hasattr(src, "index"):
                    names = list(src.index)
            if names is None:
                raise TypeError(
                    "selecting by name needs the original AnnData/DataFrame "
                    "(obs_names / index); pass integer positions instead."
                )
            pos = {n: i for i, n in enumerate(names)}
            missing = [n for n in seq if n not in pos]
            if missing:
                raise KeyError(
                    f"{len(missing)} name(s) not found in obs_names, "
                    f"e.g. {missing[:3]}."
                )
            return [pos[n] for n in seq]

        @selection.setter
        def selection(self, indices):
            idx = self._resolve_orig_indices(indices)
            perm = getattr(self, "_draw_order", None)
            if perm is not None:
                inv = getattr(self, "_inv_draw_order", None)
                if inv is None:
                    # original data index -> rendered position; a subsample omits
                    # most originals, so use a dict and drop the not-rendered ones.
                    inv = {int(o): p for p, o in enumerate(perm)}
                    self._inv_draw_order = inv
                positions = [inv[d] for d in idx if d in inv]
            else:
                positions = idx
            vp = getattr(self, "_vp", None)
            if vp is not None and "sel" in vp:
                vp["_setting"] = True          # don't let the observer shrink it to in-view
                vp["sel"] = set(idx)
                try:
                    self._selection = positions
                finally:
                    vp["_setting"] = False
            else:
                self._selection = positions

        def highlight(self, indices, color=None):
            """Persistently mark points with a crisp ring + size bump (the engine's
            selection look) — but this is **not** the selection, so it survives a
            double-click and a new lasso. ``indices`` are original data indices;
            ``color`` sets the ring colour (a hex / CSS colour). Pass ``[]`` / ``None``
            to clear. Live (``interactive=True``) only — needs the kernel link.

            Note: with ``progressive=True`` the highlight marks the currently-shown
            cells; it doesn't yet follow new cells streamed in on zoom.
            """
            idx = self._resolve_orig_indices(indices)
            perm = getattr(self, "_draw_order", None)
            if perm is not None:
                inv = getattr(self, "_inv_draw_order", None)
                if inv is None:
                    inv = {int(o): p for p, o in enumerate(perm)}
                    self._inv_draw_order = inv
                positions = [inv[d] for d in idx if d in inv]
            else:
                positions = idx
            self._highlight = idx
            send = getattr(self, "send", None)
            if callable(send):
                msg = {"type": "hl", "points": positions}
                if color is not None:
                    msg["color"] = color
                try:
                    send(msg)
                except Exception:
                    pass
            return self

        def subset(self, selection=None):
            """The source object subset to the selected cells (``adata[w.selection]``)."""
            sel = self.selection if selection is None else self._resolve_orig_indices(selection)
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

            sel = self.selection if selection is None else self._resolve_orig_indices(selection)
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

            sel = self.selection if selection is None else self._resolve_orig_indices(selection)
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
                            method="wilcoxon", key_added=None):
            """Top differential genes between two cell groups.

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
            """
            import numpy as np
            import pandas as pd

            src = getattr(self, "_source", None)
            if src is None or not hasattr(src, "X"):
                raise TypeError("diff_expression() needs an AnnData/MuData with .X.")
            n_obs = src.n_obs
            a_idx = self.selection if group_a is None else self._resolve_orig_indices(group_a)
            if not a_idx:
                raise ValueError("Group A is empty - lasso some cells first.")
            labels = np.array(["rest"] * n_obs, dtype=object)
            labels[a_idx] = "A"
            ref = "rest"
            if group_b is not None:
                labels[self._resolve_orig_indices(group_b)] = "B"
                ref = "B"

            # Auto-save when the source has .uns (AnnData). Default to scanpy's
            # 'rank_genes_groups' key; key_added overrides; key_added=False disables.
            _save_to = (None if key_added is False or not hasattr(src, "uns")
                        else (key_added or "rank_genes_groups"))

            def _save(df, uns_native=None):
                if _save_to is not None:
                    src.uns[_save_to] = uns_native if uns_native is not None else df
                return df

            # Preferred path: scanpy's rank_genes_groups (AnnData only).
            if type(src).__name__ == "AnnData":
                try:
                    import scanpy as sc

                    ad = src.copy()
                    ad.obs["_rs_grp"] = pd.Categorical(labels)
                    _key = _save_to or "rank_genes_groups"
                    sc.tl.rank_genes_groups(
                        ad, "_rs_grp", groups=["A"], reference=ref,
                        method=method, layer=layer, n_genes=n, key_added=_key,
                    )
                    df = sc.get.rank_genes_groups_df(ad, group="A", key=_key).head(n).reset_index(drop=True)
                    # store the scanpy-NATIVE uns structure on the real adata, so
                    # sc.pl.rank_genes_groups(adata) etc. work afterwards.
                    return _save(df, uns_native=ad.uns.get(_key))
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
            return _save(res.head(n).reset_index(drop=True))

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
        _filtered = traitlets.List(trait=traitlets.Int()).tag(sync=True)
        _filtered_on = traitlets.Bool(False).tag(sync=True)
        _camera = traitlets.List(trait=traitlets.Float()).tag(sync=True)

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
