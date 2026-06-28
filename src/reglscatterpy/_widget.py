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

# Keep _esm as a Path (NOT a content string): anywidget only enables live-reload
# (ANYWIDGET_HMR=1) when _esm is a FILE, and HMR is the reliable way to iterate on
# the bundle in VS Code without the webview serving a stale cached copy. (HMR also
# requires an EDITABLE install — anywidget skips watching files under site-packages.)
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

    def diff_expression_by(self, by, **kw):
        return self._w.diff_expression_by(by, selection=list(self), **kw)

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

        _warned_no_poll = False   # warn once if jupyter_ui_poll is missing

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
            """Block until the kernel has applied the LATEST front-end interaction,
            so a just-made lasso / selection / filter is reflected before we read it.

            ``w.selection`` is a synchronous Python read, but the lasso travels to
            the kernel as an async widget message; on big data the kernel falls a
            step behind (10M-point polygon tests + viewport redraws), so the read
            would return the *previous* state. Every selection/filter the front-end
            commits bumps ``_sel_gen``; the kernel echoes it to ``_applied_gen``
            once processed. Here we drain pending messages until those match —
            never stopping mid-backlog (which would land on a stale cluster). No-op
            on a static plot. Needs ``jupyter_ui_poll``; warns once if missing."""
            if not callable(getattr(self, "send", None)):
                return
            try:
                from jupyter_ui_poll import ui_events
            except Exception:
                if not _PlotAPI._warned_no_poll:
                    _PlotAPI._warned_no_poll = True
                    import warnings
                    warnings.warn(
                        "reglscatterpy: install 'jupyter-ui-poll' so w.selection / "
                        "w.filtered reflect the latest lasso immediately "
                        "(pip install jupyter-ui-poll). Without it the value can lag "
                        "by one interaction.", RuntimeWarning, stacklevel=3)
                return
            try:
                vp = getattr(self, "_vp", None)
                n = len(vp["x"]) if (vp and "x" in vp) else 0
            except Exception:
                n = 0
            # min_settle is the ARRIVAL window: the front-end bumps its request
            # counters (a tiny message sent before any heavy work), so within this
            # window the pump sees a bump and KNOWS to keep waiting — up to `cap` —
            # until the kernel acks. Without the window an immediate read sees the
            # stale (pre-bump) counters "caught up" and exits early. When nothing is
            # pending both barriers are already satisfied, so the read returns right
            # after min_settle (no over-block). cap scales with size: a 10M-point
            # select-all runs synchronously in the comm handler for seconds.
            min_settle = 0.4
            cap = min(30.0, max(3.0, (n / 1_000_000.0) * 4.0))
            import time as _t
            try:
                start = _t.monotonic()
                with ui_events() as poll:
                    while True:
                        poll(256)                 # process whatever the front-end sent
                        elapsed = _t.monotonic() - start
                        gen = int(getattr(self, "_sel_gen", 0) or 0)
                        applied = int(getattr(self, "_applied_gen", 0) or 0)
                        wreq = int(getattr(self, "_work_req", 0) or 0)
                        wdone = int(getattr(self, "_work_done", 0) or 0)
                        caught = applied >= gen and wdone >= wreq
                        if caught and elapsed >= min_settle:
                            break                 # nothing pending (or all applied)
                        if elapsed >= cap:
                            break                 # give up (kernel stuck) -> stale
                        _t.sleep(0.012)
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
            # Progressive (detail-on-zoom): the authoritative filter is the kernel's
            # ORIGINAL-cell keep over the FULL dataset (vp["filter_keep"]), NOT the
            # _filtered trait / _draw_order, which only cover the cells currently drawn
            # in the viewport. Reading those would under-count to the visible subset.
            vp = getattr(self, "_vp", None)
            if vp is not None and "x" in vp:
                import numpy as np
                keep = vp.get("filter_keep")
                if keep is not None:                       # legend/category filter active
                    return np.asarray(keep, dtype=np.int64).tolist()   # originals, ascending
                return list(range(len(vp["x"])))           # no filter -> every cell passes
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

        @property
        def colors(self):
            """The categorical colour map as ``{category: '#rrggbb'}`` (the rendered
            palette, including any in-plot recolours via the legend colorpicker).
            ``None`` for a continuous / single-colour plot. Save it scanpy-style with
            e.g. ``adata.uns['louvain_colors'] = list(w.colors.values())``."""
            lg = (getattr(self, "_spec", None) or {}).get("legend") or {}
            if lg.get("var_type") != "categorical":
                return None
            self._pump()    # pick up a just-made colorpicker edit
            # a legend colorpicker recolour (synced back) overrides the spec palette
            ec = list(getattr(self, "_legend_colors", []) or [])
            en = list(getattr(self, "_legend_names", []) or [])
            if ec and en and len(ec) == len(en):
                return dict(zip(en, ec))
            return dict(zip(lg.get("names", []), lg.get("colors", [])))

        def morph_to(self, basis, duration=1200):
            """Animate the points from their current layout to another embedding —
            e.g. ``w.morph_to('spatial')`` morphs a UMAP into the spatial layout
            (and back with ``w.morph_to('umap')``). Positions tween; colours and
            sizes stay. Live (``interactive=True``) only; not for ``progressive=True``
            (only a subset is resident). Returns ``self`` so calls can chain."""
            if not callable(getattr(self, "send", None)):
                raise AttributeError("morph_to() needs a live widget (interactive=True).")
            if getattr(self, "_vp", None) is not None:
                raise NotImplementedError(
                    "morph_to() isn't supported with progressive=True yet (only a "
                    "subset of cells is resident). Use a non-progressive plot "
                    "(e.g. max_points=) to morph between embeddings."
                )
            import numpy as np
            from ._extract import _resolve_basis
            from ._payload import _pad_range, _normalise_range, _q_u16
            src = getattr(self, "_source", None)
            if src is None or not hasattr(src, "obsm"):
                raise TypeError("morph_to() needs an AnnData-like source with .obsm embeddings.")
            key = _resolve_basis(src, basis)
            coords = np.asarray(src.obsm[key], dtype="float64")[:, :2]
            do = getattr(self, "_draw_order", None)
            if do is not None:                       # move the SAME drawn cells
                coords = coords[np.asarray(do)]
            x, y = coords[:, 0], coords[:, 1]
            xr = _pad_range(float(np.nanmin(x)), float(np.nanmax(x)), 0.05)
            yr = _pad_range(float(np.nanmin(y)), float(np.nanmax(y)), 0.05)
            xb = _q_u16(_normalise_range(x, xr[0], xr[1]))
            yb = _q_u16(_normalise_range(y, yr[0], yr[1]))
            lbl = key[2:] if key[:2].lower() == "x_" else key
            # ship the channels as memoryviews (matches the viewport binary path the
            # front-end already decodes; raw bytes can arrive in a shape the u16
            # decoder mis-reads -> a silent no-op).
            self.send(
                {"type": "morph", "duration": int(duration), "n": int(x.size),
                 "x_min": xr[0], "x_max": xr[1], "y_min": yr[0], "y_max": yr[1],
                 "xlab": "%s 1" % lbl, "ylab": "%s 2" % lbl},
                buffers=[memoryview(np.ascontiguousarray(xb)),
                         memoryview(np.ascontiguousarray(yb))],
            )
            # NB: returns None on purpose — returning self would make Jupyter render a
            # SECOND copy of the widget in the cell output.

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

        def _pick_expr_matrix(self, ad, layer, use_raw):
            """Choose ``(layer, use_raw)`` for ``sc.tl.rank_genes_groups`` so the
            logfoldchanges come out FINITE.

            logfoldchanges need NON-NEGATIVE (log-normalised) values; running on a
            SCALED / z-scored ``.X`` (which has negatives) makes scanpy emit NaN
            logFC. When the caller didn't pin a matrix and ``.X`` looks scaled,
            route to log-norm expression for the SAME genes (``adata.raw``
            restricted to ``var_names``) so logFC is finite without dragging in
            non-HVG genes; fall back to a log-norm layer, else warn. Shared by
            :meth:`diff_expression` and :meth:`diff_expression_by`.

            ``ad`` is always a throwaway copy, so for the ``adata.raw`` route we
            overwrite ``ad.X`` in place rather than stashing a synthetic layer —
            that keeps scanpy's saved ``uns[...]['params']`` clean (``layer: None``,
            ``use_raw: False``) instead of leaking a cryptic ``"_rs_lognorm"`` name.
            """
            import warnings

            import numpy as np

            _layer, _ur = layer, use_raw
            if _ur is None and _layer is None:
                def _neg(x):
                    try:
                        import scipy.sparse as _sp
                        if _sp.issparse(x):
                            return x.data.size > 0 and float(x.data.min()) < -1e-6
                    except Exception:
                        pass
                    a = np.asarray(x)
                    return a.size > 0 and float(np.nanmin(a)) < -1e-6

                if not _neg(ad.X):
                    _ur = False                                  # .X already log-norm
                elif ad.raw is not None and set(ad.var_names).issubset(set(ad.raw.var_names)):
                    ad.X = ad.raw[:, list(ad.var_names)].X       # finite logFC, same genes
                    _layer, _ur = None, False                    # clean params, no leaked layer
                else:
                    _cand = next((ln for ln in ("lognorm", "log1p", "logcounts", "data")
                                  if ln in ad.layers and not _neg(ad.layers[ln])), None)
                    if _cand is not None:
                        _layer, _ur = _cand, False
                    else:
                        _ur = False
                        warnings.warn(
                            "diff_expression: adata.X looks scaled (z-scored), so scanpy's "
                            "logfoldchanges will be NaN (the scores / p-values are still "
                            "valid). For finite logFC, set `adata.raw = adata` on the "
                            "log-normalised data before scaling, or pass a log-norm `layer=`.",
                            stacklevel=2,
                        )
            return _layer, _ur

        def _ttest_de(self, X, a_mask, b_mask, var_names, n):
            """No-scanpy fallback: Welch t-test of group A vs B over ``X`` rows.
            Returns the top-``n`` genes by |stat| as a tidy DataFrame."""
            import numpy as np
            import pandas as pd

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
                "gene": np.asarray(var_names),
                "logFC": lfc, "stat": stat, "pval": pval,
                "mean_A": ma, "mean_B": mb,
            })
            res = res.reindex(res["stat"].abs().sort_values(ascending=False).index)
            return res.head(n).reset_index(drop=True)

        def diff_expression(self, group_a=None, group_b=None, n=10, layer=None,
                            method="wilcoxon", key_added=None, use_raw=None):
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

                    # Pick the expression matrix so logFC is finite (see helper).
                    _layer, _ur = self._pick_expr_matrix(ad, layer, use_raw)

                    sc.tl.rank_genes_groups(
                        ad, "_rs_grp", groups=["A"], reference=ref,
                        method=method, layer=_layer, n_genes=n, key_added=_key,
                        use_raw=_ur,
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
            return _save(self._ttest_de(X, a_mask, b_mask, src.var_names, n))

        def diff_expression_by(self, by, group_a=None, group_b=None, selection=None,
                               n=10, layer=None, method="wilcoxon", key_added=None,
                               use_raw=None, min_cells=2):
            """Differential expression BETWEEN the levels of an ``obs`` column,
            restricted to the lasso selection.

            Lasso a group of cells, then split them by ``by`` (an ``obs`` column
            such as ``"time"`` or ``"condition"``) and compare its levels:

            * ``group_a`` **and** ``group_b`` given -> a single A-vs-B comparison
              (e.g. ``group_a="D30", group_b="Y1"``); returns one DataFrame.
            * ``group_a`` only -> that level vs the pooled rest of the selection;
              returns one DataFrame.
            * **neither** -> every present level vs the rest, run as ONE scanpy
              ``rank_genes_groups`` call (one-vs-rest, the scanpy idiom); returns a
              ``dict`` ``{level: df}``.

            Whichever the mode, the result is saved to ``adata.uns`` as a SINGLE
            scanpy-native entry — with the REAL ``by`` name and REAL level names in
            ``params`` / columns (no synthetic ``_rs_by`` / ``_rs_lognorm``), so
            ``sc.pl.rank_genes_groups(adata)`` works straight after. Default key
            ``"rank_genes_groups"``; pass ``key_added`` to choose the key, or
            ``key_added=False`` to skip saving.

            Cells default to the current lasso ``selection`` (pass ``selection=`` to
            override, e.g. integer positions / obs_names / a boolean mask); if
            nothing is selected it falls back to **all** cells. In the all-levels
            mode, levels with fewer than ``min_cells`` cells in the selection are
            skipped (with a warning); an explicit ``group_a`` / ``group_b`` is always
            honoured. Uses ``sc.tl.rank_genes_groups`` when scanpy is installed
            (same finite-logFC matrix routing as :meth:`diff_expression`), else a
            Welch t-test. AnnData/MuData only.
            """
            import warnings

            import numpy as np
            import pandas as pd

            src = getattr(self, "_source", None)
            if src is None or not hasattr(src, "X"):
                raise TypeError("diff_expression_by() needs an AnnData/MuData with .X.")
            frame = src.obs if hasattr(src, "obs") else src
            if not hasattr(frame, "columns") or by not in frame.columns:
                cols = list(getattr(frame, "columns", []))
                raise KeyError(
                    f"'{by}' is not an obs column. Available e.g.: {cols[:8]}"
                )

            # Resolve the cell set: explicit selection > lasso > all cells.
            if selection is not None:
                sel = self._resolve_orig_indices(selection)
            else:
                sel = list(self.selection)
                if not sel:
                    sel = list(range(src.n_obs))   # nothing lassoed -> whole dataset

            col = frame[by]
            grp = col.iloc[sel]
            # Levels actually present (drop NaN); keep categorical order if any.
            if hasattr(grp, "cat"):
                present = [lv for lv in grp.cat.categories if (grp == lv).any()]
            else:
                present = list(pd.unique(grp.dropna()))
            counts = {lv: int((grp == lv).sum()) for lv in present}
            usable = [lv for lv in present if counts[lv] >= min_cells]
            dropped = [lv for lv in present if counts[lv] < min_cells]

            def _coerce(val, role):
                # accept the level as given or stringified (categoricals print as str)
                for lv in present:
                    if val == lv or str(val) == str(lv):
                        return lv
                raise ValueError(
                    f"{role}={val!r} is not a level of '{by}' present in the "
                    f"selection. Present levels: {present}"
                )

            ga = None if group_a is None else _coerce(group_a, "group_a")
            gb = None if group_b is None else _coerce(group_b, "group_b")

            # Decide the comparison: which level(s) are tested, the reference,
            # which levels' cells to keep, and whether the result is a dict (multi).
            if ga is not None and gb is not None:
                groups, reference, keep_levels, multi = [ga], gb, [ga, gb], False
            elif ga is not None:
                if len(present) < 2:
                    raise ValueError(
                        f"group_a={group_a!r} vs rest needs >= 2 levels of '{by}' in "
                        f"the selection; only {present} present."
                    )
                groups, reference, keep_levels, multi = [ga], "rest", list(present), False
            else:
                if dropped:
                    warnings.warn(
                        f"diff_expression_by: skipping {len(dropped)} level(s) of "
                        f"'{by}' with < {min_cells} cells in the selection: "
                        f"{[f'{lv} (n={counts[lv]})' for lv in dropped]}",
                        stacklevel=2,
                    )
                if len(usable) < 2:
                    raise ValueError(
                        f"need >= 2 levels of '{by}' with >= {min_cells} cells in the "
                        f"selection; got {usable or present}."
                    )
                groups, reference, keep_levels, multi = list(usable), "rest", list(usable), True

            # Cells in the comparison = selected cells whose level is kept.
            keep_str = {str(lv) for lv in keep_levels}
            keep = [i for i in sel if str(col.iloc[i]) in keep_str]
            g_str = [str(lv) for lv in groups]
            ref_str = reference if reference == "rest" else str(reference)

            save_to = (None if key_added is False or not hasattr(src, "uns")
                       else (key_added or "rank_genes_groups"))

            # --- scanpy path: ONE rank_genes_groups call on the real `by` column ---
            if type(src).__name__ == "AnnData":
                try:
                    import scanpy as sc

                    ad = src[keep].copy()
                    # rebuild `by` as a clean categorical of just the kept levels,
                    # keeping the REAL column name so uns[...]['params'] reads true.
                    ad.obs[by] = pd.Categorical(
                        ad.obs[by].astype(str),
                        categories=[str(lv) for lv in keep_levels],
                    )
                    _layer, _ur = self._pick_expr_matrix(ad, layer, use_raw)
                    _key = save_to or "rank_genes_groups"
                    sc.tl.rank_genes_groups(
                        ad, by, groups=g_str, reference=ref_str,
                        method=method, layer=_layer, n_genes=n, key_added=_key,
                        use_raw=_ur,
                    )
                    if save_to is not None:
                        src.uns[save_to] = ad.uns.get(_key)   # one clean native entry
                    frames = {g: sc.get.rank_genes_groups_df(ad, group=g, key=_key)
                                   .head(n).reset_index(drop=True) for g in g_str}
                    return frames if multi else frames[g_str[0]]
                except ImportError:
                    pass  # fall through to the Welch t-test

            # --- no-scanpy fallback: Welch t-test per group vs its reference ---
            X = src.layers[layer] if layer else src.X

            def _mask(level_strs):
                m = np.zeros(src.n_obs, dtype=bool)
                for i in keep:
                    if str(col.iloc[i]) in level_strs:
                        m[i] = True
                return m

            frames = {}
            for g in g_str:
                a_mask = _mask({g})
                b_mask = (_mask({s for s in keep_str if s != g})
                          if ref_str == "rest" else _mask({ref_str}))
                if a_mask.any() and b_mask.any():
                    frames[g] = self._ttest_de(X, a_mask, b_mask, src.var_names, n)
            if not frames:
                raise ValueError("every comparison had an empty side after filtering.")
            if save_to is not None:
                src.uns[save_to] = frames if multi else frames[g_str[0]]
            return frames if multi else frames[g_str[0]]

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
        # Theme mode for the LIVE widget: "light" (default white card), "dark",
        # or "auto" (match a dark host). Rides its own trait, NOT _spec, so the
        # R-parity payload stays byte-identical.
        _theme = traitlets.Unicode("light").tag(sync=True)
        _selection = traitlets.List(trait=traitlets.Int()).tag(sync=True)
        _filtered = traitlets.List(trait=traitlets.Int()).tag(sync=True)
        _filtered_on = traitlets.Bool(False).tag(sync=True)
        _camera = traitlets.List(trait=traitlets.Float()).tag(sync=True)
        # Two barriers let a synchronous w.selection/w.filtered read wait for the
        # latest interaction WITHOUT over-blocking when idle:
        #   _sel_gen  -> bumped by the front-end on EVERY interaction (light + heavy).
        #                The observer acks it immediately and reliably (it fires on
        #                the counter itself, which always changes — unlike a
        #                _selection observer that's silent on an unchanged value).
        #                This is the liveness signal: idle => _applied_gen == _sel_gen.
        #   _work_req -> bumped BEFORE the heavy async work (a 10M-point lasso /
        #                legend filter / deselect / reset). The kernel sets _work_done
        #                only AFTER that work finishes (_ack_work), so the pump blocks
        #                for completion, not just delivery.
        _sel_gen = traitlets.Int(0).tag(sync=True)
        _applied_gen = 0
        _work_req = traitlets.Int(0).tag(sync=True)
        _work_done = 0
        # set by the legend colorpicker so w.colors reflects an in-plot recolour
        _legend_colors = traitlets.List(trait=traitlets.Unicode()).tag(sync=True)
        _legend_names = traitlets.List(trait=traitlets.Unicode()).tag(sync=True)

        @traitlets.observe("_sel_gen")
        def _ack_gen(self, change):
            self._applied_gen = int(change["new"] or 0)

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
