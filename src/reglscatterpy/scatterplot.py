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


def _is_name_list(spec) -> bool:
    """True for a list/tuple of strings (= panel names for a color-by grid),
    as opposed to a raw per-point vector (which must be an ndarray/Series)."""
    return (
        isinstance(spec, (list, tuple))
        and len(spec) > 0
        and all(isinstance(e, str) for e in spec)
    )


def _resolve_numeric(spec, data, layer=None, param="size_by"):
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
                vec, _ = _resolve_anndata_vec(data, spec, layer)
                return np.asarray(vec, dtype="float64")
            if _is_mudata(data) and ":" in spec:
                mod, key = spec.split(":", 1)
                if mod in data.mod:
                    vec, _ = _resolve_anndata_vec(data.mod[mod], key, layer)
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
    basis: Optional[Union[str, int]] = None,
    x: Optional[Union[str, int]] = None,
    y: Optional[Union[str, int]] = None,
    color_by: ColorSpec = None,
    group_by: ColorSpec = None,
    layer: Optional[str] = None,
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
    """
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
                layer=layer, dims=dims, table=table, point_size=point_size,
                opacity=opacity, point_color=point_color, size_by=size_by,
                opacity_by=opacity_by, tooltip_by=tooltip_by,
                pixel_ratio=pixel_ratio, categorical_palette=categorical_palette,
                continuous_palette=continuous_palette, custom_palette=custom_palette,
                custom_colors=custom_colors, vmin=vmin, vmax=vmax,
                center_zero=center_zero, title=(title or name), xlab=xlab, ylab=ylab,
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
        return compose(panels)

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
        layer=layer, dims=dims, table=table,
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

    size_values = _resolve_numeric(size_by, data, layer, param="size_by")
    opacity_values = _resolve_numeric(opacity_by, data, layer, param="opacity_by")
    tooltip_cols = _resolve_cols(tooltip_by, data)
    for pname, vec in (("size_by", size_values), ("opacity_by", opacity_values)):
        if vec is not None and len(vec) != n:
            raise ValueError(
                f"{pname} has length {len(vec)} but the data has {n} points."
            )

    spec = build_payload(
        pd_data,
        point_size=point_size, opacity=opacity, point_color=point_color,
        size_values=size_values, opacity_values=opacity_values,
        tooltip_cols=tooltip_cols,
        pixel_ratio=pixel_ratio,
        categorical_palette=categorical_palette, continuous_palette=continuous_palette,
        custom_palette=custom_palette, custom_colors=custom_colors,
        vmin=vmin, vmax=vmax, center_zero=center_zero,
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
    )

    w = int(width) if width else 0   # 0 => responsive (100%)
    if interactive:
        # Live, kernel-linked widget: w.selection round-trips (needs a kernel).
        from ._widget import ReglScatter

        widget = ReglScatter()
        widget._height = int(height)
        widget._width = w
        widget._source = data   # so w.annotate(...) can write back to obs/colData
        widget._spec = spec
        return widget

    # Default: a static, self-contained plot — renders with no comm and reopens
    # with no kernel (like a plotly figure), and is NOT an ipywidget, so nothing
    # is written to the notebook's widget-state.
    from ._widget import StaticPlot

    plot = StaticPlot(spec=spec, source=data, height=int(height), width=w)
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
