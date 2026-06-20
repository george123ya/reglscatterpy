"""Backend-agnostic data extraction.

This module turns any supported input object (AnnData, MuData, SpatialData,
pandas DataFrame, numpy array) into a small normalised :class:`PlotData` of
plain numpy/pandas pieces. It deliberately knows nothing about *rendering* -
so the exact same extraction feeds the current jupyter-scatter backend and any
future standalone widget, keeping behaviour identical across the two.

The single-cell resolution rules mirror the R package (``reglScatterplotR``):

* ``x`` selects an embedding (``obsm`` key); ``"umap"`` is auto-prefixed to
  the scanpy-style ``"X_umap"`` and matched case-insensitively, falling back
  to the first of umap / tsne / pca present.
* ``color_by`` / ``group_by`` resolve against ``obs`` columns first, then
  against ``var_names`` (a feature), in which case ``X`` - or ``layer`` /
  ``raw`` - is read.
"""

from __future__ import annotations

from dataclasses import dataclass
from difflib import get_close_matches
from typing import Any, Optional, Sequence, Union

import numpy as np

__all__ = ["PlotData", "extract"]

ColorSpec = Union[str, Sequence[Any], np.ndarray, None]


def _not_found(param: str, value: Any, available, *, searched: str = "") -> KeyError:
    """Build a uniform, actionable 'name not found' error.

    States the parameter, the bad value, what was searched, a 'did you mean'
    suggestion, and a (truncated) list of valid names.
    """
    avail = [str(a) for a in available]
    msg = f"{param}={value!r} not found"
    if searched:
        msg += f" (searched {searched})"
    msg += "."
    hit = get_close_matches(str(value), avail, n=1)
    if hit:
        msg += f" Did you mean {hit[0]!r}?"
    if avail:
        shown = avail[:15]
        more = "" if len(avail) <= 15 else f", … (+{len(avail) - 15} more)"
        msg += f" Available: {shown}{more}"
    return KeyError(msg)


@dataclass
class PlotData:
    """Normalised, render-agnostic plot inputs."""

    x: np.ndarray
    y: np.ndarray
    color: Optional[np.ndarray] = None
    group: Optional[np.ndarray] = None
    color_name: Optional[str] = None
    group_name: Optional[str] = None
    xlab: str = "X"
    ylab: str = "Y"

    @property
    def n(self) -> int:
        return int(self.x.shape[0])


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def _is_anndata(obj: Any) -> bool:
    return type(obj).__name__ == "AnnData" and hasattr(obj, "obsm")


def _is_mudata(obj: Any) -> bool:
    return type(obj).__name__ == "MuData" and hasattr(obj, "mod")


def _is_spatialdata(obj: Any) -> bool:
    return type(obj).__name__ == "SpatialData" and hasattr(obj, "tables")


def _densify(col: Any) -> np.ndarray:
    """Flatten a possibly-sparse single column to a 1-D numpy array."""
    if hasattr(col, "toarray"):  # scipy.sparse
        col = col.toarray()
    return np.asarray(col).ravel()


def _resolve_basis(adata: Any, x: Optional[str]) -> str:
    """Map a user ``x`` to a concrete ``obsm`` key."""
    keys = list(adata.obsm.keys())
    lower = {k.lower(): k for k in keys}

    if x is None:
        for pref in ("x_umap", "x_tsne", "x_pca", "spatial"):
            if pref in lower:
                return lower[pref]
        if keys:
            return keys[0]
        raise ValueError(
            "No embeddings found in .obsm; run e.g. sc.tl.umap() first."
        )

    cand = [x, f"X_{x}"]
    for c in cand:
        if c in adata.obsm:
            return c
        if c.lower() in lower:
            return lower[c.lower()]
    raise _not_found("basis/x", x, keys, searched=".obsm (also tried 'X_' prefix)")


def _resolve_anndata_vec(adata: Any, key: ColorSpec, layer: Optional[str]):
    """Return (vector, name) for a color/group spec against an AnnData."""
    if key is None:
        return None, None
    if not isinstance(key, str):
        return np.asarray(key), None  # raw vector passed through

    if key in adata.obs.columns:
        return adata.obs[key].to_numpy(), key

    var_names = adata.raw.var_names if (layer == "raw" and adata.raw) else adata.var_names
    if key in var_names:
        if layer == "raw" and adata.raw is not None:
            mat = adata.raw[:, key].X
        elif layer is not None:
            mat = adata[:, key].layers[layer]
        else:
            mat = adata[:, key].X
        return _densify(mat), key

    raise _not_found(
        "color_by/group_by", key,
        list(adata.obs.columns) + list(var_names),
        searched=".obs columns and .var_names",
    )


# --------------------------------------------------------------------------- #
# per-type extractors
# --------------------------------------------------------------------------- #
def _from_anndata(
    adata, x, y, color_by, group_by, layer, dims
) -> PlotData:
    basis = _resolve_basis(adata, x)
    coords = np.asarray(adata.obsm[basis])
    d0, d1 = (dims or (0, 1))
    if coords.shape[1] <= max(d0, d1):
        raise ValueError(
            f"Embedding {basis!r} has {coords.shape[1]} dims; "
            f"cannot take dims {(d0, d1)}."
        )
    color, color_name = _resolve_anndata_vec(adata, color_by, layer)
    group, group_name = _resolve_anndata_vec(adata, group_by, layer)
    label = basis[2:] if basis.startswith("X_") else basis
    return PlotData(
        x=coords[:, d0],
        y=coords[:, d1],
        color=color,
        group=group,
        color_name=color_name,
        group_name=group_name,
        xlab=f"{label.upper()} {d0 + 1}",
        ylab=f"{label.upper()} {d1 + 1}",
    )


def _from_mudata(mdata, x, y, color_by, group_by, layer, dims) -> PlotData:
    # Global embeddings live on the MuData object itself; per-modality features
    # are addressed as "modality:feature" (e.g. "rna:CD3D").
    def split_mod(spec):
        if isinstance(spec, str) and ":" in spec:
            mod, key = spec.split(":", 1)
            if mod in mdata.mod:
                return _resolve_anndata_vec(mdata.mod[mod], key, layer)
        return None

    # Embedding: support a global MuData embedding ("X_umap") or a per-modality
    # one ("rna:X_umap", stored on mdata['rna'].obsm).
    if isinstance(x, str) and ":" in x and x.split(":", 1)[0] in mdata.mod:
        mod, key = x.split(":", 1)
        sub = mdata.mod[mod]
        basis = _resolve_basis(sub, key)
        coords = np.asarray(sub.obsm[basis])
    else:
        basis = _resolve_basis(mdata, x)
        coords = np.asarray(mdata.obsm[basis])
    d0, d1 = (dims or (0, 1))

    def resolve_global(spec):
        if spec is None:
            return None, None
        sub = split_mod(spec)
        if sub is not None:
            return sub
        if isinstance(spec, str) and spec in mdata.obs.columns:
            return mdata.obs[spec].to_numpy(), spec
        if not isinstance(spec, str):
            return np.asarray(spec), None
        raise _not_found(
            "color_by/group_by", spec, list(mdata.obs.columns),
            searched=".obs (for a feature use 'modality:feature', e.g. 'rna:CD3D')",
        )

    color, color_name = resolve_global(color_by)
    group, group_name = resolve_global(group_by)
    label = basis[2:] if basis.startswith("X_") else basis
    return PlotData(
        coords[:, d0], coords[:, d1], color, group,
        color_name, group_name,
        f"{label.upper()} {d0 + 1}", f"{label.upper()} {d1 + 1}",
    )


def _from_spatialdata(
    sdata, x, y, color_by, group_by, layer, dims, table
) -> PlotData:
    tables = list(sdata.tables.keys())
    if not tables:
        raise ValueError("SpatialData object has no annotation tables.")
    tbl = table or tables[0]
    if tbl not in sdata.tables:
        raise _not_found("table", tbl, tables, searched=".tables")
    adata = sdata.tables[tbl]
    # x defaults to the spatial coordinates for a SpatialData input.
    return _from_anndata(adata, x or "spatial", y, color_by, group_by,
                         layer, dims)


def _from_dataframe(df, x, y, color_by, group_by) -> PlotData:
    if x is None or y is None:
        raise ValueError("For a DataFrame, pass x= and y= column names.")
    cols = list(df.columns)
    for param, name in (("x", x), ("y", y)):
        if name not in df.columns:
            raise _not_found(param, name, cols, searched="DataFrame columns")

    def col(param, spec):
        if spec is None:
            return None, None
        if isinstance(spec, str):
            if spec not in df.columns:
                raise _not_found(param, spec, cols, searched="DataFrame columns")
            return df[spec].to_numpy(), spec
        return np.asarray(spec), None

    color, color_name = col("color_by", color_by)
    group, group_name = col("group_by", group_by)
    return PlotData(
        df[x].to_numpy(), df[y].to_numpy(), color, group,
        color_name, group_name, str(x), str(y),
    )


def _from_array(arr, x, y, color_by, group_by) -> PlotData:
    arr = np.asarray(arr)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError("A coordinate array needs shape (n, >=2).")
    i, j = (0 if x is None else int(x)), (1 if y is None else int(y))
    color = None if color_by is None else np.asarray(color_by)
    group = None if group_by is None else np.asarray(group_by)
    return PlotData(arr[:, i], arr[:, j], color, group)


# --------------------------------------------------------------------------- #
# public entry point
# --------------------------------------------------------------------------- #
def extract(
    data: Any,
    *,
    x: Optional[Union[str, int]] = None,
    y: Optional[Union[str, int]] = None,
    color_by: ColorSpec = None,
    group_by: ColorSpec = None,
    layer: Optional[str] = None,
    dims: Optional[tuple] = None,
    table: Optional[str] = None,
) -> PlotData:
    """Normalise any supported input into :class:`PlotData`."""
    if _is_anndata(data):
        result = _from_anndata(data, x, y, color_by, group_by, layer, dims)
    elif _is_mudata(data):
        result = _from_mudata(data, x, y, color_by, group_by, layer, dims)
    elif _is_spatialdata(data):
        result = _from_spatialdata(
            data, x, y, color_by, group_by, layer, dims, table
        )
    elif type(data).__name__ == "DataFrame":
        result = _from_dataframe(data, x, y, color_by, group_by)
    elif isinstance(data, np.ndarray) or hasattr(data, "__array__"):
        result = _from_array(data, x, y, color_by, group_by)
    else:
        raise TypeError(
            f"Unsupported input type {type(data).__name__!r}. Pass an AnnData, "
            "MuData, SpatialData, pandas DataFrame, or numpy array."
        )

    # Catch a raw color/group vector whose length doesn't match the points
    # (otherwise it silently mis-renders).
    n = result.n
    for param, vec in (("color_by", result.color), ("group_by", result.group)):
        if vec is not None and len(vec) != n:
            raise ValueError(
                f"{param} has length {len(vec)} but the data has {n} points."
            )
    return result
