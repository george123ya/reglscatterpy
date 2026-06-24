"""reglscatterpy - interactive WebGL scatterplots for single-cell data in Python.

The Python companion to the R package ``reglScatterplotR``. Both wrap the
``regl-scatterplot`` WebGL engine; this package adds AnnData / MuData /
SpatialData awareness so you can go from a scanpy object to an interactive
million-point plot in one call::

    import scanpy as sc
    import reglscatterpy as rs

    adata = sc.datasets.pbmc3k_processed()
    rs.scatterplot(adata, x="X_umap", color_by="louvain")
    rs.scatterplot(adata, x="X_umap", color_by="CST3")  # a gene
"""

from . import _export
from ._compose import compose
from ._export import record_html, save_html, save_notebook_html
from ._extract import PlotData, extract
from .scatterplot import scatterplot

__all__ = [
    "scatterplot",
    "compose",
    "save_html",
    "save_notebook_html",
    "record_html",
    "extract",
    "PlotData",
    "__version__",
]
__version__ = "0.6.20"
