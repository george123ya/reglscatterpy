"""Widget-level tests for scatterplot(): basis=, color_by list grid,
gene-aware size_by/opacity_by, and clear error messages.

Needs the widget stack (anndata + ipywidgets + anywidget); skipped otherwise.
"""

import numpy as np
import pandas as pd
import pytest

import reglscatterpy as rs


def _toy_anndata(n=80, g=12):
    ad = pytest.importorskip("anndata")
    rng = np.random.RandomState(2)
    X = rng.poisson(2, size=(n, g)).astype(float)
    var = pd.DataFrame(index=[f"Gene{i}" for i in range(g)])
    obs = pd.DataFrame(
        {
            "celltype": pd.Categorical(rng.choice(["T", "B", "NK"], n)),
            "score": rng.rand(n),
        },
        index=[f"cell{i}" for i in range(n)],
    )
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm["X_umap"] = rng.randn(n, 2)
    return adata


@pytest.fixture
def adata():
    pytest.importorskip("anywidget")
    pytest.importorskip("ipywidgets")
    return _toy_anndata()


# --- Change 1: basis= ------------------------------------------------------ #
def test_basis_equals_x(adata):
    wa = rs.scatterplot(adata, basis="umap", color_by="celltype", show=False)
    wb = rs.scatterplot(adata, x="umap", color_by="celltype", show=False)
    assert wa._spec["x"] == wb._spec["x"]          # identical encoded coords
    assert wa._spec["xlab"] == "UMAP 1"


def test_basis_on_dataframe_raises():
    pytest.importorskip("anywidget")
    df = pd.DataFrame({"u1": [0.0, 1.0], "u2": [1.0, 0.0], "ct": ["A", "B"]})
    with pytest.raises(ValueError, match="basis="):
        rs.scatterplot(df, x="u1", y="u2", basis="umap", show=False)


# --- Change 2: color_by list -> linked grid -------------------------------- #
def test_color_by_list_makes_linked_grid(adata):
    w = rs.scatterplot(adata, basis="umap", color_by=["celltype", "Gene3"], show=False)
    assert type(w).__name__ == "GridBox"
    assert len(w.children) == 2
    titles = {c._spec["title"] for c in w.children}
    assert titles == {"celltype", "Gene3"}
    assert all(isinstance(c._spec["syncPlots"], list) for c in w.children)


def test_raw_vector_color_is_not_a_grid(adata):
    col = np.repeat(["a", "b"], adata.n_obs // 2)   # ndarray -> raw vector
    w = rs.scatterplot(adata, basis="umap", color_by=col, show=False)
    assert type(w).__name__ != "GridBox"


# --- Change 3: gene-aware size_by / opacity_by ----------------------------- #
def test_size_by_gene_resolves(adata):
    w = rs.scatterplot(adata, basis="umap", size_by="Gene3", show=False)
    assert w._spec["sizeBy"] is True


def test_opacity_by_gene_resolves(adata):
    w = rs.scatterplot(adata, basis="umap", opacity_by="Gene3", show=False)
    assert w._spec["opacityBy"] is True


def test_size_by_obs_column_still_works(adata):
    w = rs.scatterplot(adata, basis="umap", size_by="score", show=False)
    assert w._spec["sizeBy"] is True


# --- Change 4: clear error messages ---------------------------------------- #
def test_bad_color_name_suggests(adata):
    with pytest.raises(KeyError) as e:
        rs.scatterplot(adata, basis="umap", color_by="celltpye", show=False)  # typo
    msg = str(e.value)
    assert "Did you mean" in msg or "Available" in msg


def test_bad_size_name_lists_available(adata):
    with pytest.raises(KeyError) as e:
        rs.scatterplot(adata, basis="umap", size_by="Gene999", show=False)
    assert "size_by" in str(e.value)


def test_bad_toolbar_lists_choices(adata):
    with pytest.raises(ValueError, match="toolbar"):
        rs.scatterplot(adata, basis="umap", toolbar="sideways", show=False)


def test_bad_palette_suggests(adata):
    with pytest.raises(KeyError) as e:
        rs.scatterplot(adata, basis="umap", continuous_palette="viridiss", show=False)
    assert "continuous_palette" in str(e.value)


def test_wrong_length_raw_color_raises(adata):
    with pytest.raises(ValueError, match="length"):
        rs.scatterplot(adata, basis="umap", color_by=np.array([1.0, 2.0, 3.0]), show=False)
