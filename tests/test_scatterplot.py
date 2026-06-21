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


# --- Static-by-default render (no widget) ---------------------------------- #
def test_default_render_is_static_iframe(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype", show=False)
    bundle = w._repr_mimebundle_()
    data = bundle[0] if isinstance(bundle, tuple) else bundle
    assert "text/html" in data
    assert "application/vnd.jupyter.widget-view+json" not in data  # no live comm
    html = data["text/html"]
    assert "<iframe" in html and "srcdoc=" in html
    # srcdoc is HTML-escaped, so quotes become entities; check the unquoted stem
    assert "DecompressionStream(" in html  # self-contained bundle inlined


def test_interactive_render_is_live_widget(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype",
                       interactive=True, show=False)
    bundle = w._repr_mimebundle_()
    data = bundle[0] if isinstance(bundle, tuple) else bundle
    assert "application/vnd.jupyter.widget-view+json" in data  # live comm widget


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


# --- na_color + groups ----------------------------------------------------- #
def test_groups_greys_unlisted_categories(adata):
    w = rs.scatterplot(adata, basis="umap", color="celltype",
                       groups=["T"], na_color="#eeeeee", show=False)
    leg = w._spec["legend"]
    cols = dict(zip(leg["names"], leg["colors"]))
    assert cols["T"] != "#eeeeee"                 # kept its palette colour
    assert cols["B"] == "#eeeeee" and cols["NK"] == "#eeeeee"   # greyed out


# --- use_raw (scanpy default: gene from .raw when present) ----------------- #
def test_use_raw_default_reads_raw_not_scaled_X():
    ad = pytest.importorskip("anndata")
    rng = np.random.RandomState(7)
    n, g = 60, 6
    names = [f"Gene{i}" for i in range(g)]
    lognorm = rng.gamma(2, 1.0, (n, g))                 # .raw: 0..~6
    adata = ad.AnnData(
        X=(lognorm - lognorm.mean(0)) / (lognorm.std(0) + 1e-9),  # .X: z-scored (neg)
        obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=names),
    )
    adata.raw = ad.AnnData(X=lognorm, var=pd.DataFrame(index=names))
    adata.obsm["X_umap"] = rng.randn(n, 2)

    raw_plot = rs.scatterplot(adata, basis="umap", color="Gene0", show=False)   # default
    assert raw_plot._spec["legend"]["minVal"] >= 0          # log-normalised, non-negative
    x_plot = rs.scatterplot(adata, basis="umap", color="Gene0", use_raw=False, show=False)
    assert x_plot._spec["legend"]["minVal"] < 0             # scaled .X has negatives


# --- scanpy-aligned argument names ----------------------------------------- #
def test_color_alias_matches_color_by(adata):
    a = rs.scatterplot(adata, basis="umap", color="celltype", show=False)
    b = rs.scatterplot(adata, basis="umap", color_by="celltype", show=False)
    assert a._spec["x"] == b._spec["x"]
    assert a._spec.get("colorVar") == b._spec.get("colorVar") == "celltype"


def test_size_scalar_vs_name(adata):
    big = rs.scatterplot(adata, basis="umap", size=12, show=False)
    assert big._spec["options"]["size"] == 12
    per = rs.scatterplot(adata, basis="umap", size="Gene3", show=False)  # gene -> per-point
    assert per._spec["sizeBy"] is True


def test_components_is_one_based(adata):
    adata.obsm["X_pca"] = __import__("numpy").random.RandomState(3).randn(adata.n_obs, 4)
    a = rs.scatterplot(adata, basis="pca", components=(2, 3), show=False)
    b = rs.scatterplot(adata, basis="pca", dims=(1, 2), show=False)   # 0-based equivalent
    assert a._spec["x"] == b._spec["x"]
    assert a._spec["xlab"] == "PCA 2"


def test_cmap_and_palette_aliases(adata):
    w = rs.scatterplot(adata, basis="umap", color="Gene3", cmap="magma", show=False)
    assert w._spec is not None       # resolves without error (magma is valid)
    w2 = rs.scatterplot(adata, basis="umap", color="celltype", palette="Dark2", show=False)
    assert w2._spec is not None


def test_ncols_in_grid(adata):
    g = rs.scatterplot(adata, basis="umap", color=["celltype", "Gene3", "Gene4"],
                       ncols=3, show=False)
    cols = g.layout.grid_template_columns
    assert "repeat(3" in cols


def test_save_html(tmp_path, adata):
    out = tmp_path / "p.html"
    rs.scatterplot(adata, basis="umap", color="celltype", save=str(out), show=False)
    assert out.exists() and "<iframe" not in out.read_text()[:20]  # full standalone page


def test_save_bad_extension_raises(tmp_path, adata):
    with pytest.raises(ValueError, match="only '.html'"):
        rs.scatterplot(adata, basis="umap", color="celltype",
                       save=str(tmp_path / "p.png"), show=False)


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
