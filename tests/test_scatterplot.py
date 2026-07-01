"""Widget-level tests for scatterplot(): basis=, color_by list grid,
gene-aware size_by/opacity_by, and clear error messages.

Needs the widget stack (anndata + ipywidgets + anywidget); skipped otherwise.
"""

import importlib.util

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
    adata.obsm["X_tsne"] = rng.randn(n, 2)
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


# --- Render mode: interactive by default, static when opted out ------------ #
def test_default_render_is_interactive(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype", show=False)
    from reglscatterpy import _widget
    assert _widget.is_live_widget(w)   # interactive=True is the default now


def test_static_render_is_iframe(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype", interactive=False, show=False)
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
    # default (static) -> an HTML iframe-grid that renders without ipywidgets
    h = rs.scatterplot(adata, basis="umap", color_by=["celltype", "Gene3"], interactive=False, show=False)
    assert type(h).__name__ == "_HtmlGrid"
    assert len(h.panels) == 2 and "text/html" in h._repr_mimebundle_()
    # interactive=True -> a linked live GridBox
    w = rs.scatterplot(adata, basis="umap", color_by=["celltype", "Gene3"],
                       interactive=True, show=False)
    assert "GridBox" in [c.__name__ for c in type(w).__mro__]
    assert len(w.children) == 2
    titles = {c._spec["title"] for c in w.children}
    assert titles == {"celltype", "Gene3"}
    assert all(isinstance(c._spec["syncPlots"], list) for c in w.children)


def test_raw_vector_color_is_not_a_grid(adata):
    col = np.repeat(["a", "b"], adata.n_obs // 2)   # ndarray -> raw vector
    w = rs.scatterplot(adata, basis="umap", color_by=col, show=False)
    assert type(w).__name__ != "GridBox"


def test_basis_list_makes_grid(adata):
    h = rs.scatterplot(adata, basis=["umap", "tsne"], color_by="celltype", interactive=False, show=False)
    assert type(h).__name__ == "_HtmlGrid" and len(h.panels) == 2
    assert {p._spec["title"] for p in h.panels} == {"umap", "tsne"}


def test_basis_color_cross_product(adata):
    w = rs.scatterplot(adata, basis=["umap", "tsne"], color_by=["celltype", "Gene3"],
                       interactive=True, show=False)
    assert "GridBox" in [c.__name__ for c in type(w).__mro__]
    assert len(w.children) == 4
    assert {c._spec["title"] for c in w.children} == {
        "umap · celltype", "umap · Gene3", "tsne · celltype", "tsne · Gene3"}


def test_panel_cap_raises(adata):
    with pytest.raises(ValueError, match="exceeds"):
        rs.scatterplot(adata, basis=["umap", "tsne"],
                       color_by=[f"Gene{i}" for i in range(9)],
                       interactive=True, show=False)


def test_alpha_is_global_opacity(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype", alpha=0.3, show=False)
    assert w._spec["options"]["opacity"] == 0.3
    with pytest.raises(TypeError):
        rs.scatterplot(adata, basis="umap", alpha=[0.3, 0.5], show=False)


def test_scanpy_default_palette_and_colors(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype", show=False)
    cols = w._spec["legend"]["colors"]
    assert cols[0] == "#1f77b4"                       # scanpy/vega first colour
    assert all(c.startswith("#") and len(c) == 7 for c in cols)
    # w.colors -> {category: hex}, save-able as adata.uns['x_colors'] = list(...)
    assert w.colors == dict(zip(w._spec["legend"]["names"], cols))


def test_na_color_is_hex_not_css_name(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype", na_color="lightgray", show=False)
    # build a partially-NaN column to exercise the NA path
    obs = adata.obs.copy()
    obs["part"] = pd.Categorical(["T"] * 10 + [None] * (adata.n_obs - 10))
    adata.obs["part"] = obs["part"]
    w2 = rs.scatterplot(adata, basis="umap", color_by="part", show=False)
    lg = w2._spec["legend"]
    assert "NA" in lg["names"]
    na_hex = lg["colors"][lg["names"].index("NA")]
    assert na_hex == "#d3d3d3"                        # hex, not the raw 'lightgray'


def test_compose_grids_have_unique_ids(adata):
    g1 = rs.scatterplot(adata, basis="umap", color_by=["celltype", "Gene3"],
                        interactive=True, show=False)
    g2 = rs.scatterplot(adata, basis="umap", color_by=["celltype", "Gene4"],
                        interactive=True, show=False)
    ids1 = {c._spec["plotId"] for c in g1.children}
    ids2 = {c._spec["plotId"] for c in g2.children}
    assert ids1.isdisjoint(ids2)                      # separate grids never collide


def test_morph_to_sends_message(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype", interactive=True, show=False)
    cap = {}
    w.send = lambda content, buffers=None: cap.update(c=content, b=buffers)
    w.morph_to("tsne", duration=800)
    assert cap["c"]["type"] == "morph" and cap["c"]["duration"] == 800
    assert cap["c"]["xlab"] == "tsne 1" and cap["c"]["n"] == adata.n_obs
    assert len(cap["b"]) == 2                          # x + y buffers
    # static plot has no kernel link
    s = rs.scatterplot(adata, basis="umap", interactive=False, show=False)
    with pytest.raises(AttributeError):
        s.morph_to("tsne")


def test_add_outline_two_band(adata):
    w = rs.scatterplot(adata, basis="umap", color_by="celltype", add_outline=True,
                       outline_color=("black", "white"), outline_width=(0.3, 0.05), show=False)
    o = w._spec["spOutline"]
    assert o["width"] == 0.3 and o["gap"] == 0.05
    assert o["color"] == [0.0, 0.0, 0.0] and o["gapColor"] == [1.0, 1.0, 1.0]


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


# --- z-order / draw depth -------------------------------------------------- #
def test_sort_order_reorders_and_selection_roundtrips(adata):
    w = rs.scatterplot(adata, basis="umap", color="Gene0", interactive=True,
                       sort_order=True, show=False)
    assert w._draw_order is not None             # continuous + sort_order -> reordered
    w.selection = [3, 7, 11]                       # set in DATA indices
    assert w.selection == [3, 7, 11]               # reads back as DATA indices
    assert list(w._selection) != [3, 7, 11]        # stored internally as draw positions


def test_random_state_is_reproducible(adata):
    a = rs.scatterplot(adata, basis="umap", color="celltype", random_state=0, show=False)
    b = rs.scatterplot(adata, basis="umap", color="celltype", random_state=0, show=False)
    assert np.array_equal(a._draw_order, b._draw_order)


def test_sort_order_off_keeps_natural_order(adata):
    w = rs.scatterplot(adata, basis="umap", color="Gene0", sort_order=False, show=False)
    assert w._draw_order is None                   # no reorder -> identity


# --- max_points (downsample huge data) ------------------------------------- #
def test_max_points_subsamples_and_selection_maps(adata):
    w = rs.scatterplot(adata, basis="umap", color="celltype", max_points=20,
                       interactive=True, show=False)
    assert w._spec["n_points"] == 20                 # rendered a subsample
    assert len(w._draw_order) == 20
    target = [int(w._draw_order[5])]                  # an original-data index that's rendered
    w.selection = target
    assert w.selection == target                      # round-trips in ORIGINAL coords


def test_progressive_renders_subset_first(adata):
    # force a subset (adata has 80 cells) -> instant first paint
    w = rs.scatterplot(adata, basis="umap", color="celltype",
                       progressive=True, max_points=20, show=False)
    assert type(w).__name__ == "ReglScatter"          # progressive is live
    assert w._spec["n_points"] == 20                  # subset rendered first
    target = [int(w._draw_order[2])]                  # rendered original index
    w.selection = target
    assert w.selection == target                      # lasso still maps to original


def test_density_subsample_keeps_rare_cells():
    pytest.importorskip("anywidget")
    rng = np.random.default_rng(0)
    big = rng.normal(0, 1, (50000, 2))
    rare = rng.normal(25, 0.3, (120, 2))             # tiny far cluster (rare cell type)
    xy = np.vstack([big, rare])
    df = pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1],
                       "lab": np.r_[np.zeros(50000, int), np.ones(120, int)].astype(str)})
    d = rs.scatterplot(df, x="x", y="y", color="lab", max_points=1000,
                       subsample="density", show=False)
    r = rs.scatterplot(df, x="x", y="y", color="lab", max_points=1000,
                       subsample="random", show=False)
    rare_density = int((np.asarray(d._draw_order) >= 50000).sum())
    rare_random = int((np.asarray(r._draw_order) >= 50000).sum())
    assert rare_density > rare_random * 3             # density preserves rare cells far better


def test_max_points_noop_when_under_threshold(adata):
    w = rs.scatterplot(adata, basis="umap", color="celltype",
                       max_points=10_000, show=False)
    assert w._spec["n_points"] == adata.n_obs         # no subsample


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


def test_compose_accepts_plain_plots(adata):
    from reglscatterpy import _widget
    a = rs.scatterplot(adata, basis="umap", color="celltype", interactive=False)   # static
    b = rs.scatterplot(adata, basis="umap", color="Gene3", interactive=False)
    g = rs.compose([a, b])                                        # static -> HTML grid
    assert type(g).__name__ == "_HtmlGrid"
    assert len(g.panels) == 2 and "text/html" in g._repr_mimebundle_()
    # forcing live links upgrades to a GridBox of live widgets
    g2 = rs.compose([a, b], sync=True)
    assert "GridBox" in [c.__name__ for c in type(g2).__mro__]
    assert all(_widget.is_live_widget(c) for c in g2.children)
    assert all(c._width == 0 for c in g2.children)                # responsive in the grid


def test_single_element_color_list_is_single_plot(adata):
    w = rs.scatterplot(adata, basis="umap", color=["Gene0"], show=False)
    assert type(w).__name__ != "GridBox"   # not a 1-up full-width grid
    assert w._width == 700                  # respects the default fixed width


def test_ncols_in_grid(adata):
    g = rs.scatterplot(adata, basis="umap", color=["celltype", "Gene3", "Gene4"],
                       ncols=3, interactive=False, show=False)
    assert g._cols == 3   # the static HTML grid lays out in 3 columns


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


# --- diff_expression_by: split a selection by an obs column ----------------- #
def _de_anndata(n=180, g=30):
    ad = pytest.importorskip("anndata")
    rng = np.random.RandomState(3)
    counts = rng.poisson(2, size=(n, g)).astype(float)
    lognorm = np.log1p(counts)
    cond = rng.choice(["D30", "Y1", "Y2"], n)
    # plant a real signal in Gene0 for D30 (so it ranks top, logFC finite)
    lognorm[cond == "D30", 0] += 5
    Xs = (lognorm - lognorm.mean(0)) / (lognorm.std(0) + 1e-9)   # z-scored .X
    names = [f"Gene{i}" for i in range(g)]
    obs = pd.DataFrame(
        {"condition": pd.Categorical(cond)},
        index=[f"c{i}" for i in range(n)],
    )
    A = ad.AnnData(X=Xs, obs=obs, var=pd.DataFrame(index=names))
    A.raw = ad.AnnData(X=lognorm, var=pd.DataFrame(index=names))   # log-norm for finite logFC
    A.obsm["X_umap"] = rng.randn(n, 2)
    return A


def test_de_by_returns_scanpy_native_and_one_clean_uns():
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    res = w.diff_expression_by("condition", n=5)
    # native scanpy result: a dict of params + rec.arrays (NOT a dict of frames),
    # and the SAME object that lands in adata.uns
    assert isinstance(res, dict) and "params" in res and "names" in res
    assert res is A.uns["rank_genes_groups"]
    # REAL column + level names (no _rs_by leak); one column per level
    assert res["params"]["groupby"] == "condition"
    assert set(res["names"].dtype.names) == {"D30", "Y1", "Y2"}
    # raw-routing gives finite logFC even though .X is z-scored, params stay clean
    assert res["params"]["layer"] is None
    assert np.isfinite(np.asarray(res["logfoldchanges"]["D30"], float)).all()
    assert res["names"]["D30"][0] == "Gene0"


def test_de_by_pair_returns_native_single_group_and_real_reference():
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    res = w.diff_expression_by("condition", group_a="D30", group_b="Y1",
                               n=5, key_added="d30_v_y1")
    # coherent: same native shape, just one group column
    assert isinstance(res, dict) and set(res["names"].dtype.names) == {"D30"}
    assert res["params"]["reference"] == "Y1"
    assert A.uns["d30_v_y1"]["params"]["reference"] == "Y1"


def test_diff_expression_returns_native():
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    res = w.diff_expression(group_a=list(range(40)), n=5)
    assert isinstance(res, dict) and "params" in res and "names" in res
    assert res is A.uns["rank_genes_groups"]
    assert "A" in res["names"].dtype.names


def test_de_by_bad_column_lists_available():
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    with pytest.raises(KeyError, match="not an obs column"):
        w.diff_expression_by("nope")


def test_de_by_unknown_level_is_rejected():
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    with pytest.raises(ValueError, match="group_a"):
        w.diff_expression_by("condition", group_a="NOPE")


def test_de_by_exposed_on_compose_grid():
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w1 = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    w2 = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    grid = rs.compose([w1, w2], sync=True)
    res = grid.diff_expression_by("condition", n=5)
    assert isinstance(res, dict) and set(res["names"].dtype.names) == {"D30", "Y1", "Y2"}


def test_compose_lasso_on_any_panel_reads_back():
    # Regression: in a linked color_by-list grid, a lasso on a NON-first panel
    # must reach grid.selection / composition. Each panel is its own anywidget
    # model, so the JS cross-panel sync only mirrors the *look* (preventEvent) to
    # the other canvases; the group now reads the last-lassoed panel's model.
    ad = pytest.importorskip("anndata")
    rng = np.random.RandomState(0)
    n, g = 120, 6
    A = ad.AnnData(
        X=rng.rand(n, g),
        obs=pd.DataFrame({"Sample": pd.Categorical(rng.choice(["s1", "s2"], n))},
                         index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"G{i}" for i in range(g)]),
    )
    A.obsm["X_umap"] = rng.randn(n, 2)
    grid = rs.scatterplot(A, basis="X_umap",
                          color_by=["Sample", "G0", "G1", "G2", "G3", "G4"],
                          ncols=3, show=False, interactive=True)
    panels = grid.panels
    assert len(panels) == 6 and grid._synced

    def _orig(panel, pos):
        perm = getattr(panel, "_draw_order", None)
        return sorted(int(perm[p]) for p in pos) if perm is not None else sorted(pos)

    # simulate a lasso on the 4th panel (== JS `model.set("_selection", positions)`)
    panels[3]._selection = [2, 5, 9, 14]
    assert grid._active_sel_panel is panels[3]
    assert sorted(grid.selection) == _orig(panels[3], [2, 5, 9, 14])
    # composition of that lasso comes from the right panel too
    assert int(grid.composition("Sample")["count"].sum()) == 4

    # a later lasso on another panel takes over
    panels[1]._selection = [0, 1]
    assert grid._active_sel_panel is panels[1]
    assert sorted(grid.selection) == _orig(panels[1], [0, 1])


def test_de_engine_scanpy_actually_ran():
    # engine='scanpy' -> real scanpy result: its BH correction is recorded in
    # params, which is how you confirm the Welch fallback did NOT run.
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    res = w.diff_expression_by("condition", n=5, engine="scanpy")
    assert res["params"]["corr_method"] != "none"       # scanpy's BH, not the fallback
    assert res is A.uns["rank_genes_groups"]


def test_de_engine_ttest_forces_welch_fallback():
    # engine='ttest' skips scanpy even when it's installed; the fallback stamps a
    # recognisable params signature.
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    res = w.diff_expression(group_a=list(range(40)), n=5, engine="ttest")
    assert res["params"]["method"] == "t-test"
    assert res["params"]["corr_method"] == "none"
    assert set(res["names"].dtype.names) == {"A"}        # still the native shape


def test_de_engine_invalid_is_rejected():
    pytest.importorskip("scanpy")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    with pytest.raises(ValueError, match="engine="):
        w.diff_expression(group_a=list(range(40)), engine="bogus")


def test_de_engine_gpu_matches_cpu():
    # GPU Welch t-test (cupy) must agree with the CPU Welch t-test to fp, and
    # return the same scanpy-native shape. Skipped when cupy isn't installed.
    pytest.importorskip("scanpy")
    if importlib.util.find_spec("cupy") is None:
        pytest.skip("cupy not installed - GPU path unavailable")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    ga = list(np.where(A.obs["condition"].to_numpy() == "D30")[0])
    cpu = w.diff_expression(group_a=ga, n=8, engine="ttest")
    gpu = w.diff_expression(group_a=ga, n=8, engine="gpu")
    assert gpu["params"]["method"] == "t-test"          # not scanpy, a t-test
    assert set(gpu["names"].dtype.names) == {"A"}
    assert gpu["names"]["A"][0] == cpu["names"]["A"][0]  # same top gene
    np.testing.assert_allclose(
        np.asarray(gpu["scores"]["A"], float),
        np.asarray(cpu["scores"]["A"], float), rtol=1e-5, atol=1e-6,
    )


def test_de_engine_gpu_without_cupy_raises_clearly():
    pytest.importorskip("scanpy")
    if importlib.util.find_spec("cupy") is not None:
        pytest.skip("cupy is installed; the GPU path would run")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    with pytest.raises(ImportError, match="cupy"):
        w.diff_expression(group_a=list(range(40)), engine="gpu")


def test_de_engine_rapids_without_gpu_raises_clearly():
    # No CUDA/rapids-singlecell here -> a clear ImportError naming the fix,
    # never a silent CPU fallback.
    pytest.importorskip("scanpy")
    if importlib.util.find_spec("rapids_singlecell") is not None:
        pytest.skip("rapids-singlecell is installed; the GPU path would run")
    A = _de_anndata()
    w = rs.scatterplot(A, x="X_umap", color_by="condition", show=False, interactive=True)
    with pytest.raises(ImportError, match="rapids-singlecell"):
        w.diff_expression(group_a=list(range(40)), engine="rapids")
