"""save_html / to_html: standalone offline HTML export."""

from __future__ import annotations

import base64
import json
import pathlib

import numpy as np
import pandas as pd
import pytest

import reglscatterpy as rs


def _widget():
    np.random.seed(0)
    df = pd.DataFrame(
        {
            "x": np.random.rand(120),
            "y": np.random.rand(120),
            "ct": np.random.choice(list("ABC"), 120),
        }
    )
    return rs.scatterplot(df, x="x", y="y", color_by="ct")


def test_save_html_writes_self_contained_file(tmp_path):
    w = _widget()
    out = tmp_path / "plot.html"
    rs.save_html(w, out)
    assert out.exists()
    html = out.read_text(encoding="utf-8")
    # Self-contained: no external script/style URLs (no CDN, no network).
    assert "http://" not in html and "https://" not in html
    # The bundle is inlined gzip+base64 and decompressed in-browser; the plot
    # state is inlined too.
    assert "window.__rsLoad(" in html
    assert "DecompressionStream('gzip')" in html
    assert "JSON.parse(atob(" in html


def test_save_html_embeds_the_plot_data(tmp_path):
    w = _widget()
    out = tmp_path / "plot.html"
    w.to_html(out, title="my plot")  # method form
    html = out.read_text(encoding="utf-8")
    assert "<title>my plot</title>" in html
    # Decode the embedded state and confirm the spec round-trips.
    marker = 'JSON.parse(atob("'
    b64 = html.split(marker, 1)[1].split('"', 1)[0]
    state = json.loads(base64.b64decode(b64))
    assert "_spec" in state and state["_spec"]
    assert state["_spec"]["n_points"] == 120


def test_save_html_returns_path(tmp_path):
    w = _widget()
    out = tmp_path / "p.html"
    assert rs.save_html(w, out) == str(out)


def test_save_html_rejects_empty(tmp_path):
    class Dummy:
        _spec = {}

    with pytest.raises(ValueError):
        rs.save_html(Dummy(), tmp_path / "x.html")


def test_record_mode_emits_self_contained_html():
    rs.record_html()
    try:
        w = _widget()
        bundle = w._repr_mimebundle_()
        assert "text/html" in bundle
        assert "application/vnd.jupyter.widget-view+json" not in bundle
        frag = bundle["text/html"]
        # Self-contained: carries the bundle so a reopened notebook works.
        assert "DecompressionStream('gzip')" in frag
        assert "window.__rsBundleURL" in frag
    finally:
        rs.record_html(False)
    # Off again -> the default is now the static iframe snapshot (no live comm).
    off = _widget()._repr_mimebundle_()
    data = off[0] if isinstance(off, tuple) else off
    assert "application/vnd.jupyter.widget-view+json" not in data
    assert "<iframe" in data["text/html"]
    # The live, kernel-linked widget view is opt-in via interactive=True.
    wi = rs.scatterplot(
        pd.DataFrame({"x": [0.0, 1.0], "y": [1.0, 0.0], "ct": ["a", "b"]}),
        x="x", y="y", color_by="ct", interactive=True,
    )
    di = wi._repr_mimebundle_()
    di = di[0] if isinstance(di, tuple) else di
    assert "application/vnd.jupyter.widget-view+json" in di


def _make_nb(tmp_path, record):
    import nbformat as nbf

    setup = "import numpy as np, pandas as pd, reglscatterpy as rs\n"
    if record:
        setup += "rs.record_html()\n"
    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("# Report"),
        nbf.v4.new_code_cell(setup),
        nbf.v4.new_code_cell(
            "df = pd.DataFrame({'x': np.arange(30), 'y': np.arange(30),\n"
            "                   'ct': ['a','b','c']*10})\n"
            "rs.scatterplot(df, x='x', y='y', color_by='ct')"
        ),
    ]
    p = tmp_path / "r.ipynb"
    nbf.write(nb, str(p))
    return p


def test_save_notebook_html_execute(tmp_path):
    pytest.importorskip("nbconvert")
    pytest.importorskip("ipykernel")
    nb_path = _make_nb(tmp_path, record=False)
    out = rs.save_notebook_html(nb_path, tmp_path / "r.html", execute=True)
    html = pathlib.Path(out).read_text(encoding="utf-8")
    # Bundle embedded exactly once; the plot is baked in as a static fragment.
    assert html.count(rs._export._bundle_gz_b64()[:60]) == 1
    assert "window.__rsBundleURL).then" in html
    assert 'id="rs_' in html


def test_save_notebook_html_no_execute_uses_recorded_outputs(tmp_path):
    pytest.importorskip("nbconvert")
    pytest.importorskip("ipykernel")
    from nbconvert.preprocessors import ExecutePreprocessor
    import nbformat as nbf

    # Run the recorded notebook ONCE (as a user would), then export with NO
    # re-execution and confirm the plot is still baked in (de-duped to one bundle).
    nb_path = _make_nb(tmp_path, record=True)
    nb = nbf.read(str(nb_path), as_version=4)
    ExecutePreprocessor(timeout=120, kernel_name="python3").preprocess(
        nb, {"metadata": {"path": str(tmp_path)}}
    )
    nbf.write(nb, str(nb_path))

    out = rs.save_notebook_html(nb_path, tmp_path / "r.html")  # execute=False
    html = pathlib.Path(out).read_text(encoding="utf-8")
    assert html.count(rs._export._bundle_gz_b64()[:60]) == 1
    assert 'id="rs_' in html
