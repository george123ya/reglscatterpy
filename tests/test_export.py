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


def test_save_notebook_html(tmp_path):
    pytest.importorskip("nbconvert")
    pytest.importorskip("nbformat")
    pytest.importorskip("ipykernel")
    import nbformat as nbf

    nb = nbf.v4.new_notebook()
    nb.cells = [
        nbf.v4.new_markdown_cell("# Report"),
        nbf.v4.new_code_cell(
            "import numpy as np, pandas as pd, reglscatterpy as rs\n"
            "df = pd.DataFrame({'x': np.arange(30), 'y': np.arange(30),\n"
            "                   'ct': ['a','b','c']*10})\n"
            "rs.scatterplot(df, x='x', y='y', color_by='ct')"
        ),
    ]
    nb_path = tmp_path / "r.ipynb"
    nbf.write(nb, str(nb_path))

    out = rs.save_notebook_html(nb_path, tmp_path / "r.html")
    html = pathlib.Path(out).read_text(encoding="utf-8")
    # The shared bundle is embedded exactly once, and the plot is baked in as a
    # static fragment that imports it (not a live widget view).
    assert html.count('window.__rsBundleURL = window.__rsLoad(') == 1
    assert "window.__rsBundleURL).then" in html
    assert 'id="rs_' in html  # the inlined plot div fragment
