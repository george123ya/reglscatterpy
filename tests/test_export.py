"""save_html / to_html: standalone offline HTML export."""

from __future__ import annotations

import base64
import json

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
    # The widget bundle and the plot state are both inlined.
    assert 'import("data:text/javascript;base64,' in html
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
