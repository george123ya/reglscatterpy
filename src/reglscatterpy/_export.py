"""Export a plot to a self-contained, offline HTML file.

Jupyter widgets need a live kernel to render, so a *reopened* notebook can't show
the plot until you re-run the cell, and ``nbconvert`` doesn't bake the widget in
either. This mirrors R's ``htmlwidgets::saveWidget``: it writes a standalone
``.html`` that **inlines the widget bundle and the plot's data**, so the file
opens in any browser with no kernel and no internet connection.

Usage::

    w = rs.scatterplot(adata, x="X_umap", color_by="leiden")
    rs.save_html(w, "umap.html")     # or: w.to_html("umap.html")

The whole point is persistence/sharing — the saved file is fully interactive
(pan/zoom, legend, lasso, tooltips, export) but it is a *snapshot*: it has no
kernel, so ``w.selection`` / ``w.annotate`` round-trips are not available in it.
"""

from __future__ import annotations

import base64
import html as _html
import json
import pathlib

__all__ = ["save_html"]

_STATIC = pathlib.Path(__file__).parent / "static" / "widget.js"

# A tiny harness: a shim "model" backed by the snapshot state, then the anywidget
# ESM bundle (inlined as a data: URL module) is asked to render into #root - the
# exact same render({model, el}) entrypoint the live kernel uses, so the offline
# plot is identical to the notebook one.
_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  html, body { margin: 0; padding: 0; background: __PAGEBG__; }
  #root { width: 100%; height: __HEIGHT__px; }
</style>
</head>
<body>
<div id="root"></div>
<script type="module">
const STATE = JSON.parse(atob("__STATE_B64__"));
const listeners = {};
// Minimal stand-in for the ipywidgets model the widget expects. No kernel, so
// save_changes() is a no-op and there is no Python back-channel.
const model = {
  get: (k) => STATE[k],
  set: (k, v) => { STATE[k] = v; (listeners["change:" + k] || []).forEach((f) => f()); },
  save_changes: () => {},
  on: (ev, cb) => { (listeners[ev] = listeners[ev] || []).push(cb); },
  off: (ev, cb) => { if (listeners[ev]) listeners[ev] = listeners[ev].filter((f) => f !== cb); },
};
import("data:text/javascript;base64,__BUNDLE_B64__").then((mod) => {
  const widget = mod.default || mod;
  widget.render({ model, el: document.getElementById("root") });
}).catch((e) => {
  document.getElementById("root").textContent = "reglscatterpy: failed to load (" + e + ")";
});
</script>
</body>
</html>
"""


def save_html(widget, path, title="reglscatterpy plot"):
    """Write a plot to a standalone, offline HTML file (like R's ``saveWidget``).

    Parameters
    ----------
    widget
        A widget returned by :func:`reglscatterpy.scatterplot`.
    path
        Destination ``.html`` path.
    title
        ``<title>`` of the page (defaults to ``"reglscatterpy plot"``).

    Returns
    -------
    str
        The path written, as a string.
    """
    spec = dict(getattr(widget, "_spec", {}) or {})
    if not spec:
        raise ValueError(
            "This widget has no plot spec to export "
            "(was it created by reglscatterpy.scatterplot?)."
        )
    height = int(getattr(widget, "_height", 500) or 500)
    state = {
        "_spec": spec,
        "_selection": [int(i) for i in (getattr(widget, "_selection", []) or [])],
        "_width": int(getattr(widget, "_width", 0) or 0),
        "_height": height,
    }

    bundle_b64 = base64.b64encode(_STATIC.read_bytes()).decode("ascii")
    state_b64 = base64.b64encode(
        json.dumps(state, separators=(",", ":")).encode("utf-8")
    ).decode("ascii")
    page_bg = spec.get("backgroundColor") or "#ffffff"

    page = (
        _TEMPLATE.replace("__TITLE__", _html.escape(str(title)))
        .replace("__PAGEBG__", str(page_bg))
        .replace("__HEIGHT__", str(height))
        .replace("__STATE_B64__", state_b64)
        .replace("__BUNDLE_B64__", bundle_b64)
    )

    out = pathlib.Path(path)
    out.write_text(page, encoding="utf-8")
    return str(out)
