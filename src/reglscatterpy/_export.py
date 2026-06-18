"""Export plots to self-contained, offline HTML — one plot or a whole notebook.

Jupyter widgets need a live kernel to render, so a *reopened* notebook can't show
the plot until you re-run the cell, and plain ``jupyter nbconvert --to html``
produces blank plots (the widget bundle isn't saved in the notebook). The
functions here work around that the way R's htmlwidgets do — by **inlining the
widget bundle and the plot's data** into the HTML — with NO R involved; it's pure
Python (gzip + base64 + templating).

Three entry points:

* :func:`save_html` — one plot to a standalone ``.html`` (like ``saveWidget``).
* :func:`record_html` — flip this on at the top of a notebook; every plot then
  bakes a static, interactive copy into its *own cell output* as you run it. The
  plots survive reopening the notebook and ``jupyter nbconvert --to html``
  **with no re-execution**. (Trade-off: a recorded plot is one-way — pan/zoom/
  lasso work, but ``w.selection`` no longer round-trips to Python. Turn it off
  for the live round-trip.)
* :func:`save_notebook_html` — convert a whole notebook to one HTML report.
  Defaults to ``execute=False``: it uses the existing (e.g. recorded) outputs,
  so a heavy notebook is **not re-run**. Pass ``execute=True`` to re-run a
  notebook that wasn't recorded.
"""

from __future__ import annotations

import base64
import functools
import gzip
import html as _html
import json
import os
import pathlib
import uuid
import warnings

__all__ = ["save_html", "save_notebook_html", "record_html"]

_STATIC = pathlib.Path(__file__).parent / "static" / "widget.js"

_WIDGET_VIEW_MIME = "application/vnd.jupyter.widget-view+json"

# Static-render modes. _REPORT_REPR is flipped only inside an export kernel
# (save_notebook_html(execute=True)); _RECORD is flipped by the user via
# record_html() in a live notebook. Either makes ReglScatter emit a
# self-contained text/html snapshot instead of the live widget view. The normal
# interactive widget is unaffected unless one of these is on.
_REPORT_REPR = False
_RECORD = False


def _enable_report_repr():
    """Make widgets render as shared-bundle HTML snapshots (export kernel only)."""
    global _REPORT_REPR
    _REPORT_REPR = True


def _report_repr_enabled():
    return _REPORT_REPR


def record_html(enable=True):
    """Bake a static, offline copy of every plot into its cell output.

    Call once near the top of a notebook::

        import reglscatterpy as rs
        rs.record_html()

    From then on each ``scatterplot(...)`` displays a self-contained interactive
    plot that is stored in the notebook output — so reopening the notebook or
    running ``jupyter nbconvert --to html`` shows the plots with **no
    re-execution**. Recorded plots are one-way (no ``w.selection`` round-trip);
    call ``rs.record_html(False)`` to return to the live, kernel-linked widget.
    """
    global _RECORD
    _RECORD = bool(enable)


def _record_enabled():
    return _RECORD


@functools.lru_cache(maxsize=1)
def _bundle_gz_b64():
    # gzip then base64: the ~1.3 MB ESM bundle drops to ~0.5 MB inlined, vs
    # ~1.75 MB for plain base64. Decompressed in-browser via DecompressionStream
    # (baseline in all modern browsers since 2023).
    return base64.b64encode(gzip.compress(_STATIC.read_bytes(), 9)).decode("ascii")


# Helper: turn the gzip+base64 bundle into a blob URL (a Promise) so import()
# can load it. Idempotent + shared across all plots in one document.
_LOADER_JS = (
    "window.__rsLoad = window.__rsLoad || function(b64){"
    "var bin = Uint8Array.from(atob(b64), function(c){return c.charCodeAt(0);});"
    "var s = new Blob([bin]).stream().pipeThrough(new DecompressionStream('gzip'));"
    "return new Response(s).arrayBuffer().then(function(buf){"
    "return URL.createObjectURL(new Blob([buf], {type: 'text/javascript'}));});};"
)


def _bundle_call():
    """JS that decompresses the bundle once into the shared ``__rsBundleURL``."""
    return (
        'window.__rsBundleURL = window.__rsBundleURL || window.__rsLoad("'
        + _bundle_gz_b64()
        + '");'
    )


def _state(widget):
    spec = dict(getattr(widget, "_spec", {}) or {})
    if not spec:
        raise ValueError(
            "This widget has no plot spec to export "
            "(was it created by reglscatterpy.scatterplot?)."
        )
    height = int(getattr(widget, "_height", 500) or 500)
    return {
        "_spec": spec,
        "_selection": [int(i) for i in (getattr(widget, "_selection", []) or [])],
        "_width": int(getattr(widget, "_width", 0) or 0),
        "_height": height,
    }, height


def _fragment(widget, div_id=None):
    """An interactive plot as an HTML fragment that imports ``__rsBundleURL``."""
    state, height = _state(widget)
    div_id = div_id or ("rs_" + uuid.uuid4().hex[:12])
    state_b64 = base64.b64encode(
        json.dumps(state, separators=(",", ":")).encode("utf-8")
    ).decode("ascii")
    return (
        f'<div id="{div_id}" style="width:100%;height:{height}px"></div>\n'
        '<script type="module">\n'
        "(function(){\n"
        f'  const STATE = JSON.parse(atob("{state_b64}"));\n'
        "  const listeners = {};\n"
        "  const model = {\n"
        "    get: (k) => STATE[k],\n"
        '    set: (k, v) => { STATE[k] = v; (listeners["change:" + k] || []).forEach((f) => f()); },\n'
        "    save_changes: () => {},\n"
        "    on: (ev, cb) => { (listeners[ev] = listeners[ev] || []).push(cb); },\n"
        "    off: (ev, cb) => { if (listeners[ev]) listeners[ev] = listeners[ev].filter((f) => f !== cb); },\n"
        "  };\n"
        f"  Promise.resolve(window.__rsBundleURL).then((u) => import(u)).then((mod) => {{\n"
        f'    (mod.default || mod).render({{ model, el: document.getElementById("{div_id}") }});\n'
        "  });\n"
        "})();\n"
        "</script>"
    )


def report_fragment(widget):
    """Shared-bundle fragment (the bundle is injected once by the exporter)."""
    return _fragment(widget)


def record_fragment(widget):
    """Self-contained fragment: carries the bundle so it works on reopen too.

    Each recorded output inlines the (gzipped) bundle, but a runtime guard means
    it is only decompressed once per page; :func:`save_notebook_html` later
    de-dupes the repeated copies down to one in the report file.
    """
    setup = "<script>" + _LOADER_JS + _bundle_call() + "</script>\n"
    return setup + _fragment(widget)


_PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>html, body { margin: 0; padding: 0; background: __PAGEBG__; }</style>
<script>__LOADER__</script>
</head>
<body>
__FRAGMENT__
</body>
</html>
"""


def save_html(widget, path, title="reglscatterpy plot"):
    """Write one plot to a standalone, offline HTML file (like ``saveWidget``)."""
    spec = dict(getattr(widget, "_spec", {}) or {})
    page = (
        _PAGE.replace("__TITLE__", _html.escape(str(title)))
        .replace("__PAGEBG__", str(spec.get("backgroundColor") or "#ffffff"))
        .replace("__LOADER__", _LOADER_JS + _bundle_call())
        .replace("__FRAGMENT__", _fragment(widget))
    )
    out = pathlib.Path(path)
    out.write_text(page, encoding="utf-8")
    return str(out)


def _strip_widget_views(nb):
    """Drop widget-view outputs so nbconvert renders our text/html instead.

    Returns the number of plot outputs that had a widget view but no recorded
    text/html (those will be blank — the notebook wasn't run with record_html).
    """
    unrecorded = 0
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        for out in cell.get("outputs", []):
            data = out.get("data")
            if not isinstance(data, dict) or _WIDGET_VIEW_MIME not in data:
                continue
            if "text/html" in data:
                del data[_WIDGET_VIEW_MIME]  # keep our static plot
            else:
                unrecorded += 1
    return unrecorded


def save_notebook_html(
    notebook,
    out_path=None,
    execute=False,
    kernel_name="python3",
    timeout=600,
):
    """Convert a whole notebook to a self-contained HTML **report**.

    Unlike plain ``jupyter nbconvert --to html`` (which leaves reglscatterpy
    plots blank), this bakes every plot in as an interactive, kernel-free figure,
    sharing **one** copy of the bundle across all plots.

    Parameters
    ----------
    notebook
        Path to the ``.ipynb``.
    out_path
        Destination ``.html`` (defaults to the notebook name with ``.html``).
    execute
        Re-run the notebook before exporting. **Default ``False``** — it uses the
        notebook's existing outputs, so a heavy notebook is *not* re-run. This
        works when the notebook was run with :func:`record_html` (recommended).
        Pass ``True`` to re-run a notebook that wasn't recorded.
    kernel_name, timeout
        Passed to the execute step when ``execute=True``.

    Returns
    -------
    str
        The path written.
    """
    # Recursion guard: when this runs inside an export kernel (execute=True),
    # any save_notebook_html call in the notebook itself is a safe no-op.
    if os.environ.get("REGLSCATTERPY_EXPORTING") == "1":
        return ""

    try:
        import nbformat
        from nbconvert import HTMLExporter
    except ModuleNotFoundError as exc:  # pragma: no cover - optional dep
        raise ModuleNotFoundError(
            "save_notebook_html() needs nbconvert + nbformat (and, to execute, "
            "nbclient + ipykernel). Install with: pip install 'reglscatterpy[report]'"
        ) from exc

    nb_path = pathlib.Path(notebook)
    nb = nbformat.read(str(nb_path), as_version=4)

    if execute:
        from nbconvert.preprocessors import ExecutePreprocessor

        setup = nbformat.v4.new_code_cell(
            "import reglscatterpy as _rs; _rs._export._enable_report_repr()"
        )
        setup.metadata["tags"] = ["rs-report-setup"]
        nb.cells.insert(0, setup)
        os.environ["REGLSCATTERPY_EXPORTING"] = "1"
        try:
            ExecutePreprocessor(timeout=timeout, kernel_name=kernel_name).preprocess(
                nb, {"metadata": {"path": str(nb_path.parent)}}
            )
        finally:
            os.environ.pop("REGLSCATTERPY_EXPORTING", None)
        nb.cells = [
            c
            for c in nb.cells
            if "rs-report-setup" not in c.get("metadata", {}).get("tags", [])
        ]
    else:
        unrecorded = _strip_widget_views(nb)
        if unrecorded:
            warnings.warn(
                f"{unrecorded} plot(s) have no recorded HTML and will be blank. "
                "Run the notebook with reglscatterpy.record_html() at the top, "
                "or call save_notebook_html(..., execute=True).",
                stacklevel=2,
            )

    body, _ = HTMLExporter().from_notebook_node(nb)
    # De-dupe: collapse every inlined bundle copy to the one shared __rsBundleURL,
    # then embed that single copy once.
    gz_call = _bundle_call()
    body = body.replace("<script>" + _LOADER_JS + gz_call + "</script>\n", "")
    inject = "<script>" + _LOADER_JS + gz_call + "</script>"
    if "</head>" in body:
        body = body.replace("</head>", inject + "\n</head>", 1)
    else:  # pragma: no cover - HTMLExporter always emits a head
        body = inject + body

    out = pathlib.Path(out_path) if out_path else nb_path.with_suffix(".html")
    out.write_text(body, encoding="utf-8")
    return str(out)
