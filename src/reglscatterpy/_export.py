"""Export plots to self-contained, offline HTML — one plot or a whole notebook.

Jupyter widgets need a live kernel to render, so a *reopened* notebook can't show
the plot until you re-run the cell, and plain ``jupyter nbconvert --to html``
produces blank plots (the widget bundle isn't captured in the saved state). The
functions here work around that the way R's htmlwidgets do — by **inlining the
widget bundle and the plot's data** into the HTML — with NO R involved; it's pure
Python (base64 + templating).

* :func:`save_html` — one plot to a standalone ``.html`` (like ``saveWidget``).
* :func:`save_notebook_html` — a whole executed notebook to an HTML **report**,
  with every reglscatterpy plot baked in as an interactive, kernel-free figure.

A saved plot is a *snapshot*: it has no kernel, so the Python round-trips
(``w.selection`` / ``w.annotate`` …) are only live in the notebook itself.
"""

from __future__ import annotations

import base64
import functools
import gzip
import html as _html
import json
import pathlib
import uuid

__all__ = ["save_html", "save_notebook_html"]

_STATIC = pathlib.Path(__file__).parent / "static" / "widget.js"

# When True (only ever flipped inside an export kernel), ReglScatter emits a
# self-contained text/html snapshot instead of the live widget view, so
# nbconvert bakes interactive plots into a static report. The live notebook is
# never affected.
_REPORT_REPR = False


def _enable_report_repr():
    """Make widgets render as static HTML snapshots (export kernels only)."""
    global _REPORT_REPR
    _REPORT_REPR = True


def _report_repr_enabled():
    return _REPORT_REPR


@functools.lru_cache(maxsize=1)
def _bundle_gz_b64():
    # gzip then base64: the ~1.3 MB ESM bundle drops to ~0.5 MB inlined, vs
    # ~1.75 MB for plain base64. Decompressed in-browser via DecompressionStream
    # (baseline in all modern browsers since 2023).
    return base64.b64encode(gzip.compress(_STATIC.read_bytes(), 9)).decode("ascii")


# Helper injected once per file: turn the gzip+base64 bundle into a blob URL
# (a Promise) so import() can load it. Shared across all plots in a report.
_LOADER_JS = (
    "window.__rsLoad = window.__rsLoad || function(b64){"
    "var bin = Uint8Array.from(atob(b64), function(c){return c.charCodeAt(0);});"
    "var s = new Blob([bin]).stream().pipeThrough(new DecompressionStream('gzip'));"
    "return new Response(s).arrayBuffer().then(function(buf){"
    "return URL.createObjectURL(new Blob([buf], {type: 'text/javascript'}));});};"
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


def _fragment(widget, bundle_js, div_id=None):
    """An interactive plot as a standalone HTML fragment (div + module script).

    ``bundle_js`` is the JS expression used as the dynamic-import argument: an
    inline ``data:`` URL for a single-file export, or ``window.__rsBundleURL``
    for a report (where the ~1.3 MB bundle is embedded once and shared).
    """
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
        f"  Promise.resolve({bundle_js}).then((u) => import(u)).then((mod) => {{\n"
        f'    (mod.default || mod).render({{ model, el: document.getElementById("{div_id}") }});\n'
        "  });\n"
        "})();\n"
        "</script>"
    )


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
    """Write one plot to a standalone, offline HTML file (like ``saveWidget``).

    Parameters
    ----------
    widget
        A widget returned by :func:`reglscatterpy.scatterplot`.
    path
        Destination ``.html`` path.
    title
        ``<title>`` of the page.

    Returns
    -------
    str
        The path written.
    """
    spec = dict(getattr(widget, "_spec", {}) or {})
    # Decompress the shared bundle once; the fragment imports the resulting URL.
    bundle_js = 'window.__rsLoad("' + _bundle_gz_b64() + '")'
    page = (
        _PAGE.replace("__TITLE__", _html.escape(str(title)))
        .replace("__PAGEBG__", str(spec.get("backgroundColor") or "#ffffff"))
        .replace("__LOADER__", _LOADER_JS)
        .replace("__FRAGMENT__", _fragment(widget, bundle_js))
    )
    out = pathlib.Path(path)
    out.write_text(page, encoding="utf-8")
    return str(out)


def report_fragment(widget):
    """The report-mode fragment for ``widget`` (shared-bundle form).

    Used by ReglScatter._repr_mimebundle_ when report rendering is enabled;
    pairs with the ``window.__rsBundleURL`` definition injected by
    :func:`save_notebook_html`.
    """
    return _fragment(widget, "window.__rsBundleURL")


def save_notebook_html(
    notebook,
    out_path=None,
    execute=True,
    kernel_name="python3",
    timeout=600,
):
    """Convert a whole notebook to a self-contained HTML **report**.

    Unlike plain ``jupyter nbconvert --to html`` — which leaves reglscatterpy
    plots blank because the widget bundle isn't saved in the notebook — this
    re-executes the notebook with plots rendered as inlined, interactive,
    kernel-free figures, so the report is fully self-contained and offline. The
    ~1.3 MB widget bundle is embedded **once** and shared across all plots.

    Parameters
    ----------
    notebook
        Path to the ``.ipynb``.
    out_path
        Destination ``.html`` (defaults to the notebook name with ``.html``).
    execute
        Re-run the notebook before exporting (default ``True``). This is what
        makes the plots appear; with ``False`` the plots may be blank unless the
        notebook was already executed in report mode.
    kernel_name
        Jupyter kernel to execute with (default ``"python3"``).
    timeout
        Per-cell execution timeout in seconds.

    Returns
    -------
    str
        The path written.
    """
    try:
        import nbformat
        from nbconvert import HTMLExporter
    except ModuleNotFoundError as exc:  # pragma: no cover - optional dep
        raise ModuleNotFoundError(
            "save_notebook_html() needs nbconvert + nbformat (and, to execute, "
            "nbclient + ipykernel). Install with: pip install nbconvert ipykernel"
        ) from exc

    nb_path = pathlib.Path(notebook)
    nb = nbformat.read(str(nb_path), as_version=4)

    if execute:
        from nbconvert.preprocessors import ExecutePreprocessor

        # Inject a hidden setup cell that flips report mode ON inside the export
        # kernel only, so widgets emit static HTML snapshots during this run.
        setup = nbformat.v4.new_code_cell(
            "import reglscatterpy as _rs; _rs._export._enable_report_repr()"
        )
        setup.metadata["tags"] = ["rs-report-setup"]
        nb.cells.insert(0, setup)
        ExecutePreprocessor(timeout=timeout, kernel_name=kernel_name).preprocess(
            nb, {"metadata": {"path": str(nb_path.parent)}}
        )
        # Drop the injected cell so it doesn't show in the report.
        nb.cells = [
            c
            for c in nb.cells
            if "rs-report-setup" not in c.get("metadata", {}).get("tags", [])
        ]

    body, _ = HTMLExporter().from_notebook_node(nb)
    # Embed the shared (gzipped) bundle once; every plot fragment imports the one
    # decompressed blob URL via the shared __rsBundleURL promise.
    inject = (
        "<script>"
        + _LOADER_JS
        + 'window.__rsBundleURL = window.__rsLoad("'
        + _bundle_gz_b64()
        + '");</script>'
    )
    if "</head>" in body:
        body = body.replace("</head>", inject + "\n</head>", 1)
    else:  # pragma: no cover - HTMLExporter always emits a head
        body = inject + body

    out = pathlib.Path(out_path) if out_path else nb_path.with_suffix(".html")
    out.write_text(body, encoding="utf-8")
    return str(out)
