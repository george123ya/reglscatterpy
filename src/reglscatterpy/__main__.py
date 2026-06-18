"""Command-line report export:  ``reglscatterpy-report notebook.ipynb``.

Convert a notebook to a self-contained, offline HTML report with every plot
baked in. By default it does NOT re-run the notebook (it uses the existing
outputs), so a heavy notebook isn't re-executed — run it once with
``reglscatterpy.record_html()`` at the top, then export. Use ``--execute`` to
re-run a notebook that wasn't recorded.
"""

from __future__ import annotations

import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="reglscatterpy-report",
        description="Convert a notebook to a self-contained, offline HTML report "
        "with interactive reglscatterpy plots baked in.",
    )
    parser.add_argument("notebook", help="path to the .ipynb to export")
    parser.add_argument(
        "-o", "--output", default=None, help="output .html (default: notebook name)"
    )
    parser.add_argument(
        "--execute",
        action="store_true",
        help="re-run the notebook before exporting (default: use existing outputs)",
    )
    parser.add_argument(
        "--kernel", default="python3", help="kernel to use with --execute"
    )
    parser.add_argument(
        "--timeout", type=int, default=600, help="per-cell timeout with --execute"
    )
    args = parser.parse_args(argv)

    from ._export import save_notebook_html

    out = save_notebook_html(
        args.notebook,
        args.output,
        execute=args.execute,
        kernel_name=args.kernel,
        timeout=args.timeout,
    )
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
