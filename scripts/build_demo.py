#!/usr/bin/env python3
"""Build the live GitHub Pages demo from REAL pbmc3k data (not synthetic).

    python scripts/build_demo.py

Writes a self-contained, interactive WebGL page to ``docs/demo_plot.html`` — the
WebGL bundle is baked in, so it runs offline and on GitHub Pages with no kernel.
The docs workflow (.github/workflows/docs.yml) copies docs/** into the site, so a
push to main publishes it at https://george123ya.github.io/reglscatterpy/demo_plot.html
"""
import pathlib

import scanpy as sc

import reglscatterpy as rs

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "demo_plot.html"


def main():
    adata = sc.datasets.pbmc3k_processed()   # downloads once
    # A single, full-width interactive embedding of real PBMCs: pan / zoom / lasso
    # and toggle cell types in the legend. (A single panel renders reliably even
    # in headless SwiftShader; a linked compose() grid works in a real browser but
    # nests iframes that some headless screenshotters can't paint — so for the
    # public demo we keep one robust panel.)
    w = rs.scatterplot(adata, basis="umap", color_by="louvain", show=False)
    rs.save_html(w, OUT, title="reglscatterpy — live demo (pbmc3k, 2,638 PBMCs, louvain)")
    print(f"wrote {OUT}  ({OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
