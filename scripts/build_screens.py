#!/usr/bin/env python3
"""Render REAL-data still images for the README from pbmc3k (no synthetic data).

    python scripts/build_screens.py

For each panel it writes a self-contained plot with ``save_html`` and screenshots
it with headless Chromium (SwiftShader WebGL), then trims the whitespace. Outputs
land in ``assets/``. Needs a ``chromium`` (or ``chromium-browser``) on PATH and
Pillow; both are dev-only, so this is a maintainer tool, not a runtime dependency.
"""
import pathlib
import shutil
import subprocess
import tempfile

import scanpy as sc
from PIL import Image

import reglscatterpy as rs

ROOT = pathlib.Path(__file__).resolve().parent.parent
ASSETS = ROOT / "assets"

CHROME = shutil.which("chromium") or shutil.which("chromium-browser") or shutil.which("google-chrome-stable")

# (output filename, scatterplot kwargs)
PANELS = [
    ("umap-categorical.png", dict(basis="umap", color_by="louvain")),
    ("umap-continuous.png", dict(basis="umap", color_by="NKG7")),  # NK / CD8 cytotoxic marker
]


def shoot(html: pathlib.Path, png: pathlib.Path, w=1180, h=620):
    if CHROME is None:
        raise SystemExit("no Chromium on PATH — install chromium to render stills.")
    subprocess.run(
        [CHROME, "--headless=new", "--no-sandbox", "--hide-scrollbars",
         "--enable-unsafe-swiftshader", "--use-gl=angle", "--use-angle=swiftshader",
         "--ignore-gpu-blocklist", "--run-all-compositor-stages-before-draw",
         f"--window-size={w},{h}", "--virtual-time-budget=25000",
         f"--screenshot={png}", html.as_uri()],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )


def autocrop(png: pathlib.Path, pad=12, bg=(255, 255, 255)):
    im = Image.open(png).convert("RGB")
    from PIL import ImageChops
    diff = ImageChops.difference(im, Image.new("RGB", im.size, bg))
    box = diff.getbbox()
    if box:
        x0, y0, x1, y1 = box
        x0, y0 = max(0, x0 - pad), max(0, y0 - pad)
        x1, y1 = min(im.width, x1 + pad), min(im.height, y1 + pad)
        im = im.crop((x0, y0, x1, y1))
    im.save(png)
    return im.size


def main():
    ASSETS.mkdir(exist_ok=True)
    adata = sc.datasets.pbmc3k_processed()   # downloads once
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        for name, kw in PANELS:
            w = rs.scatterplot(adata, show=False, **kw)
            html = td / (name + ".html")
            rs.save_html(w, html, title=name)
            png = ASSETS / name
            shoot(html, png)
            size = autocrop(png)
            print(f"wrote {png.relative_to(ROOT)}  {size[0]}x{size[1]}")


if __name__ == "__main__":
    main()
