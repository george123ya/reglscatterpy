#!/usr/bin/env python3
"""Generate the README's animated GIFs from REAL pbmc3k — fully headless.

A lasso can't be faked headlessly, but the plot's public API can: this writes a
self-contained ``save_html`` plot, then drives ``zoomToPoints`` / ``select`` on
the live regl-scatterplot instance while screenshotting each frame (Chromium +
SwiftShader WebGL via puppeteer), and assembles the frames into a GIF with
ffmpeg. Maintainer tool — needs a ``chromium`` on PATH, the puppeteer-core in the
sibling ``reglScatterplotR/js`` checkout, ffmpeg, and Pillow-free.

    python scripts/build_gifs.py            # -> assets/demo.gif (+ others)
"""
import json
import pathlib
import subprocess
import tempfile

import numpy as np
import scanpy as sc

import reglscatterpy as rs

ROOT = pathlib.Path(__file__).resolve().parent.parent
ASSETS = ROOT / "assets"
CAPTURE = ROOT / "scripts" / "_capture_frames.mjs"

# Each shot: output gif, the cluster to zoom into + lasso-highlight, and timing.
SHOTS = [
    dict(name="demo.gif", color_by="louvain", cluster="NK cells",
         w=720, h=560, fps=14, width=680,
         plan=dict(holdStart=5, zoomIn=16, holdSel=10, zoomOut=13, holdEnd=4, padding=0.35)),
]


def build_html(adata, color_by, tmp):
    # sort_order=False -> draw order is identity, so obs-row indices == the
    # rendered positions the plot's select()/zoomToPoints() expect.
    w = rs.scatterplot(adata, basis="umap", color_by=color_by,
                       sort_order=False, show=False, width=760)
    html = tmp / f"{color_by}.html"
    rs.save_html(w, html, title=color_by)
    return html


def capture(html, outdir, params):
    pj = outdir / "params.json"
    pj.write_text(json.dumps(params))
    subprocess.run(["node", str(CAPTURE), str(html), str(outdir), str(pj)],
                   check=True, cwd=str(ROOT))


def assemble(outdir, gif, fps, width):
    frames = sorted(outdir.glob("f*.png"))
    if not frames:
        raise SystemExit(f"no frames captured in {outdir}")
    pal = outdir / "pal.png"
    vf = f"scale={width}:-1:flags=lanczos"
    subprocess.run(["ffmpeg", "-hide_banner", "-loglevel", "error", "-y",
                    "-framerate", str(fps), "-i", str(outdir / "f%03d.png"),
                    "-vf", f"{vf},palettegen=stats_mode=diff", str(pal)], check=True)
    subprocess.run(["ffmpeg", "-hide_banner", "-loglevel", "error", "-y",
                    "-framerate", str(fps), "-i", str(outdir / "f%03d.png"),
                    "-i", str(pal),
                    "-lavfi", f"{vf}[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=3",
                    str(gif)], check=True)
    return len(frames), gif.stat().st_size


def main():
    ASSETS.mkdir(exist_ok=True)
    adata = sc.datasets.pbmc3k_processed()   # downloads once
    with tempfile.TemporaryDirectory() as td:
        td = pathlib.Path(td)
        for shot in SHOTS:
            html = build_html(adata, shot["color_by"], td)
            sel = np.where(adata.obs[shot["color_by"]].to_numpy() == shot["cluster"])[0]
            params = dict(w=shot["w"], h=shot["h"], select=sel.tolist(), **shot["plan"])
            outdir = td / shot["name"].replace(".", "_")
            outdir.mkdir()
            capture(html, outdir, params)
            n, size = assemble(outdir, ASSETS / shot["name"], shot["fps"], shot["width"])
            print(f"wrote {(ASSETS / shot['name']).relative_to(ROOT)}  "
                  f"{n} frames, {size // 1024} KB  (lasso={shot['cluster']}, {len(sel)} cells)")


if __name__ == "__main__":
    main()
