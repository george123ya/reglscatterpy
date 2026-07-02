# Maintainer scripts

Tools for regenerating the README figures and the live docs demo from real
pbmc3k data (downloaded via `sc.datasets.pbmc3k_processed()`), not synthetic toy
data. None of these are runtime dependencies; they need the dev extras plus a
`chromium` on `PATH` and (for stills) Pillow.

## `build_demo.py` — the live GitHub Pages demo

```bash
python scripts/build_demo.py        # -> docs/demo_plot.html  (~0.6 MB, self-contained)
```

`save_html()` bakes the plot **and** the WebGL engine into one offline page. The
docs workflow (`.github/workflows/docs.yml`) copies `docs/**` into the site, so a
push to `main` publishes it at
`https://george123ya.github.io/reglscatterpy/demo_plot.html` (embedded in
`docs/demo.md`, linked from the nav as *Live demo*).

## `build_screens.py` — real-data still PNGs

```bash
python scripts/build_screens.py     # -> assets/umap-categorical.png, assets/umap-continuous.png
```

Renders each panel headlessly (Chromium + SwiftShader WebGL) and trims the
whitespace. Edit the `PANELS` list to add or recolour panels.

## `build_gifs.py` — animated GIFs, fully headless

```bash
python scripts/build_gifs.py        # -> assets/demo.gif
```

A lasso can't be *drawn* headlessly, but the plot's public API can be driven: this
writes a `save_html` plot, then animates `zoomToPoints` + `select` on the live
regl-scatterplot instance (via the puppeteer-core in the sibling
`reglScatterplotR/js` checkout) while screenshotting each frame, and stitches them
with ffmpeg. The hero GIF zooms into a real pbmc3k cluster and highlights it (the
selection ring), so it shows the lasso-to-selection story from real data with no
human in the loop. Needs `chromium`, `node` + that sibling checkout, and ffmpeg.
`_capture_frames.mjs` is the browser-side capturer it calls.

**Can't be automated** (need a live browser recording — see below): the
`compose` **linked grid** (nested sandboxed iframes don't init WebGL headless) and
**`morph_to`** (the target embedding isn't baked into a single-panel export).

## `record_to_gif.sh` — hand-recorded GIFs (morph, linked grid, filters)

For the shots `build_gifs.py` can't do, record by hand from a real browser and
optimise here (`tmp_record_demos.ipynb` stages each shot).

1. Open the live demo (`docs/demo_plot.html`, or the published page) — or run the
   `notebooks/reglscatterpy_tour.ipynb` differential-expression cells for a
   lasso → `diff_expression_by` story.
2. Record the screen region with whatever you use (OBS, the omarchy
   screen-record hotkey, …) to an `.mp4`. A good ~10 s loop: **pan → zoom into a
   cluster → lasso a population → toggle a couple of legend entries.**
3. Convert + optimise:

   ```bash
   scripts/record_to_gif.sh ~/Videos/capture.mp4 assets/demo.gif 760 15
   #                         input                output          width fps
   # optional trim:  scripts/record_to_gif.sh in.mp4 assets/demo.gif 760 15 3 9   # 3s..12s
   ```

`assets/demo.gif` is what the top of the main `README.md` shows.
