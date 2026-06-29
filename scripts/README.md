# Maintainer scripts

Tools for regenerating the README figures and the live docs demo from **real
pbmc3k** data (`notebooks/data/pbmc3k_processed.h5ad`), not synthetic toy data.
None of these are runtime dependencies — they need the dev extras plus a
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

## `record_to_gif.sh` — the README hero GIF

Headless capture can't drive a **lasso** gesture, so the animated GIF is recorded
by hand from a real browser and then optimised here.

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
