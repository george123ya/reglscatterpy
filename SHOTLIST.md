# Media to capture for the README

Drop the files in `assets/` with these exact names (the README already links to
them). Use the bundled-style data so anyone can reproduce them, e.g.
`scanpy.datasets.pbmc3k_processed()` (a UMAP with a `louvain` column and genes
like `CST3`, `MS4A1`, `NKG7`).

## Stills (PNG)

| File | What to show | How |
|------|--------------|-----|
| `assets/umap-categorical.png` | UMAP coloured by a cluster column, **frosted legend visible** in a corner | `rs.scatterplot(adata, x="X_umap", color_by="louvain")` |
| `assets/umap-continuous.png` | Same UMAP coloured by a gene, **colour bar visible** | `rs.scatterplot(adata, x="X_umap", color_by="CST3", continuous_palette="viridis", vmax="p99")` |
| `assets/filter-sliders.png` | The `filter_by` panel: a histogram with the dual-handle range brush, some points dimmed | `rs.scatterplot(adata, x="X_umap", color_by="louvain", filter_by=["n_genes"])` then drag a handle in |
| `assets/linked-grid.png` | Two embeddings side by side, one zoomed (to prove the camera is synced) | `compose([scatterplot(adata, x="X_umap", color_by="louvain"), scatterplot(adata, x="X_pca", color_by="louvain")])` |

Capture the canvas region only (not the whole browser chrome). ~1400 px wide is
plenty; PNG, not JPG, so the points stay crisp.

## Hero animation (GIF)

`assets/demo.gif` — one ~8–12 s clip, in this order:

1. Pan and zoom (scroll) around the UMAP.
2. Drag the legend to another corner, then click a category to filter it out;
   shift-click a second to extend.
3. Switch to the lasso tool and circle a cluster.

Keep it short and loopable. Target ≤ ~4 MB so it loads fast on PyPI/GitHub.

### Recording → GIF (Linux/Wayland)

```bash
# record a region to mp4 (pick the plot area)
wf-recorder -g "$(slurp)" -f demo.mp4
# trim if needed:  ffmpeg -ss 2 -t 10 -i demo.mp4 -c copy demo_cut.mp4
# mp4 -> small GIF (high-quality palette, ~12 fps, 760 px wide to match README)
ffmpeg -i demo.mp4 -vf "fps=12,scale=760:-1:flags=lanczos,palettegen" -y pal.png
ffmpeg -i demo.mp4 -i pal.png -lavfi "fps=12,scale=760:-1:flags=lanczos[x];[x][1:v]paletteuse" -y assets/demo.gif
```

> Tip: an `.mp4` is far smaller than a GIF. If you'd rather embed video, drag the
> `.mp4` into the GitHub README via the web editor (GitHub hosts it) — but PyPI
> only renders images, so keep `demo.gif` for the PyPI page.
