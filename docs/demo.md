# Live demo

Real **pbmc3k** (2,638 PBMCs) in an interactive WebGL embedding — **pan, zoom,
lasso**, and toggle cell types in the legend. It runs entirely in your browser
with no kernel: the plot and the WebGL engine are baked into one self-contained
page, exactly what [`save_html()`](api.md#reglscatterpy.save_html) produces.

[Open the demo full screen ↗](demo_plot.html){ .md-button .md-button--primary }

<iframe src="../demo_plot.html" title="reglscatterpy live demo — pbmc3k"
        loading="lazy"
        style="width:100%; height:640px; border:1px solid var(--md-default-fg-color--lightest); border-radius:8px; margin-top:1em;">
</iframe>

!!! tip "Reproduce it"
    The page is generated from the bundled pbmc3k file by
    [`scripts/build_demo.py`](https://github.com/george123ya/reglscatterpy/blob/main/scripts/build_demo.py):

    ```bash
    python scripts/build_demo.py        # -> docs/demo_plot.html
    ```

    A linked **`compose([...])` grid** (lasso in one embedding, highlight in the
    other) works the same way in a live browser — see the *User guide*.
