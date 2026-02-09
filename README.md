# Interactive Volcano Plot

A polished, interactive **volcano plot** built with **D3.js v7** for visualizing effect size vs. statistical significance. Suitable as a portfolio or resume demo. Runs entirely in the browser with no backend.

## Overview

A volcano plot shows:

- **X-axis:** log₂ fold change (effect size)
- **Y-axis:** −log₁₀(p-value) (statistical significance)

Points in the upper left are significant and “down”-regulated; points in the upper right are significant and “up”-regulated. The plot uses **FDR (False Discovery Rate)** for significance, computed with the **Benjamini–Hochberg** procedure.

## Features

- **Synthetic data:** 1,200 points (id, log2FC, pval, FDR, −log10(p)) generated in JS; FDR computed via Benjamini–Hochberg.
- **Thresholds:**
  - **log₂ FC threshold** (slider 0–3): vertical lines at ±threshold.
  - **FDR threshold** (slider 0–0.2): horizontal line at the −log10(p) value that corresponds to that FDR level.
- **Color categories:**
  - **Significant up:** log2FC ≥ FC threshold and FDR ≤ FDR threshold (green).
  - **Significant down:** log2FC ≤ −FC threshold and FDR ≤ FDR threshold (red).
  - **Not significant:** otherwise (gray).
- **Interactivity:**
  - **Tooltip** on hover: id, log2FC, pval, FDR.
  - **Click** a point to pin its label; click again to unpin. Click on empty plot to clear pins and selection.
  - **Box selection:** drag on the plot to select points; selection count and **Export selected CSV** button appear.
- **Controls (left panel):**
  - FC threshold slider, FDR threshold slider.
  - Top N significant labels (number input + “Show labels” toggle).
  - Search box: type an id and click **Go** to highlight and zoom to that point.
  - Optional Y-axis log scale toggle (default remains −log10(p)).
  - **Regenerate data** to create a new synthetic dataset (new random seed).
- **Responsive** layout and **smooth transitions** when changing thresholds.
- **Accessibility:** keyboard focus for controls, ARIA labels on sliders and search.

## How thresholds work

- **FC threshold:** Only points with |log2FC| ≥ this value can be “significant” in the colored sense; the vertical lines mark the boundaries.
- **FDR threshold:** Points with FDR ≤ this value are considered significant. The **horizontal line** is drawn at the p-value that corresponds to this FDR level (i.e., the −log10 of that p-value). So points **above** the horizontal line have p-value below that cutoff; together with the FDR condition they drive the green/red coloring.

### How FDR is computed

FDR is computed with the **Benjamini–Hochberg (BH) procedure**:

1. **Sort p-values**  
   `p_(1) ≤ p_(2) ≤ … ≤ p_(n)`

2. **Compute adjusted values**  
   For each rank `k`, compute:  
   `q_k = min(1, p_(k) * n / k)`

3. **Assign adjusted values and threshold**  
   The adjusted value `q_k` is assigned back to the original index of `p_(k)`, with monotonicity enforced using the standard BH step-up procedure.  

The horizontal threshold line is set to match the p-value cutoff implied by the chosen FDR level for this dataset.

## How to run locally

No build step. Serve the project over HTTP (required for ES modules):

```bash
# From the project root:
python -m http.server 8080
```

Then open **http://localhost:8080** in a browser.

Other options:

- **Node:** `npx serve .` then open the URL shown.
- **PHP:** `php -S localhost:8080`

Do not open `index.html` as a file (e.g. `file:///...`); the D3 import will fail.

## Deploy to GitHub Pages

2. **Settings → Pages** → Source: **Deploy from a branch**.
3. Branch: **main** (or your default), folder: **/ (root)**.
4. Save. The site will be at `https://mrkchoe.github.io/interactive-volcano-plot/`.

## Performance

- **Rendering:** The demo uses **SVG** (one `<g>` of circles, no heavy filters). For ~1,200 points this stays smooth on modern devices.
- **Tradeoff:** For much larger datasets (e.g. 10k+ points), consider switching to **canvas** (e.g. draw circles in a canvas layer) or downsampling for the main view. This version is tuned for portfolio-sized data and avoids a build step; see the code for the single group of circles and minimal DOM.
- **Transitions:** Threshold and point updates use short D3 transitions (≈200–250 ms) for recolor and movement.

## File structure

```
/
├── index.html    # Entry point, structure, controls
├── styles.css    # Layout, theme, controls, tooltip
├── app.js        # D3 plot, BH FDR, data generation, interactivity
├── data/         # Optional: add CSV here if you switch from synthetic data
└── README.md     # This file
```



