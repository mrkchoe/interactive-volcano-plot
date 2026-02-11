/**
 * Interactive Volcano Plot — D3.js v7
 * Effect size (log2FC) vs statistical significance (-log10 p-value).
 */

import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7/+esm";

// --- Constants ---
const DEFAULT_FC_THRESHOLD = 1;
const DEFAULT_FDR_THRESHOLD = 0.05;
const DEFAULT_TOP_N = 10;
const N_POINTS = 1200;
const MARGIN = { top: 24, right: 24, bottom: 40, left: 48 };

// --- Benjamini–Hochberg FDR ---
/**
 * Computes FDR using the Benjamini–Hochberg procedure.
 * @param {number[]} pvals - Raw p-values (will not be mutated)
 * @returns {number[]} FDR values in same order as pvals
 */
function benjaminiHochberg(pvals) {
  const n = pvals.length;
  const indexed = pvals.map((p, i) => ({ p, i })).sort((a, b) => a.p - b.p);
  const fdr = new Array(n);
  for (let k = 0; k < n; k++) {
    const q = (indexed[k].p * n) / (k + 1);
    fdr[indexed[k].i] = Math.min(1, q);
  }
  return fdr;
}

// --- Seeded RNG (simple LCG) ---
function createRng(seed) {
  let s = seed;
  return function () {
    s = (s * 1664525 + 1013904223) >>> 0;
    return s / 0xffffffff;
  };
}

// --- Synthetic data (use real gene symbols so UniProt hover works) ---
const GENE_SYMBOLS = [
  "TP53", "BRCA1", "EGFR", "MYC", "AKT1", "PTEN", "KRAS", "ERBB2", "VEGFA", "IL6",
  "TNF", "MAPK1", "JUN", "FOS", "STAT3", "NFKB1", "CDKN1A", "BCL2", "CASP3", "ESR1",
  "AR", "INS", "INSR", "IGF1", "CTNNB1", "APC", "SMAD4", "TGFB1", "CDK1", "CCND1",
  "RB1", "E2F1", "MDM2", "CDKN2A", "GAPDH", "ACTB", "HSP90AA1", "HSPA8", "TUBB", "LMNA",
  "SOD1", "CAT", "GPX1", "NFE2L2", "HIF1A", "VHL", "MTOR", "PIK3CA", "GSK3B", "NOTCH1",
  "WNT1", "DVL1", "AXIN1", "LEF1", "TCF7L2", "MYCN", "FLT1", "KDR", "PDGFRA", "MET",
  "RET", "BRAF", "RAF1", "MAP2K1", "MAPK3", "ELK1", "CREB1", "ATF2", "JAK2", "SOCS1",
  "IL10", "IFNG", "CD4", "CD8A", "CD19", "CD34", "KIT", "FLT3", "NPM1", "CEBPA",
  "RUNX1", "GATA1", "TPO", "EPO", "VEGFB", "FGF2", "PDGFA", "EGF", "TGFB2", "BMP4",
  "WNT3A", "SHH", "DLL1", "JAG1", "HES1", "HEY1", "SNAI1", "TWIST1", "ZEB1", "CDH1",
  "VIM", "FN1", "COL1A1", "MMP2", "MMP9", "TIMP1", "SERPINE1", "PLAU", "CXCL12", "CCL2",
  "IL1B", "IL8", "COX2", "PTGS2", "NOS2", "ARG1", "IDO1", "CD274", "PDCD1", "CTLA4",
  "CD80", "CD86", "IL2", "IL12A", "TGFB3", "BMP2", "WNT5A", "FZD1", "LRP5", "DKK1",
];

function generateData(seed = Date.now()) {
  const rng = createRng(seed);
  const data = [];
  for (let i = 0; i < N_POINTS; i++) {
    const geneSymbol = GENE_SYMBOLS[i % GENE_SYMBOLS.length];
    const id = `${geneSymbol}_${i + 1}`;
    const log2FC = (rng() - 0.5) * 6;
    const z = Math.abs(log2FC) + rng() * 2;
    const pval = Math.max(1e-20, 2 * (1 - normalCdf(Math.abs(z))));
    data.push({ id, geneSymbol, log2FC, pval });
  }
  const pvals = data.map((d) => d.pval);
  const fdr = benjaminiHochberg(pvals);
  data.forEach((d, i) => {
    d.fdr = fdr[i];
    d.negLog10P = -Math.log10(d.pval);
  });
  return data;
}

function normalCdf(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.SQRT2;
  const t = 1.0 / (1.0 + p * x);
  const y = 1.0 - ((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return 0.5 * (1.0 + sign * y);
}

// --- Classification ---
function getCategory(d, fcThreshold, fdrThreshold) {
  if (d.fdr <= fdrThreshold && d.log2FC >= fcThreshold) return "sig_up";
  if (d.fdr <= fdrThreshold && d.log2FC <= -fcThreshold) return "sig_down";
  return "not_sig";
}

const COLOR = { sig_up: "#3fb950", sig_down: "#f85149", not_sig: "#484f58" };

// --- UniProt (https://www.uniprot.org/) ---
const UNIPROT_CACHE = new Map();
const UNIPROT_DESC_MAX = 220;

function fetchUniProtDescription(geneId) {
  const key = String(geneId).trim();
  if (!key) return Promise.resolve(null);
  if (UNIPROT_CACHE.has(key)) return Promise.resolve(UNIPROT_CACHE.get(key));

  const url =
    "https://rest.uniprot.org/uniprotkb/search?" +
    "query=" + encodeURIComponent("(gene:" + key + ")") +
    "&format=json&size=1" +
    "&fields=protein_name,gene_names,organism_name,cc_function";

  return fetch(url)
    .then((res) => (res.ok ? res.json() : null))
    .then((data) => {
      let proteinName = "";
      let description = "";
      const entry = data?.results?.[0];
      if (entry?.proteinDescription?.recommendedName?.fullName?.value) {
        proteinName = entry.proteinDescription.recommendedName.fullName.value;
      }
      const funcComment = entry?.comments?.find((c) => c.commentType === "FUNCTION");
      if (funcComment?.texts?.[0]?.value) {
        description = funcComment.texts[0].value;
        if (description.length > UNIPROT_DESC_MAX) {
          description = description.slice(0, UNIPROT_DESC_MAX).trim() + "…";
        }
      }
      const result = proteinName || description ? { proteinName, description } : null;
      UNIPROT_CACHE.set(key, result);
      return result;
    })
    .catch(() => {
      UNIPROT_CACHE.set(key, null);
      return null;
    });
}

// --- State ---
let state = {
  data: [],
  fcThreshold: DEFAULT_FC_THRESHOLD,
  fdrThreshold: DEFAULT_FDR_THRESHOLD,
  topN: DEFAULT_TOP_N,
  showLabels: true,
  pinned: new Set(),
  selected: new Set(),
  searchHighlightId: null,
  pvalAtFdr: null,
  zoomDomain: null, // { x: [min, max], y: [min, max] } when zoomed to a point
};

function pvalueAtFdrThreshold(data, fdrThreshold) {
  const sorted = [...data].sort((a, b) => a.pval - b.pval);
  const m = sorted.length;
  for (let k = 0; k < m; k++) {
    const q = (sorted[k].pval * m) / (k + 1);
    if (q > fdrThreshold) {
      return k === 0 ? sorted[0].pval : sorted[k - 1].pval;
    }
  }
  return sorted[m - 1].pval;
}

// --- DOM refs ---
let container, svg, gPlot, xScale, yScale;
let gPoints, gThresholds, gLabels, gSelectionBox;
let overlay, resizeObserver;

function initContainer() {
  container = d3.select("#plot-container");
  if (container.empty()) return;
  const rect = container.node().getBoundingClientRect();
  const width = rect.width > 0 ? rect.width : 800;
  const height = Math.max(400, rect.height || 500);
  svg = container
    .append("svg")
    .attr("viewBox", `0 0 ${width} ${height}`)
    .attr("width", "100%")
    .attr("height", "100%");
  gPlot = svg.append("g").attr("transform", `translate(${MARGIN.left},${MARGIN.top})`);
  overlay = gPlot.append("rect").attr("fill", "none").attr("pointer-events", "all").style("cursor", "crosshair");
  overlay.on("click", (e) => {
    if (e.defaultPrevented) return;
    state.pinned.clear();
    state.selected.clear();
    updateSelectionUI();
    redraw();
  });
  gThresholds = gPlot.append("g").attr("class", "thresholds");
  gPoints = gPlot.append("g").attr("class", "points");
  gLabels = gPlot.append("g").attr("class", "labels");
  gSelectionBox = gPlot.append("rect").attr("class", "selection-box").attr("visibility", "hidden");
  gPlot.append("g").attr("class", "x-axis");
  gPlot.append("g").attr("class", "y-axis");
  gPlot.append("text").attr("class", "axis-label x-label").attr("text-anchor", "middle").attr("fill", "#8b949e").text("log₂ FC");
  gPlot.append("text").attr("class", "axis-label y-label").attr("text-anchor", "middle").attr("fill", "#8b949e").text("−log₁₀(p)");
  resizeObserver = new ResizeObserver(() => {
    const w = container.node().getBoundingClientRect().width;
    const h = Math.max(400, container.node().getBoundingClientRect().height);
    svg.attr("viewBox", `0 0 ${w} ${h}`);
    updateScalesAndAxes(w, h);
    redraw();
  });
  resizeObserver.observe(container.node());
  overlay.on("mousedown", startBox);
  d3.select("body").on("mousemove", moveBox).on("mouseup", endBox);
  return { width, height };
}

function updateScalesAndAxes(svgWidth, svgHeight) {
  const width = svgWidth - MARGIN.left - MARGIN.right;
  const height = svgHeight - MARGIN.top - MARGIN.bottom;
  const extentX = d3.extent(state.data, (d) => d.log2FC);
  const padX = Math.max(0.5, (extentX[1] - extentX[0]) * 0.05);
  const domainX = state.zoomDomain ? state.zoomDomain.x : [extentX[0] - padX, extentX[1] + padX];
  xScale = d3.scaleLinear().domain(domainX).range([0, width]);
  const maxY = d3.max(state.data, (d) => d.negLog10P);
  const minY = 0;
  const domainY = state.zoomDomain ? state.zoomDomain.y : [minY, Math.max(maxY * 1.05, 2)];
  yScale = d3.scaleLinear().domain(domainY).range([height, 0]);
  gPlot.select(".x-axis").attr("transform", `translate(0,${height})`).call(d3.axisBottom(xScale).ticks(8));
  gPlot.select(".y-axis").call(d3.axisLeft(yScale).ticks(8));
  gPlot.select(".x-label").attr("x", width / 2).attr("y", height + 36);
  gPlot.select(".y-label").attr("x", -36).attr("y", height / 2).attr("transform", `rotate(-90, -36, ${height / 2})`);
  overlay.attr("width", width).attr("height", height);
  return { width, height };
}

function drawThresholdLines(width, height) {
  const fc = state.fcThreshold;
  const fdr = state.fdrThreshold;
  state.pvalAtFdr = pvalueAtFdrThreshold(state.data, fdr);
  const yVal = -Math.log10(state.pvalAtFdr);
  const clipY = Math.max(0, Math.min(height, yScale(yVal)));
  const lines = [
    { type: "v", x1: xScale(fc), y1: 0, x2: xScale(fc), y2: height },
    { type: "v", x1: xScale(-fc), y1: 0, x2: xScale(-fc), y2: height },
    { type: "h", x1: 0, y1: clipY, x2: width, y2: clipY },
  ];
  const sel = gThresholds.selectAll("line").data(lines);
  sel.join("line")
    .attr("stroke", "#8b949e")
    .attr("stroke-width", 1)
    .attr("stroke-dasharray", "4 2")
    .attr("opacity", 0.8)
    .attr("x1", (d) => d.x1)
    .attr("y1", (d) => d.y1)
    .attr("x2", (d) => d.x2)
    .attr("y2", (d) => d.y2)
    .transition()
    .duration(200)
    .attr("x1", (d) => d.x1)
    .attr("y1", (d) => d.y1)
    .attr("x2", (d) => d.x2)
    .attr("y2", (d) => d.y2);
}

function redraw() {
  if (state.data.length === 0) return;
  const el = container.node();
  const rect = el ? el.getBoundingClientRect() : { width: 0, height: 0 };
  const w = rect.width > 0 ? rect.width : 800;
  const h = rect.height > 0 ? Math.max(400, rect.height) : 500;
  const { width, height } = updateScalesAndAxes(w, h);
  drawThresholdLines(width, height);

  const fc = state.fcThreshold;
  const fdr = state.fdrThreshold;
  state.data.forEach((d) => {
    d._category = getCategory(d, fc, fdr);
  });

  const points = gPoints.selectAll("circle").data(state.data, (d) => d.id);
  points
    .join("circle")
    .attr("class", (d) => {
      let c = "point";
      if (state.pinned.has(d.id)) c += " pinned";
      if (state.selected.has(d.id)) c += " selected";
      if (state.searchHighlightId === d.id) c += " highlight-search";
      return c;
    })
    .attr("r", (d) => (state.selected.has(d.id) || state.pinned.has(d.id) ? 5 : 3.5))
    .attr("cx", (d) => xScale(d.log2FC))
    .attr("cy", (d) => yScale(d.negLog10P))
    .attr("fill", (d) => COLOR[d._category])
    .attr("stroke", "transparent")
    .attr("stroke-width", 2)
    .style("cursor", "pointer")
    .on("mouseenter", (e, d) => showTooltip(e, d))
    .on("mousemove", (e) => moveTooltip(e))
    .on("mouseleave", hideTooltip)
    .on("click", (e, d) => {
      e.preventDefault();
      if (state.selected.size > 0) return;
      if (state.pinned.has(d.id)) state.pinned.delete(d.id);
      else state.pinned.add(d.id);
      redraw();
    });
  points
    .transition()
    .duration(250)
    .attr("cx", (d) => xScale(d.log2FC))
    .attr("cy", (d) => yScale(d.negLog10P))
    .attr("fill", (d) => COLOR[d._category]);

  const topSignificant = state.data
    .filter((d) => d._category !== "not_sig")
    .sort((a, b) => a.pval - b.pval)
    .slice(0, state.topN);
  const toLabel = new Set(
    state.showLabels ? topSignificant.map((d) => d.id) : []
  );
  state.pinned.forEach((id) => toLabel.add(id));

  const labels = gLabels.selectAll("text").data(state.data.filter((d) => toLabel.has(d.id)), (d) => d.id);
  labels
    .join("text")
    .attr("class", "point-label")
    .attr("x", (d) => xScale(d.log2FC))
    .attr("y", (d) => yScale(d.negLog10P))
    .attr("dy", -10)
    .attr("text-anchor", "middle")
    .attr("font-size", "10px")
    .attr("fill", "#e6edf3")
    .text((d) => d.id)
    .clone(true)
    .lower()
    .attr("fill", "none")
    .attr("stroke", "#0f1419")
    .attr("stroke-width", 3);
}

function showTooltip(e, d) {
  const tip = d3.select("#tooltip");
  const baseHtml =
    `<div class="row"><span class="label">id</span> ${d.id}</div>` +
    `<div class="row"><span class="label">log2FC</span> ${d.log2FC.toFixed(3)}</div>` +
    `<div class="row"><span class="label">pval</span> ${d.pval.toExponential(2)}</div>` +
    `<div class="row"><span class="label">fdr</span> ${d.fdr.toExponential(2)}</div>` +
    `<div class="row uniprot-row"><span class="label">UniProt</span> <span class="uniprot-text">Loading…</span></div>`;
  tip.html(baseHtml).classed("visible", true).attr("aria-hidden", "false");
  moveTooltip(e);

  const geneId = d.geneSymbol != null ? d.geneSymbol : d.id;
  const textEl = tip.select(".uniprot-text");
  fetchUniProtDescription(geneId).then((info) => {
    if (!textEl.node() || !tip.classed("visible")) return;
    const currentId = tip.node().getAttribute("data-current-id");
    if (currentId !== geneId) return;
    if (info?.proteinName || info?.description) {
      let html = "";
      if (info.proteinName) html += `<strong>${info.proteinName}</strong>`;
      if (info.description) html += (html ? " — " : "") + info.description;
      textEl.html(html);
    } else {
      textEl.text("No description found.");
    }
  });
  tip.node().setAttribute("data-current-id", geneId);
}

function moveTooltip(e) {
  d3.select("#tooltip")
    .style("left", `${e.pageX + 12}px`)
    .style("top", `${e.pageY + 12}px`);
}

function hideTooltip() {
  d3.select("#tooltip").classed("visible", false).attr("aria-hidden", "true");
}

// --- Box selection ---
let boxStart = null;
function startBox(e) {
  if (e.button !== 0) return;
  const pt = d3.pointer(e, gPlot.node());
  boxStart = { x: pt[0], y: pt[1] };
  gSelectionBox
    .attr("x", pt[0])
    .attr("y", pt[1])
    .attr("width", 0)
    .attr("height", 0)
    .attr("visibility", "visible");
}

function moveBox(e) {
  if (!boxStart) return;
  const pt = d3.pointer(e, gPlot.node());
  const x = Math.min(boxStart.x, pt[0]);
  const y = Math.min(boxStart.y, pt[1]);
  const w = Math.abs(pt[0] - boxStart.x);
  const h = Math.abs(pt[1] - boxStart.y);
  gSelectionBox.attr("x", x).attr("y", y).attr("width", w).attr("height", h);
}

function endBox(e) {
  if (!boxStart || e.button !== 0) return;
  const pt = d3.pointer(e, gPlot.node());
  const x = Math.min(boxStart.x, pt[0]);
  const y = Math.min(boxStart.y, pt[1]);
  const w = Math.abs(pt[0] - boxStart.x);
  const h = Math.abs(pt[1] - boxStart.y);
  gSelectionBox.attr("visibility", "hidden");
  if (w > 4 && h > 4) {
    const [x0, x1] = [x, x + w].map((px) => xScale.invert(px));
    const [y0, y1] = [y, y + h].map((py) => yScale.invert(py));
    state.selected.clear();
    state.data.forEach((d) => {
      if (
        d.log2FC >= Math.min(x0, x1) &&
        d.log2FC <= Math.max(x0, x1) &&
        d.negLog10P >= Math.min(y0, y1) &&
        d.negLog10P <= Math.max(y0, y1)
      )
        state.selected.add(d.id);
    });
    updateSelectionUI();
    redraw();
  }
  boxStart = null;
}

// Box selection listeners attached in init() after overlay exists

function updateSelectionUI() {
  const info = document.getElementById("selection-info");
  const countEl = document.getElementById("selection-count");
  if (!info || !countEl) return;
  const n = state.selected.size;
  if (n === 0) {
    info.hidden = true;
  } else {
    info.hidden = false;
    countEl.textContent = n;
  }
}

function exportSelectedCsv() {
  if (state.selected.size === 0) return;
  const rows = state.data.filter((d) => state.selected.has(d.id));
  const header = "id,log2FC,pval,fdr,negLog10P\n";
  const body = rows
    .map((d) => `${d.id},${d.log2FC},${d.pval},${d.fdr},${d.negLog10P}`)
    .join("\n");
  const blob = new Blob([header + body], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "volcano_selected.csv";
  a.click();
  URL.revokeObjectURL(a.href);
}

// --- Controls ---
function bindControls() {
  const fcInput = document.getElementById("fc-threshold");
  const fdrInput = document.getElementById("fdr-threshold");
  const fcValue = document.getElementById("fc-value");
  const fdrValue = document.getElementById("fdr-value");
  const topNInput = document.getElementById("top-n");
  const showLabelsCb = document.getElementById("show-labels");
  const searchId = document.getElementById("search-id");
  const searchBtn = document.getElementById("search-btn");
  const regenerateBtn = document.getElementById("regenerate");
  const exportBtn = document.getElementById("export-csv");

  function updateFc() {
    const v = parseFloat(fcInput.value);
    state.fcThreshold = v;
    if (fcValue) fcValue.textContent = v.toFixed(1);
    fcInput.setAttribute("aria-valuenow", v);
    fcInput.setAttribute("aria-valuetext", v.toFixed(1));
    redraw();
  }
  function updateFdr() {
    const v = parseFloat(fdrInput.value);
    state.fdrThreshold = v;
    if (fdrValue) fdrValue.textContent = v.toFixed(3);
    fdrInput.setAttribute("aria-valuetext", v.toFixed(3));
    redraw();
  }

  if (fcInput) fcInput.addEventListener("input", updateFc);
  if (fdrInput) fdrInput.addEventListener("input", updateFdr);

  if (topNInput)
    topNInput.addEventListener("change", () => {
      state.topN = Math.max(0, Math.min(50, parseInt(topNInput.value, 10) || 10));
      redraw();
    });
  if (showLabelsCb)
    showLabelsCb.addEventListener("change", () => {
      state.showLabels = showLabelsCb.checked;
      redraw();
    });

  function doSearch() {
    const q = (searchId?.value || "").trim().toLowerCase();
    if (!q) {
      state.searchHighlightId = null;
      redraw();
      return;
    }
    const found = state.data.find((d) => d.id.toLowerCase().includes(q));
    if (found) {
      state.searchHighlightId = found.id;
      const pad = 1.2;
      state.zoomDomain = {
        x: [found.log2FC - pad, found.log2FC + pad],
        y: [Math.max(0, found.negLog10P - pad), found.negLog10P + pad],
      };
      redraw();
    } else {
      state.searchHighlightId = null;
      state.zoomDomain = null;
      redraw();
    }
  }
  if (searchBtn) searchBtn.addEventListener("click", doSearch);
  if (searchId) searchId.addEventListener("keydown", (e) => { if (e.key === "Enter") doSearch(); });
  const resetZoomBtn = document.getElementById("reset-zoom-btn");
  if (resetZoomBtn)
    resetZoomBtn.addEventListener("click", () => {
      state.zoomDomain = null;
      state.searchHighlightId = null;
      redraw();
    });

  const graphSizeInput = document.getElementById("graph-size");
  const graphSizeValue = document.getElementById("graph-size-value");
  const plotWrapper = document.querySelector(".plot-wrapper");
  if (graphSizeInput && plotWrapper) {
    function updateGraphSize() {
      const pct = parseInt(graphSizeInput.value, 10);
      const fraction = pct / 100;
      plotWrapper.style.setProperty("--graph-size", String(fraction));
      if (graphSizeValue) graphSizeValue.textContent = pct + "%";
      graphSizeInput.setAttribute("aria-valuenow", pct);
      graphSizeInput.setAttribute("aria-valuetext", pct + "%");
      if (typeof redraw === "function") redraw();
    }
    graphSizeInput.addEventListener("input", updateGraphSize);
    updateGraphSize();
  }

  if (regenerateBtn)
    regenerateBtn.addEventListener("click", () => {
      state.data = generateData();
      state.pinned.clear();
      state.selected.clear();
      state.searchHighlightId = null;
      state.zoomDomain = null;
      updateSelectionUI();
      redraw();
    });

  if (exportBtn) exportBtn.addEventListener("click", exportSelectedCsv);
}

// --- Resizable plot (drag handle to scale plot area) ---
const PANEL_MIN = 200;
const PANEL_MAX = 480;

function setupResizeHandle() {
  const handle = document.getElementById("resize-handle");
  const panel = document.querySelector(".panel");
  if (!handle || !panel) return;

  let startX = 0;
  let startWidth = 0;

  function move(e) {
    const dx = e.clientX - startX;
    let w = startWidth + dx;
    w = Math.max(PANEL_MIN, Math.min(PANEL_MAX, w));
    panel.style.width = `${w}px`;
    handle.setAttribute("aria-valuenow", String(Math.round(w)));
  }

  function stop() {
    document.removeEventListener("mousemove", move);
    document.removeEventListener("mouseup", stop);
    document.body.style.cursor = "";
    document.body.style.userSelect = "";
  }

  handle.addEventListener("mousedown", (e) => {
    if (e.button !== 0) return;
    e.preventDefault();
    startX = e.clientX;
    startWidth = panel.getBoundingClientRect().width;
    document.body.style.cursor = "col-resize";
    document.body.style.userSelect = "none";
    document.addEventListener("mousemove", move);
    document.addEventListener("mouseup", stop);
  });

  handle.setAttribute("tabindex", "0");
}

// --- Init ---
function init() {
  state.data = generateData();
  const dims = initContainer();
  if (dims) {
    const w = dims.width > 0 ? dims.width : 800;
    const h = dims.height > 0 ? dims.height : 500;
    updateScalesAndAxes(w, h);
    redraw();
  }
  bindControls();
  setupResizeHandle();
}

if (document.readyState === "loading") {
  document.addEventListener("DOMContentLoaded", init);
} else {
  init();
}
