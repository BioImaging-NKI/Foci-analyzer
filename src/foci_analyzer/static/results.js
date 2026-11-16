// The Results page: graphs of the tables of the last analysis, tests between conditions, a list per image
// and the cells to check. A click on a point or a row opens that cell or image in the viewer.

import { $, el, api, status, fmt, pvalue, helpIcon, isDark, store } from "./util.js";

const PALETTE = ["#3b82c4", "#e07b39", "#3a9d6e", "#b3589a", "#8a8f2a", "#5c6bc0", "#c0504d", "#2a9aa8"];
const SKIP = new Set(["cell", "frame", "focus", "channel", "x", "y", "z", "track", "parent", "voxels"]);
const HELP = {
  table: "Per cell: one value per nucleus or cell (for example the number of foci). Per focus: one value per focus (for example its size).",
  measure: "What is shown. foci chN = number of foci of channel N in the cell; total and mean intensity are above the local background; fraction chA with chB = share of the A foci that colocalize with a B focus; pearson and manders are intensity-based colocalization inside the cell; height = peak of the focus in threshold units.",
  coloc: "Per pair of foci channels: the share of the foci of one channel that colocalize with a focus of the other, per cell (boxes) with the mean per image (dots). Grey: the share expected when the foci of the other channel lay at random places in the same cell (the chance control). A share clearly above the grey box means real colocalization.",
  plot: "Points and box: every cell is a point, the box shows median and quartiles, the large dots are the means per image. Histogram: the distribution per condition. Scatter: two measures against each other. Over time: the mean per frame of a time-lapse. Tracks: one line per tracked cell over time (needs tracking to be on), with the mean per condition as a thick line.",
  unit: "What counts as one independent value in the test. Cells of one image are not independent of each other, so a test per cell easily gives very small p-values. Per image or per folder (one folder = one experiment) tests the means, which is the careful choice. Per cell uses Mann-Whitney, means use the Welch t-test; all pairwise p-values are Holm-corrected.",
};
const R = { cells: null, foci: null, meta: {}, loaded: false };
const C = Object.assign({ table: "cells", measure: "", plot: "points", x: "", unit: "image", channel: 0, log: false, off: [] }, store("fa.results") || {});
let hooks = {};
const put = (node, ...kids) => node.replaceChildren(...kids.flat().filter((k) => k !== null && k !== undefined));

function columns(t) {  // split-oriented JSON → { name: array }
  if (!t) return null;
  const out = { n: t.data.length };
  t.columns.forEach((name, j) => { out[name] = t.data.map((row) => row[j]); });
  return out;
}
const nice = (name) => name.replaceAll("_", " ").replace(/\bch(\d)/g, "ch $1");
function table() { return C.table === "foci" ? R.foci : R.cells; }
function measures(t) { return Object.keys(t).filter((k) => k !== "n" && !SKIP.has(k) && t[k].some((v) => typeof v === "number")); }
function conditions() { return (R.meta.conditions || []).filter((c) => R.cells.condition.includes(c)).concat([...new Set(R.cells.condition)].filter((c) => !(R.meta.conditions || []).includes(c))); }
function rows(t) {  // indices of the rows that are shown
  const out = [];
  for (let i = 0; i < t.n; i++) if (!C.off.includes(t.condition[i]) && (C.table !== "foci" || !C.channel || t.channel[i] === C.channel)) out.push(i);
  return out;
}

export function initResults(h) { hooks = h; document.addEventListener("themechange", () => { if (R.loaded && !$("page-results").hidden) draw(); }); }

export async function loadResults() {
  try {
    const r = await (await api("/api/results")).json();
    R.cells = columns(r.cells); R.foci = columns(r.foci); R.meta = r.meta || {}; R.loaded = true;
    render();
  } catch (e) { status(e.message, true); }
}

function render() {
  const side = $("results-side"), banner = $("results-banner");
  if (!R.cells) {
    side.replaceChildren(el("h2", { textContent: "Results" }), el("p", { class: "hint", textContent: "There are no results yet. Choose 'Analyze all images' on the Analyze page." }));
    for (const id of ["plot", "stats-card", "images-card", "check-card"]) $(id).replaceChildren();
    banner.hidden = true; return;
  }
  const notes = [];
  if (!R.meta.current) notes.push("The settings have changed since this analysis: the results below belong to the earlier settings.");
  for (const e of R.meta.errors || []) notes.push(e);
  banner.hidden = !notes.length; banner.textContent = notes.join("\n"); banner.style.whiteSpace = "pre-line";

  if (C.table === "foci" && !R.foci) C.table = "cells";
  const t = table(), ms = measures(t);
  if (!ms.includes(C.measure)) C.measure = ms.find((m) => m.startsWith("foci_ch")) || ms[0];
  if (!ms.includes(C.x)) C.x = ms.find((m) => m === "area") || ms[0];
  const timelapse = R.cells.frame.some((f) => f > 1), coloc = C.table === "cells" && pairs().length > 0;
  const tracked = timelapse && C.table === "cells" && !!R.cells.track && R.cells.track.some((v) => typeof v === "number" && v > 0);
  if (C.plot === "time" && !timelapse) C.plot = "points";
  if (C.plot === "tracks" && !tracked) C.plot = timelapse ? "time" : "points";
  if (C.plot === "coloc" && !coloc) C.plot = "points";
  if (C.plot === "coloc" && !/^fraction_/.test(C.measure)) C.measure = pairs()[0][0];  // the table below then tests that share
  const select = (key, options, help) => {
    const s = el("select", { onchange: (e) => { C[key] = key === "channel" ? +e.target.value : e.target.value; changed(); } },
      options.map(([v, text]) => el("option", { value: v, textContent: text })));
    s.value = C[key]; return [s, helpIcon(help)];
  };
  const field = (label, [input, help]) => el("div", { class: "field" }, el("label", { textContent: label }), help, el("div", { class: "control" }, input));
  const channels = R.foci ? [...new Set(R.foci.channel)].sort() : [];
  put(side,
    el("h2", { textContent: "Results" }),
    el("p", { class: "hint small", textContent: `${R.meta.when || ""} · ${R.cells.n} cells${R.foci ? ` · ${R.foci.n} foci` : ""}` }),
    field("Table", select("table", [["cells", "Per cell"], ...(R.foci ? [["foci", "Per focus"]] : [])], HELP.table)),
    C.table === "foci" && channels.length > 1 ? field("Channel", select("channel", [[0, "all"], ...channels.map((c) => [c, `channel ${c}`])], "")) : null,
    field("Measure", select("measure", ms.map((m) => [m, nice(m)]), HELP.measure)),
    field("Graph", select("plot", [["points", "Points and box"], ["violin", "Violin"], ["histogram", "Histogram"], ["scatter", "Scatter"], ...(timelapse ? [["time", "Over time"]] : []), ...(tracked ? [["tracks", "Tracks"]] : []), ...(coloc ? [["coloc", "Colocalization"]] : [])], C.plot === "coloc" ? HELP.coloc : HELP.plot)),
    C.plot === "scatter" ? field("Against", select("x", ms.map((m) => [m, nice(m)]), "The measure on the horizontal axis.")) : null,
    field("Test per", select("unit", [["cell", C.table === "foci" ? "focus" : "cell"], ["image", "image"], ["folder", "folder"]], HELP.unit)),
    el("label", { class: "small" }, el("input", { type: "checkbox", checked: C.log, onchange: (e) => { C.log = e.target.checked; changed(); } }), " logarithmic axis"),
    el("h2", { textContent: "Conditions", style: "margin-top:12px" }),
    el("div", { class: "cond-list" }, conditions().map((c, i) => el("label", {},
      el("input", { type: "checkbox", checked: !C.off.includes(c), onchange: (e) => { C.off = e.target.checked ? C.off.filter((x) => x !== c) : [...C.off, c]; changed(); } }),
      el("span", { style: `color:${PALETTE[i % PALETTE.length]}`, textContent: " ● " }), c))),
    el("div", { class: "row", style: "margin-top:12px" },
      el("button", { textContent: "Save graph data", title: "The values in this graph as a CSV file, one row per cell or focus", onclick: saveGraphData }),
      el("button", { textContent: "Open results folder", title: R.meta.folder || "", onclick: async () => { try { await api("/api/results/open", {}); } catch (e) { status(e.message, true); } } }),
      helpIcon("The folder holds cells.csv (one row per cell), foci.csv (one row per focus), settings.toml (the settings of this analysis) and run.json. The tables open in Excel, R or Python.")));
  draw(); drawStats(); drawImages(); drawCheck();
}
function changed() { store("fa.results", C); render(); }

// the pairs of foci channels with a colocalization share per cell: [column, column of the chance control or null, name]
function pairs() {
  return Object.keys(R.cells).filter((k) => /^fraction_ch\d+_with_ch\d+$/.test(k)).map((k) => {
    const [, a, b] = k.match(/^fraction_ch(\d+)_with_ch(\d+)$/), chance = `chance_ch${a}_with_ch${b}`;
    return [k, chance in R.cells ? chance : null, `ch ${a} foci with a ch ${b} focus`];
  });
}
function drawColoc(t, idx, conds, color, ink, grid, layout) {
  const traces = [], names = pairs();
  const perImage = (col, c) => {  // mean per image of the cells of condition c
    const m = new Map();
    for (const i of idx) if (t.condition[i] === c && typeof t[col][i] === "number") { const key = `${t.folder[i]}|${t.image[i]}`, a = m.get(key) || [0, 0, t.image[i]]; a[0] += t[col][i]; a[1]++; m.set(key, a); }
    return [...m.values()];
  };
  conds.forEach((c) => {
    const x = [], y = [], cx = [], cy = [], mx = [], my = [], mt = [];
    for (const [col, chance, name] of names) {
      for (const i of idx) if (t.condition[i] === c) {
        if (typeof t[col][i] === "number") { x.push(name); y.push(t[col][i]); }
        if (chance && typeof t[chance][i] === "number") { cx.push(name); cy.push(t[chance][i]); }
      }
      for (const m of perImage(col, c)) { mx.push(name); my.push(m[0] / m[1]); mt.push(`${m[2]} (mean of ${m[1]} cells)`); }
    }
    traces.push({ type: "box", x, y, name: c, legendgroup: c, offsetgroup: c, boxpoints: false, line: { color: color(c) }, hoverinfo: "y" });
    if (cx.length) traces.push({ type: "box", x: cx, y: cy, name: `${c}: by chance`, legendgroup: c, offsetgroup: `${c} chance`, boxpoints: false,
      line: { color: "#9aa0ad", width: 1.2 }, fillcolor: "rgba(154,160,173,.25)", hoverinfo: "y" });
    traces.push({ type: "box", x: mx, y: my, text: mt, name: `${c}: images`, legendgroup: c, offsetgroup: c, showlegend: false, boxpoints: "all", pointpos: 0, jitter: 0.5,
      line: { color: "rgba(0,0,0,0)" }, fillcolor: "rgba(0,0,0,0)", marker: { color: color(c), size: 9, line: { color: ink, width: 1.2 } }, hovertemplate: "%{text}<br>%{y:.3g}<extra></extra>" });
  });
  Object.assign(layout, { boxmode: "group", showlegend: true, legend: { orientation: "h", y: -0.18 } });
  layout.yaxis = { gridcolor: grid, zeroline: false, range: [-0.02, 1.02], title: { text: "share of the foci that colocalize" } };
  return traces;
}

// ---- the graph
function draw() {
  const t = table(), idx = rows(t), conds = conditions().filter((c) => !C.off.includes(c)), all = conditions();
  const color = (c) => PALETTE[all.indexOf(c) % PALETTE.length], y = t[C.measure];
  const dark = isDark(), ink = dark ? "#e6e8ee" : "#1d2330", grid = dark ? "#2c323d" : "#e3e6ec";
  const by = new Map(conds.map((c) => [c, []]));
  for (const i of idx) if (typeof y[i] === "number" && by.has(t.condition[i])) by.get(t.condition[i]).push(i);
  const hover = (i) => `${t.image[i]}${t.frame && t.frame[i] > 1 ? ` · frame ${t.frame[i]}` : ""} · cell ${t.cell[i]}`;
  const traces = [], layout = {
    paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)", font: { color: ink, family: "system-ui, Segoe UI, sans-serif", size: 13 },
    margin: { l: 60, r: 20, t: 20, b: 50 }, showlegend: false, hovermode: "closest",
    xaxis: { gridcolor: grid, zeroline: false }, yaxis: { gridcolor: grid, zeroline: false, title: { text: nice(C.measure) }, type: C.log ? "log" : "linear" },
  };
  if (C.plot === "coloc") traces.push(...drawColoc(t, idx, conds, color, ink, grid, layout));
  else if (C.plot === "points" || C.plot === "violin") {
    conds.forEach((c, k) => {
      const ii = by.get(c), vals = ii.map((i) => y[i]);
      if (C.plot === "violin") traces.push({ type: "violin", x0: k, y: vals, name: c, line: { color: color(c) }, box: { visible: true }, meanline: { visible: true }, points: false, hoverinfo: "y" });
      else {
        const step = Math.max(1, Math.ceil(ii.length / 4000)), shown = ii.filter((_, j) => j % step === 0);
        traces.push({ type: "scatter", mode: "markers", x: shown.map((_, j) => k + ((j * 0.6180339) % 1 - 0.5) * 0.5), y: shown.map((i) => y[i]), customdata: shown,
          text: shown.map(hover), hovertemplate: "%{text}<br>%{y}<extra></extra>", marker: { color: color(c), size: 4, opacity: 0.45 } });
        traces.push({ type: "box", x0: k, y: vals, boxpoints: false, width: 0.55, line: { color: ink, width: 1.2 }, fillcolor: "rgba(0,0,0,0)", hoverinfo: "y" });
      }
      const means = new Map();  // the mean per image, the large dots
      for (const i of ii) { const key = `${t.folder[i]}|${t.image[i]}`, m = means.get(key) || [0, 0, t.image[i]]; m[0] += y[i]; m[1]++; means.set(key, m); }
      const mm = [...means.values()];
      traces.push({ type: "scatter", mode: "markers", x: mm.map((_, j) => k + (mm.length > 1 ? (j / (mm.length - 1) - 0.5) * 0.3 : 0)), y: mm.map((m) => m[0] / m[1]),
        text: mm.map((m) => `${m[2]} (mean of ${m[1]})`), hovertemplate: "%{text}<br>%{y:.3g}<extra></extra>",
        marker: { color: color(c), size: 10, line: { color: ink, width: 1.5 } } });
    });
    Object.assign(layout.xaxis, { tickvals: conds.map((_, k) => k), ticktext: conds, range: [-0.6, conds.length - 0.4] });
  } else if (C.plot === "histogram") {
    conds.forEach((c) => traces.push({ type: "histogram", x: by.get(c).map((i) => y[i]), name: c, opacity: 0.6, histnorm: "probability", marker: { color: color(c) } }));
    Object.assign(layout, { barmode: "overlay", showlegend: true });
    layout.xaxis.title = { text: nice(C.measure) }; layout.yaxis = { gridcolor: grid, title: { text: "fraction" } };
  } else if (C.plot === "scatter") {
    const x = t[C.x];
    conds.forEach((c) => { const ii = by.get(c).filter((i) => typeof x[i] === "number");
      traces.push({ type: "scattergl", mode: "markers", x: ii.map((i) => x[i]), y: ii.map((i) => y[i]), customdata: ii, name: c, text: ii.map(hover),
        hovertemplate: "%{text}<br>%{x:.3g}, %{y:.3g}<extra></extra>", marker: { color: color(c), size: 5, opacity: 0.55 } }); });
    layout.showlegend = true; layout.xaxis.title = { text: nice(C.x) };
  } else if (C.plot === "tracks") {
    conds.forEach((c) => {
      const tracks = new Map(), per = new Map();
      for (const i of by.get(c)) {
        if (!(t.track[i] > 0)) continue;
        const key = `${t.path[i]}|${t.track[i]}`, a = tracks.get(key) || []; a.push(i); tracks.set(key, a);
        const f = per.get(t.frame[i]) || []; f.push(y[i]); per.set(t.frame[i], f);
      }
      [...tracks.values()].slice(0, 400).forEach((ii) => {
        ii.sort((a, b) => t.frame[a] - t.frame[b]);
        traces.push({ type: "scatter", mode: "lines", x: ii.map((i) => t.frame[i]), y: ii.map((i) => y[i]), customdata: ii, showlegend: false, opacity: 0.5,
          text: ii.map((i) => `${t.image[i]} · track ${t.track[i]}`), hovertemplate: "%{text}<br>frame %{x}: %{y:.3g}<extra></extra>", line: { color: color(c), width: 1 } });
      });
      const fr = [...per.keys()].sort((a, b) => a - b);
      traces.push({ type: "scatter", mode: "lines", x: fr, y: fr.map((f) => per.get(f).reduce((a, b) => a + b, 0) / per.get(f).length), name: c, line: { color: color(c), width: 3.5 } });
    });
    layout.showlegend = true; layout.xaxis.title = { text: "frame" };
  } else {
    conds.forEach((c) => {
      const per = new Map();
      for (const i of by.get(c)) { const a = per.get(t.frame[i]) || []; a.push(y[i]); per.set(t.frame[i], a); }
      const fr = [...per.keys()].sort((a, b) => a - b), mean = fr.map((f) => per.get(f).reduce((a, b) => a + b, 0) / per.get(f).length);
      const sem = fr.map((f, j) => { const v = per.get(f); return v.length > 1 ? Math.sqrt(v.reduce((a, b) => a + (b - mean[j]) ** 2, 0) / (v.length - 1) / v.length) : 0; });
      traces.push({ type: "scatter", mode: "lines+markers", x: fr, y: mean, name: c, error_y: { type: "data", array: sem, width: 0 }, line: { color: color(c) } });
    });
    layout.showlegend = true; layout.xaxis.title = { text: "frame" };
  }
  const plot = $("plot");
  window.Plotly.react(plot, traces, layout, { displaylogo: false, responsive: true, modeBarButtonsToRemove: ["lasso2d", "select2d"] });
  plot.removeAllListeners && plot.removeAllListeners("plotly_click");
  plot.on("plotly_click", (ev) => { const p = ev.points[0]; if (p && Number.isInteger(p.customdata)) openRow(t, p.customdata); });
}
// the rows behind the graph as a CSV file
function saveGraphData() {
  const t = table(), idx = rows(t);
  const cols = ["condition", "folder", "image", "frame", "cell", ...(t.focus ? ["focus", "channel"] : []), ...(t.track ? ["track"] : []),
    ...(C.plot === "coloc" ? pairs().flatMap(([a, b]) => (b ? [a, b] : [a])) : [...(C.plot === "scatter" ? [C.x] : []), C.measure])];
  const cell = (v) => v === null || v === undefined ? "" : /[",\n]/.test(String(v)) ? `"${String(v).replaceAll('"', '""')}"` : String(v);
  const text = [cols.join(","), ...idx.map((i) => cols.map((c) => cell(t[c][i])).join(","))].join("\n");
  const a = el("a", { href: URL.createObjectURL(new Blob([text], { type: "text/csv" })), download: `foci_analyzer_${C.plot === "coloc" ? "colocalization" : C.measure}.csv` });
  document.body.append(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}
function pathOf(t, i) {
  if (t.path) return t.path[i];
  const c = R.cells, j = c.image.findIndex((im, k) => im === t.image[i] && c.folder[k] === t.folder[i]);
  return j >= 0 ? c.path[j] : null;
}
function openRow(t, i) {
  const path = pathOf(t, i), info = (R.meta.images || {})[path], img = (R.meta.settings || {}).image || {};
  if (!path) return;
  if (!info || typeof t.x[i] !== "number") return hooks.openCell(path, t.frame[i]);
  const k = (img.xy_binning || 1) / info.pixel_size, b = img.crop_border || 0;
  hooks.openCell(path, t.frame[i], t.x[i] * k + b, t.y[i] * k + b, t.area ? 4 * Math.sqrt(t.area[i] / Math.PI) * k : 150);
}

// ---- tests
async function drawStats() {
  const card = $("stats-card");
  try {
    const s = await (await api("/api/stats", { table: C.table, measure: C.measure, unit: C.unit, channel: C.channel,
      conditions: conditions().filter((c) => !C.off.includes(c)) })).json();
    const what = { cell: C.table === "foci" ? "foci" : "cells", image: "images", folder: "folders" }[C.unit];
    put(card, el("h2", {}, `${nice(C.measure)} per condition `, helpIcon(HELP.unit)),
      el("table", {}, el("thead", {}, el("tr", {}, ["Condition", `n (${what})`, "Mean", "SD", "Median"].map((h, j) => el("th", { class: j ? "num" : "", textContent: h })))),
        el("tbody", {}, s.summary.map((r) => el("tr", {}, el("td", { textContent: r.condition }), [r.n, r.mean, r.sd, r.median].map((v) => el("td", { class: "num", textContent: fmt(v) })))))),
      s.pairs.length ? el("table", { style: "margin-top:12px" },
        el("thead", {}, el("tr", {}, [`${s.test}`, "p", "p (Holm)"].map((h, j) => el("th", { class: j ? "num" : "", textContent: h })))),
        el("tbody", {}, s.pairs.map((p) => el("tr", {}, el("td", { textContent: `${p.a}  vs  ${p.b}` }), el("td", { class: "num", textContent: pvalue(p.p) }),
          el("td", { class: "num", textContent: pvalue(p.p_holm), style: p.p_holm !== null && p.p_holm < 0.05 ? "font-weight:650" : "" }))))) : null,
      el("p", { class: "hint small", style: "margin:8px 0 0", textContent: [s.overall ? `${s.overall.test}: p = ${pvalue(s.overall.p)}.` : "", s.note].filter(Boolean).join(" ") }));
  } catch (e) { card.replaceChildren(el("p", { class: "hint", textContent: e.message })); }
}

// ---- per image, and the cells to check
function drawImages() {
  const c = R.cells, counts = Object.keys(c).filter((k) => /^foci_ch\d+$/.test(k)), per = new Map();
  for (let i = 0; i < c.n; i++) {
    const m = per.get(c.path[i]) || { i, n: 0, sum: counts.map(() => 0) };
    m.n++; counts.forEach((k, j) => { m.sum[j] += c[k][i] || 0; }); per.set(c.path[i], m);
  }
  $("images-card").replaceChildren(el("h2", {}, "Images ", helpIcon("One row per image: the number of cells and the mean number of foci per cell. Click a row to open the image in the viewer.")),
    el("table", {}, el("thead", {}, el("tr", {}, ["Condition", "Image", "Cells", ...counts.map((k) => `${nice(k)} per cell`)].map((h, j) => el("th", { class: j > 1 ? "num" : "", textContent: h })))),
      el("tbody", {}, [...per.entries()].map(([path, m]) => el("tr", { class: "click", onclick: () => hooks.openCell(path) },
        el("td", { textContent: c.condition[m.i] }), el("td", { textContent: c.image[m.i] }), el("td", { class: "num", textContent: m.n }),
        m.sum.map((s) => el("td", { class: "num", textContent: fmt(s / m.n) })))))));
}
function drawCheck() {
  const c = R.cells, list = [];
  for (let i = 0; i < c.n && list.length < 200; i++) if (c.check && c.check[i]) list.push(i);
  $("check-card").replaceChildren(el("h2", {}, `Cells to check (${list.length}${list.length === 200 ? "+" : ""}) `,
    helpIcon("Cells with something that may make their numbers less reliable, for example saturated foci or a position on the image edge. They are never left out automatically. Click a row to look at the cell.")),
    list.length ? el("table", {}, el("tbody", {}, list.map((i) => el("tr", { class: "click", onclick: () => openRow(c, i) },
      el("td", { textContent: c.image[i] }), el("td", { textContent: `cell ${c.cell[i]}` }), el("td", { textContent: c.check[i] })))))
      : el("p", { class: "hint", style: "margin:0", textContent: "None." }));
}
