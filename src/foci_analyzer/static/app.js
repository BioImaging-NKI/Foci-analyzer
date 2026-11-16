// Foci Analyzer window: pages, the experiment, the settings form (drawn from the schema), the viewer with
// the live preview, and the batch run. The Results page is in results.js, the viewer in viewer.js.

import { $, el, api, status, fmt, helpIcon, store } from "./util.js";
import { Viewer } from "./viewer.js";
import { initResults, loadResults } from "./results.js";

const S = { schema: [], settings: {}, conditions: [], images: [], remarks: [], output: "", models: {}, batch: {} };

// ---------------------------------------------------------------- pages
function showPage(name) {
  for (const b of $("tabs").children) b.classList.toggle("on", b.dataset.page === name);
  for (const p of document.querySelectorAll(".page")) p.hidden = p.id !== `page-${name}`;
  if (name === "analyze") viewer.requestDraw();
  if (name === "results") loadResults();
}
$("tabs").addEventListener("click", (e) => { if (e.target.dataset.page) showPage(e.target.dataset.page); });

// ---------------------------------------------------------------- experiment
function renderConditions() {
  const box = $("conditions"); box.replaceChildren();
  S.conditions.forEach((c, i) => {
    const folders = el("div");
    c.folders.forEach((f, j) => folders.append(el("div", { class: "folder" },
      el("input", { type: "text", value: f, placeholder: "Folder or image file", "aria-label": "Folder", oninput: (e) => { c.folders[j] = e.target.value; } }),
      el("button", { class: "link", textContent: "Remove", onclick: () => { c.folders.splice(j, 1); renderConditions(); } }))));
    box.append(el("div", { class: "condition" },
      el("div", { class: "row" },
        el("input", { type: "text", class: "name", value: c.name, placeholder: "Condition name", "aria-label": "Condition name", oninput: (e) => { c.name = e.target.value; } }),
        el("label", {}, "Files named ", el("input", { type: "text", class: "pattern", value: c.pattern, title: "Pattern for the file name, * = anything", oninput: (e) => { c.pattern = e.target.value; } })),
        el("label", {}, el("input", { type: "checkbox", checked: c.recursive, onchange: (e) => { c.recursive = e.target.checked; } }), " with sub-folders"),
        el("button", { class: "link", textContent: "Remove condition", onclick: () => { S.conditions.splice(i, 1); renderConditions(); } })),
      folders,
      el("div", { class: "row" }, el("button", { textContent: "Add folder", onclick: () => addPaths(c, "folder") }),
        el("button", { textContent: "Add images", onclick: () => addPaths(c, "files") }))));
  });
}
async function pick(kind) {
  try { return (await (await api("/api/pick", { kind })).json()).paths; }
  catch (e) { status(e.message, true); return null; }
}
async function addPaths(c, kind) {
  const paths = await pick(kind);
  if (paths === null) c.folders.push("");  // no dialog: an empty box to type or paste the folder in
  else c.folders.push(...paths);
  renderConditions();
}
function renderImages() {
  const body = $("images").tBodies[0]; body.replaceChildren();
  for (const im of S.images) body.append(el("tr", { class: "click", onclick: () => { showPage("analyze"); openImage(im.id); } },
    el("td", { textContent: im.name }), el("td", { textContent: im.condition }), el("td", { class: "path", textContent: im.folder, title: im.folder })));
  $("image-count").textContent = S.images.length ? `(${S.images.length})` : "";
  $("remarks").replaceChildren(...S.remarks.map((r) => el("li", { textContent: r })));
  const sel = $("image-select"); sel.replaceChildren(el("option", { value: "", textContent: "Choose an image…" }));
  for (const im of S.images) sel.append(el("option", { value: im.id, textContent: `${im.condition} · ${im.name}` }));
  sel.value = S.images.some((im) => im.id === V.id) ? V.id : "";
  $("output").value = S.output || "";
}
async function scan() {
  try {
    status("Looking for images…");
    Object.assign(S, await (await api("/api/experiment", { conditions: S.conditions, output: $("output").value })).json());
    renderConditions(); renderImages(); status(`${S.images.length} images`);
    return true;
  } catch (e) { status(e.message, true); return false; }
}
$("add-condition").onclick = () => { S.conditions.push({ name: `Condition ${S.conditions.length + 1}`, folders: [], pattern: "*", recursive: false }); renderConditions(); };
$("scan").onclick = scan;
$("output").onchange = scan;
$("pick-output").onclick = async () => { const p = await pick("folder"); if (p && p[0]) { $("output").value = p[0]; scan(); } };
async function quick(kind) {
  const paths = await pick(kind);
  if (!paths || !paths.length) return;
  const empty = S.conditions.every((c) => !c.folders.length);
  if (empty) S.conditions = [{ name: "Images", folders: paths, pattern: "*", recursive: false }];
  else {
    let c = S.conditions.find((x) => x.name === "Images");
    if (!c) { c = { name: "Images", folders: [], pattern: "*", recursive: false }; S.conditions.push(c); }
    c.folders.push(...paths);
  }
  if (await scan() && S.images.length) {
    const first = S.images.find((im) => paths.some((p) => im.path.startsWith(p))) || S.images[0];
    showPage("analyze"); openImage(first.id);
  }
}
$("quick-files").onclick = () => quick("files");
$("quick-folder").onclick = () => quick("folder");

// ---------------------------------------------------------------- settings form
let changeTimer = 0;
function settingsChanged(delay = 300) {
  clearTimeout(changeTimer);
  changeTimer = setTimeout(async () => {
    updateChannelColors();
    if (V.id && $("live").checked) return runPreview();  // the preview request also saves the settings
    try { S.settings = await (await api("/api/settings", { settings: S.settings })).json(); status("Settings saved"); }
    catch (e) { status(e.message, true); }
  }, delay);
}
function control(f, group, sectionKey, index) {
  const set = (v, delay) => { group[f.key] = v; settingsChanged(delay); if (f.type === "choice") renderSettings(); };
  const id = `f${Math.random().toString(36).slice(2)}`, v = group[f.key];
  const number = (e) => { const n = f.type === "int" ? parseInt(e.target.value, 10) : parseFloat(e.target.value); return Number.isFinite(n) ? n : null; };
  let input, range = null;
  if (f.type === "bool") input = el("input", { id, type: "checkbox", checked: v, onchange: (e) => set(e.target.checked) });
  else if (f.type === "choice") { input = el("select", { id, onchange: (e) => set(e.target.value) }, f.choices.map((c) => el("option", { value: c, textContent: c }))); input.value = v; }
  else if (f.type === "str") input = el("input", { id, type: "text", value: v, onchange: (e) => set(e.target.value) });
  else {
    input = el("input", { id, type: "number", value: v, min: f.min ?? "", max: f.max ?? "", step: f.step ?? (f.type === "int" ? 1 : "any"),
      onchange: (e) => { const n = number(e); if (n === null) e.target.value = group[f.key]; else { set(n); if (range) range.value = n; } } });
    if (f.slider) range = el("input", { type: "range", min: f.min, max: f.max, step: f.step, value: v, "aria-label": f.label,
      oninput: (e) => { const n = number(e); input.value = n; liveThreshold(sectionKey, index, f.key, n); set(n, 120); } });
  }
  return el("div", { class: `field${f.type === "bool" ? " check" : ""}${range ? " slider" : ""}` },
    el("label", { htmlFor: id, textContent: f.label }), helpIcon(f.help),
    el("div", { class: "control" }, input, f.unit ? el("span", { class: "unit", textContent: f.unit }) : null), range);
}
function fieldsOf(section, group, index) {
  const shown = section.fields.filter((f) => !f.when || Object.entries(f.when).every(([k, values]) => values.includes(group[k])));
  const basic = shown.filter((f) => !f.advanced), adv = shown.filter((f) => f.advanced);
  const open = store(`fa.adv.${section.key}`);
  return [...basic.map((f) => control(f, group, section.key, index)),
    adv.length ? el("details", { class: "advanced", open: !!open, ontoggle: (e) => store(`fa.adv.${section.key}`, e.target.open) },
      el("summary", { textContent: "Advanced" }), adv.map((f) => control(f, group, section.key, index))) : null];
}
function renderSettings() {
  const box = $("settings"); box.replaceChildren();
  for (const sec of S.schema) {
    const div = el("div", { class: "section" });
    if (!sec.repeat) { div.append(el("h3", { textContent: sec.title }), ...fieldsOf(sec, S.settings[sec.key], 0).filter(Boolean)); box.append(div); continue; }
    const list = S.settings[sec.key];
    div.append(el("h3", {}, `${sec.title}s`, el("button", { textContent: "Add", title: "Add a foci channel", onclick: () => {
      const g = Object.fromEntries(sec.fields.map((f) => [f.key, f.default]));
      if ("channel" in g) g.channel = Math.max(...list.map((x) => x.channel)) + 1;
      list.push(g); renderSettings(); settingsChanged(); } })));
    list.forEach((g, i) => div.append(el("div", { class: "group" },
      el("div", { class: "group-head" }, `${sec.title} ${i + 1}`,
        list.length > 1 ? el("button", { class: "link", textContent: "Remove", onclick: () => { list.splice(i, 1); renderSettings(); settingsChanged(); } }) : null),
      ...fieldsOf(sec, g, i).filter(Boolean))));
    box.append(div);
  }
}

// ---------------------------------------------------------------- viewer
const viewer = new Viewer($("canvas"), $("markers"), (url) => api(url));
const V = { id: "", info: null, mode: "mip", z: 0, t: 0, ch: [], preview: null, prevKey: "", cand: [], load: 0 };
const NUCLEI = [0, 0.5, 1], CELLS = [1, 0.5, 0], OTHER = [[0, 1, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]];
const OUTLINE = [[1, 0.88, 0.3], [0.3, 0.88, 1], [1, 0.6, 0.2], [0.8, 0.6, 1]];  // foci channel 1, 2, …
const hex = (c) => `#${c.map((v) => Math.round(v * 255).toString(16).padStart(2, "0")).join("")}`;
const unhex = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);

function channelColor(c) {
  const saved = (store("fa.colors") || {})[c];
  if (saved) return unhex(saved);
  if (c + 1 === S.settings.image.nuclei_channel) return NUCLEI;
  if (c + 1 === S.settings.image.cell_channel) return CELLS;
  const foci = S.settings.foci.map((f) => f.channel - 1), i = foci.indexOf(c);
  return OTHER[(i >= 0 ? i : foci.length + c) % OTHER.length];
}
function updateChannelColors() { V.ch.forEach((ch, c) => { ch.color = channelColor(c); }); renderChannels(); syncLayers(); }

async function openImage(id) {
  const load = ++V.load;
  $("image-select").value = id;
  if (!id) return;
  try {
    status("Reading image…");
    const info = await (await api(`/api/image/info?path=${encodeURIComponent(id)}`)).json();
    if (load !== V.load) return;
    Object.assign(V, { id, info, z: Math.floor(info.slices / 2), t: 0, preview: null, prevKey: "", cand: [] });
    const top = info.dtype === "uint8" ? 255 : 65535;
    V.ch = info.ranges.map(([lo, hi], c) => ({ on: true, lo, hi, auto: [lo, hi], top: Math.min(top, Math.max(hi * 3, lo + 10)), color: channelColor(c) }));
    if (info.slices === 1) V.mode = "mip";
    $("mode").hidden = info.slices === 1;
    $("z").max = info.slices - 1; $("z").value = V.z; $("t").max = info.frames - 1; $("t").value = 0;
    $("viewer-empty").hidden = true; $("cell-info").hidden = true;
    viewer.setImage(info.width, info.height); viewer.pixelSize = info.pixel_size; viewer.unit = info.unit;
    renderChannels(); syncMode(); syncLayers(); viewer.fit(); status("");
    $("preview-info").textContent = "";
    if ($("live").checked) runPreview();
  } catch (e) { status(e.message, true); }
}
function syncLayers() {
  if (!V.info) { viewer.layers = []; viewer.markers = []; viewer.requestDraw(); return; }
  const z = V.mode === "mip" ? "mip" : V.z, base = `path=${encodeURIComponent(V.id)}&t=${V.t}&z=${z}`;
  const region = (extra) => (level, x, y, w, h) => `${extra}&${base}&level=${level}&x=${x}&y=${y}&w=${w}&h=${h}`;
  const layers = V.ch.map((ch, c) => ({ type: "image", key: `i|${V.id}|${V.t}|${c}|${z}`, url: region(`/api/image/region?c=${c}`),
    lo: ch.lo, hi: ch.hi, color: ch.color, on: ch.on }));
  viewer.markers = [];
  if (V.preview) {
    const k = V.preview.key, pk = V.prevKey;
    layers.push({ type: "labels", key: `l|cells|${z}|${k}`, prevKey: pk && `l|cells|${z}|${pk}`, url: region("/api/overlay?kind=cells"),
      color: [1, 1, 1, 0.85], on: $("show-cells").checked });
    V.preview.channels.forEach((ch, i) => {
      const col = OUTLINE[i % OUTLINE.length];
      layers.push({ type: "labels", key: `l|foci${ch.channel}|${z}|${k}`, prevKey: pk && `l|foci${ch.channel}|${z}|${pk}`,
        url: region(`/api/overlay?kind=foci&ch=${ch.channel}`), color: [...col, 1], on: $("show-foci").checked });
      viewer.markers.push({ pts: V.cand[i], color: hex(col), threshold: ch.live ?? ch.threshold, on: $("show-markers").checked });
    });
  }
  viewer.layers = layers; viewer.zslice = V.mode === "mip" ? null : V.z; viewer.requestDraw();
}
viewer.onview = () => {
  const i = V.info; if (!i) return;
  const cal = i.pixel_size ? ` · ${fmt(i.pixel_size)} ${i.unit}/px` : " · not calibrated";
  $("canvas").title = ""; $("image-select").title = `${i.width} × ${i.height} px · ${i.channels} channels · ${i.slices} slices${cal} · zoom ${Math.round(viewer.view.s * 100)} %`;
};

// ---- preview: the analysis of the shown image with the current settings
let previewBusy = false, previewAgain = false;
async function runPreview() {
  if (!V.id) return;
  if (previewBusy) { previewAgain = true; return; }
  previewBusy = true; status("Analyzing…");
  const id = V.id, t = V.t;
  try {
    const r = await (await api("/api/preview", { path: id, t, settings: S.settings })).json();
    const cand = await Promise.all(r.channels.map(async (ch) =>
      new Float32Array(await (await api(`/api/candidates?path=${encodeURIComponent(id)}&t=${t}&ch=${ch.channel}`)).arrayBuffer())));
    if (id === V.id && t === V.t) {
      V.prevKey = V.preview ? V.preview.key : ""; V.preview = r; V.cand = cand;
      syncLayers(); previewInfo(); status("");
    }
  } catch (e) { status(e.message, true); }
  previewBusy = false;
  if (previewAgain) { previewAgain = false; runPreview(); }
}
function countAbove(pts, thr) {  // heights are sorted from high to low
  let lo = 0, hi = pts.length / 6;
  while (lo < hi) { const m = (lo + hi) >> 1; if (pts[m * 6 + 3] >= thr) lo = m + 1; else hi = m; }
  return lo;
}
function previewInfo() {
  const p = V.preview; if (!p) { $("preview-info").textContent = ""; return; }
  const parts = [el("b", { textContent: `${p.cells} ${p.cells === 1 ? "cell" : "cells"}` })];
  p.channels.forEach((ch, i) => {
    const live = ch.live !== undefined, n = live ? countAbove(V.cand[i], ch.live) : ch.foci;
    const thr = ch.height_unit === "× noise" ? `noise ${fmt(ch.noise)}, threshold ${fmt((ch.live ?? ch.threshold) * ch.noise)}`
      : ch.height_unit === "filtered" ? `threshold ${fmt(ch.threshold)} (noise ${fmt(ch.noise)})` : `probability ≥ ${fmt(ch.live ?? ch.threshold)}`;
    parts.push(" · ", el("span", { style: `color:${hex(OUTLINE[i % OUTLINE.length])}`, textContent: "● " }),
      el("b", { textContent: `channel ${ch.channel}: ${live ? "≈ " : ""}${n} foci` }), ` (${fmt(n / Math.max(1, p.cells))} per cell; ${thr})`);
  });
  if (p.remarks.length) parts.push(` · ${p.remarks.join("; ")}`);
  $("preview-info").replaceChildren(...parts);
}
// the slider moves: the markers and the count follow at once, the outlines when the server has answered
function liveThreshold(section, index, key, value) {
  if (section !== "foci" || !V.preview || !V.preview.channels[index]) return;
  const ch = V.preview.channels[index];
  if ((key === "sensitivity" && ch.height_unit === "× noise") || (key === "probability" && ch.height_unit === "probability")) {
    ch.live = value; if (viewer.markers[index]) viewer.markers[index].threshold = value;
    viewer.requestDraw(); previewInfo();
  }
}
viewer.onclick = (x, y) => {
  const p = V.preview, box = $("cell-info");
  if (!p || !p.table.cell) { box.hidden = true; return; }
  const k = p.binning / p.pixel_size; let best = -1, bestD = Infinity;
  p.table.cell.forEach((_, i) => {
    const d = Math.hypot(p.table.x[i] * k + p.border - x, p.table.y[i] * k + p.border - y);
    if (d < bestD && d < 1.2 * Math.sqrt(p.table.area[i] / Math.PI) * k) { best = i; bestD = d; }
  });
  if (best < 0) { box.hidden = true; return; }
  const rows = Object.entries(p.table).filter(([name, col]) => !["x", "y"].includes(name) && col[best] !== null && col[best] !== "")
    .map(([name, col]) => el("tr", {}, el("td", { textContent: name.replaceAll("_", " ") }), el("td", { textContent: fmt(col[best]) })));
  box.replaceChildren(el("table", {}, rows)); box.hidden = false;
};

function syncMode() {
  for (const b of $("mode").children) b.classList.toggle("on", b.dataset.mode === V.mode);
  const i = V.info;
  $("z-wrap").hidden = V.mode !== "slice" || !i || i.slices === 1;
  $("t-wrap").hidden = !i || i.frames === 1;
  $("z-label").textContent = i ? `${V.z + 1} / ${i.slices}` : ""; $("t-label").textContent = i ? `${V.t + 1} / ${i.frames}` : "";
}
function setMode(mode) { if (!V.info || V.info.slices === 1) return; V.mode = mode; syncMode(); syncLayers(); }
function setZ(z) { if (!V.info) return; V.z = Math.min(V.info.slices - 1, Math.max(0, z)); $("z").value = V.z; syncMode(); syncLayers(); }
function setT(t) {
  if (!V.info) return;
  t = Math.min(V.info.frames - 1, Math.max(0, t)); if (t === V.t) return;
  V.t = t; $("t").value = t; V.preview = null; V.cand = []; syncMode(); syncLayers(); if ($("live").checked) runPreview();
}
$("mode").addEventListener("click", (e) => { if (e.target.dataset.mode) setMode(e.target.dataset.mode); });
$("z").addEventListener("input", (e) => setZ(+e.target.value));
$("t").addEventListener("change", (e) => setT(+e.target.value));
$("fit").onclick = () => viewer.fit();
$("image-select").onchange = (e) => openImage(e.target.value);
for (const id of ["show-cells", "show-foci", "show-markers"]) $(id).onchange = syncLayers;
$("live").onchange = (e) => { if (e.target.checked) runPreview(); };

function renderChannels() {
  $("channels").replaceChildren(...V.ch.map((ch, c) => {
    const lo = el("input", { type: "range", min: 0, max: ch.top, value: ch.lo, title: "Darkest value shown", "aria-label": `Minimum of channel ${c + 1}`,
      oninput: (e) => { ch.lo = Math.min(+e.target.value, ch.hi - 1); syncLayers(); } });
    const hi = el("input", { type: "range", min: 1, max: ch.top, value: ch.hi, title: "Brightest value shown", "aria-label": `Maximum of channel ${c + 1}`,
      oninput: (e) => { ch.hi = Math.max(+e.target.value, ch.lo + 1); syncLayers(); } });
    return el("div", { class: "channel" },
      el("input", { type: "checkbox", id: `ch${c}`, checked: ch.on, title: `Show channel ${c + 1} (key ${c + 1})`, onchange: (e) => { ch.on = e.target.checked; syncLayers(); } }),
      el("input", { type: "color", value: hex(ch.color), title: "Colour of this channel", "aria-label": `Colour of channel ${c + 1}`,
        oninput: (e) => { ch.color = unhex(e.target.value); store("fa.colors", { ...(store("fa.colors") || {}), [c]: e.target.value }); syncLayers(); } }),
      el("label", { htmlFor: `ch${c}`, textContent: `Ch ${c + 1}` }),
      el("div", { class: "ranges" }, lo, hi),
      el("button", { textContent: "Auto", title: "Contrast from the image (0.1 to 99.9 percentile)",
        onclick: () => { [ch.lo, ch.hi] = ch.auto; lo.value = ch.lo; hi.value = ch.hi; syncLayers(); } }));
  }));
}

document.addEventListener("keydown", (e) => {
  if ($("page-analyze").hidden || e.ctrlKey || e.metaKey || e.altKey) return;
  const tag = e.target.tagName, type = e.target.type;
  if (tag === "SELECT" || tag === "TEXTAREA" || (tag === "INPUT" && !["checkbox", "range", "color"].includes(type))) return;
  const toggle = (id) => { $(id).checked = !$(id).checked; $(id).onchange({ target: $(id) }); };
  const k = e.key.toLowerCase();
  if (k >= "1" && k <= "9" && V.ch[+k - 1]) { V.ch[+k - 1].on = !V.ch[+k - 1].on; renderChannels(); syncLayers(); }
  else if (k === "c") toggle("show-cells"); else if (k === "f") toggle("show-foci"); else if (k === "m") toggle("show-markers");
  else if (k === "l") toggle("live"); else if (k === "p") setMode(V.mode === "mip" ? "slice" : "mip");
  else if (k === ",") setZ(V.z - 1); else if (k === ".") setZ(V.z + 1);
  else if (k === "[") setT(V.t - 1); else if (k === "]") setT(V.t + 1); else if (k === "0") viewer.fit();
  else return;
  e.preventDefault();
});

// ---------------------------------------------------------------- batch
let pollTimer = 0;
function showProgress(b) {
  S.batch = b;
  $("progress").hidden = !b.running; $("run").disabled = $("run-one").disabled = b.running;
  $("progress-bar").style.width = `${b.total ? 100 * b.done / b.total : 0}%`;
  $("progress-text").textContent = b.running ? `${b.done} of ${b.total} images${b.current ? ` · ${b.current}` : ""}`
    : b.errors && b.errors.length ? b.errors.join("\n") : "";
}
async function poll() {
  clearTimeout(pollTimer);
  try {
    const b = await (await api("/api/batch")).json(); showProgress(b);
    if (b.running) pollTimer = setTimeout(poll, 600);
    else if (b.total) { status(`Analyzed ${b.done - (b.errors || []).length} of ${b.total} images`, !!(b.errors || []).length); if (!(b.errors || []).length) showPage("results"); }
  } catch (e) { status(e.message, true); }
}
async function run(images) {
  try { showProgress(await (await api("/api/batch/start", { settings: S.settings, images })).json()); poll(); }
  catch (e) { status(e.message, true); }
}
$("run").onclick = () => run(null);
$("run-one").onclick = () => { if (V.id) run([V.id]); else status("Choose an image first", true); };

$("quit").onclick = async () => {
  if (S.batch.running && !window.confirm("An analysis is running. Stop Foci Analyzer anyway?")) return;
  try { await api("/api/quit", {}); } catch { /* it is stopping */ }
  document.body.replaceChildren(el("p", { style: "margin:60px auto;font-size:16px", textContent: "Foci Analyzer has stopped. You can close this tab." }));
};

// ---------------------------------------------------------------- start
window.FA = { viewer, V, S };  // for looking at the state from the browser console
initResults({
  openCell: async (id, frame, x, y, width) => {  // from a clicked point on the Results page
    showPage("analyze");
    if (V.id !== id) await openImage(id);
    if (frame) setT(frame - 1);
    if (x !== undefined) viewer.centerOn(x, y, width);
  },
});
(async () => {
  try {
    Object.assign(S, await (await api("/api/state")).json());
    if (!S.conditions.length) S.conditions.push({ name: "Condition 1", folders: [], pattern: "*", recursive: false });
    renderConditions(); renderImages(); renderSettings(); syncMode(); showProgress(S.batch);
    if (S.batch.running) poll();
  } catch (e) { status(e.message, true); }
})();
