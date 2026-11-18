// Foci Analyzer window: pages, the experiment, the settings form (drawn from the schema), the viewer with
// the live preview, and the batch run. The Results page is in results.js, the viewer in viewer.js.

import { $, el, api, status, fmt, helpIcon, store } from "./util.js";
import { Viewer, parsePolys, inside } from "./viewer.js";
import { initResults, loadResults } from "./results.js";

const S = { schema: [], settings: {}, conditions: [], images: [], remarks: [], output: "", models: {}, batch: {},
  catalog: { models: {}, docs: {} }, installed: null, analyzed: { ids: [], current: false } };

// ---------------------------------------------------------------- pages
function showPage(name) {
  for (const b of $("tabs").children) b.classList.toggle("on", b.dataset.page === name);
  for (const p of document.querySelectorAll(".page")) p.hidden = p.id !== `page-${name}`;
  if (name === "analyze") viewer.requestDraw();
  if (name === "results") loadResults();
  if (name === "settings") loadProgram();
}
$("tabs").addEventListener("click", (e) => { if (e.target.dataset.page) showPage(e.target.dataset.page); });

// ---------------------------------------------------------------- experiment
function renderConditions() {
  const box = $("conditions"); box.replaceChildren();
  S.conditions.forEach((c, i) => {
    const folders = el("div");
    c.folders.forEach((f, j) => folders.append(el("div", { class: "folder" },
      el("input", { type: "text", value: f, placeholder: "Folder or image file", "aria-label": "Folder", oninput: (e) => { c.folders[j] = e.target.value; } }),
      el("button", { class: "link", textContent: "Remove", onclick: () => { c.folders.splice(j, 1); renderConditions(); scan(); } }))));
    box.append(el("div", { class: "condition" },
      el("div", { class: "row" },
        el("input", { type: "text", class: "name", value: c.name, placeholder: "Condition name", "aria-label": "Condition name", oninput: (e) => { c.name = e.target.value; }, onchange: () => scan() }),
        el("label", {}, "Files named ", el("input", { type: "text", class: "pattern", value: c.pattern, title: "Pattern for the file name, * = anything", oninput: (e) => { c.pattern = e.target.value; }, onchange: () => scan() })),
        el("label", {}, el("input", { type: "checkbox", checked: c.recursive, onchange: (e) => { c.recursive = e.target.checked; scan(); } }), " with sub-folders"),
        el("button", { class: "link", textContent: "Remove condition", onclick: () => { S.conditions.splice(i, 1); renderConditions(); scan(); } })),
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
  if (paths === null) { c.folders.push(""); renderConditions(); return; }  // no dialog: an empty box to type or paste the folder in
  if (!paths.length) return;
  c.folders.push(...paths);
  renderConditions(); scan();
}
function renderImages() {
  const body = $("images").tBodies[0]; body.replaceChildren();
  for (const im of S.images) body.append(el("tr", { class: "click", onclick: () => { showPage("analyze"); openImage(im.id); } },
    el("td", { class: analyzedMark(im.id), textContent: analyzedMark(im.id) ? "✓" : "", title: analyzedTitle(im.id) }), el("td", { textContent: im.name }), el("td", { textContent: im.condition }), el("td", { class: "path", textContent: im.folder, title: im.folder })));
  $("image-count").textContent = S.images.length ? `(${S.images.length})` : "";
  $("remarks").replaceChildren(...S.remarks.map((r) => el("li", { textContent: r })));
  const sel = $("image-select"); sel.replaceChildren(el("option", { value: "", textContent: "Choose an image…" }));
  for (const im of S.images) {
    const m = analyzedMark(im.id);
    sel.append(el("option", { value: im.id, textContent: `${m ? "✓ " : ""}${im.condition} · ${im.name}`, title: analyzedTitle(im.id),
      style: { done: "color:#2e9e5b", viewed: "color:#2f7fc8", old: "color:#888", "": "" }[m] }));
  }
  if (V.id && !S.images.some((im) => im.id === V.id)) closeImage();  // it was removed from the experiment
  sel.value = V.id;
  $("output").value = S.output || "";
}
// green: in the results folder with the current settings (Analyze all images); blue: analyzed in the viewer
// with the current settings (Analyze this image); grey: analyzed, but with other settings
const viewed = {};  // image id: the settings it was analyzed with in the viewer (this session)
function analyzedMark(id) {
  if (S.analyzed.ids.includes(id) && S.analyzed.current) return "done";
  if (viewed[id] && viewed[id] === JSON.stringify(S.settings)) return "viewed";
  return S.analyzed.ids.includes(id) || viewed[id] ? "old" : "";
}
function analyzedTitle(id) {
  return { done: "In the results folder, analyzed with the current settings (Analyze all images)",
    viewed: "Analyzed in the viewer with the current settings (not in the results folder)",
    old: "Analyzed with other settings than the current ones", "": "" }[analyzedMark(id)];
}
async function refreshAnalyzed() {
  try { S.analyzed = await (await api("/api/analyzed")).json(); renderImages(); } catch { /* not important */ }
}
async function scan() {
  try {
    status("Looking for images…", false, -1);
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
  if (await scan() && S.images.length) status(`${S.images.length} images · click an image in the list below to open it on the Analyze page`);
}
$("quick-files").onclick = () => quick("files");
$("quick-folder").onclick = () => quick("folder");

// ---------------------------------------------------------------- settings form
let changeTimer = 0;
function settingsChanged(delay = 300) {
  clearTimeout(changeTimer);
  changeTimer = setTimeout(async () => {
    updateChannelColors(); renderOverlays(); syncMode();
    if (V.id && V.preview && $("live").checked) return runPreview();  // the preview request also saves the settings
    try {
      S.settings = await (await api("/api/settings", { settings: S.settings })).json();
      status(V.id ? "Settings saved · press Analyze this image (A) to see their effect" : "Settings saved");
      refreshAnalyzed();
    }
    catch (e) { status(e.message, true); }
  }, delay);
}
const installed = (net) => !S.installed || !!S.installed[net];
const OWN = "\u0001own";
function modelControl(f, group, id) {
  const net = group.method && S.catalog.models[group.method] ? group.method : "spotiflow", v = group[f.key];
  const list = S.catalog.models[net] || [], known = list.some((m) => m.name === v);
  const save = (value) => { group[f.key] = value; settingsChanged(); };
  const path = el("input", { type: "text", value: known ? "" : v, placeholder: "The folder or file of your model", "aria-label": "Path of your own model", onchange: (e) => save(e.target.value.trim()) });
  const browse = el("button", { textContent: "Browse…", onclick: async () => {
    const p = await pick(net === "cellpose" || net === "instanseg" ? "files" : "folder");
    if (p && p[0]) { path.value = p[0]; save(p[0]); } } });
  const own = el("div", { class: "row own-model" }, path, browse);
  own.hidden = known;
  const sel = el("select", { id, onchange: (e) => {
    if (e.target.value === OWN) { own.hidden = false; path.focus(); } else { own.hidden = true; save(e.target.value); renderSettings(); } } },
  list.map((m) => el("option", { value: m.name, disabled: !installed(m.network) && m.name !== v,
    textContent: m.label + (installed(m.network) ? (m.name === "cpsam_v2" && S.gpu && S.gpu.cellpose4 === false ? " · slow without a modern graphics card" : "")
      : m.network === "cellpose3" && installed("cellpose") ? " · needs Cellpose 3" : " · not installed (Settings page)") })),
  el("option", { value: OWN, textContent: "Own model…" }));
  sel.value = known ? v : OWN;
  return el("div", { class: "field" }, el("label", { htmlFor: id, textContent: f.label }), helpIcon(f.help, S.catalog.docs[net]),
    el("div", { class: "control" }, sel), own);
}
function control(f, group, sectionKey, index) {
  const set = (v, delay) => {
    group[f.key] = v;
    if (sectionKey === "segmentation" && f.key === "dimension" && v === "3D") S.settings.foci_general.dimension = "3D";
    settingsChanged(delay); if (f.type === "choice" || f.type === "bool") renderSettings();
  };
  const id = `f${Math.random().toString(36).slice(2)}`, v = group[f.key];
  if (f.catalog) return modelControl(f, group, id);
  // cells in 3D: the foci are found in 3D too (a projected focus has no slice to belong to)
  const lock3d = sectionKey === "foci_general" && f.key === "dimension" && S.settings.segmentation.dimension === "3D";
  const number = (e) => { const n = f.type === "int" ? parseInt(e.target.value, 10) : parseFloat(e.target.value); return Number.isFinite(n) ? n : null; };
  let input, range = null;
  // shown, but not used with the current choices (for example the border shift when foci count in the whole image)
  const off = !!f.off_when && Object.entries(f.off_when).every(([k, values]) => values.includes(group[k]));
  if (f.type === "bool") input = el("input", { id, type: "checkbox", checked: v, onchange: (e) => set(e.target.checked) });
  else if (f.type === "choice") {
    input = el("select", { id, onchange: (e) => set(e.target.value) },
      f.choices.map((c) => el("option", { value: c, textContent: c, disabled: lock3d && c !== "3D" })));
    input.value = v;
  }
  else if (f.type === "str") input = el("input", { id, type: "text", value: v, onchange: (e) => set(e.target.value) });
  else {
    input = el("input", { id, type: "number", value: v, min: f.min ?? "", max: f.max ?? "", step: f.step ?? (f.type === "int" ? 1 : "any"),
      onchange: (e) => { const n = number(e); if (n === null) e.target.value = group[f.key]; else { set(n); if (range) range.value = n; } } });
    if (f.slider) range = el("input", { type: "range", min: f.min, max: f.max, step: f.step, value: v, "aria-label": f.label,
      oninput: (e) => { const n = number(e); input.value = n; liveThreshold(sectionKey, index, f.key, n); set(n, 120); } });
  }
  if (off) { input.disabled = true; if (range) range.disabled = true; }
  return el("div", { class: `field${f.type === "bool" ? " check" : ""}${range ? " slider" : ""}${off ? " off" : ""}` },
    el("label", { htmlFor: id, textContent: f.label }), helpIcon(f.help),
    el("div", { class: "control" }, input, f.unit ? el("span", { class: "unit", textContent: f.unit }) : null), range,
    lock3d ? el("div", { class: "sub", textContent: "The cells are segmented in 3D, so the foci are detected in 3D too. Set 'Segment in' to 2D to use a projection." }) : null);
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
    list.forEach((g, i) => div.append(el("div", sec.key === "foci" ? { class: "group tinted", style: `--g:${hex(OUTLINE[i % OUTLINE.length])}` } : { class: "group" },
      el("div", { class: "group-head" }, `${sec.title} ${i + 1}`,
        list.length > 1 ? el("button", { class: "link", textContent: "Remove", onclick: () => { list.splice(i, 1); renderSettings(); settingsChanged(); } }) : null),
      ...fieldsOf(sec, g, i).filter(Boolean))));
    box.append(div);
  }
  renderOverlays();
}

// ---------------------------------------------------------------- viewer
const viewer = new Viewer($("canvas"), $("markers"), (url) => api(url));
const V = { id: "", info: null, mode: "mip", z: 0, t: 0, ch: [], preview: null, cand: [], load: 0, slots: {}, fociOn: {}, lastDim: "",
  polys: {}, polysPrev: {}, sel: 0, selPt: null, multi: new Map(), vr: [0, 1] };
const rgba = (c, a) => `rgba(${c.map((v) => Math.round(v * 255)).join(",")},${a})`;
// the projection the foci are detected on, when it is not the max projection: [name in a request, name on the button]
const DET = { "extended depth of focus": ["edf", "EDF"], "sum projection": ["sum", "Sum"] };
const NUCLEI = [0, 0.5, 1], CELLS = [1, 0.5, 0], OTHER = [[0, 1, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]];
const OUTLINE = [[1, 0.88, 0.3], [0.3, 0.88, 1], [1, 0.6, 0.2], [0.8, 0.6, 1]];  // foci channel 1, 2, …
const hex = (c) => `#${c.map((v) => Math.round(v * 255).toString(16).padStart(2, "0")).join("")}`;
const unhex = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16) / 255);

function channelColor(c) {
  const saved = (store("fa.colors") || {})[c];
  if (saved) return unhex(saved);
  if (c + 1 === S.settings.image.nuclei_channel) return NUCLEI;
  if (c + 1 === S.settings.image.cell_channel) return CELLS;
  // the foci channels first, then the rest, each the next colour: no two channels alike while colours last
  const n = V.info ? V.info.channels : c + 1, special = [S.settings.image.nuclei_channel - 1, S.settings.image.cell_channel - 1];
  const foci = [...new Set(S.settings.foci.map((f) => f.channel - 1))].filter((x) => !special.includes(x));
  const rest = [...Array(n).keys()].filter((x) => !special.includes(x) && !foci.includes(x));
  return OTHER[Math.max(0, [...foci, ...rest].indexOf(c)) % OTHER.length];
}
function updateChannelColors() { V.ch.forEach((ch, c) => { ch.color = channelColor(c); }); renderChannels(); syncLayers(); }

async function openImage(id, analyze = false) {
  const load = ++V.load;
  $("image-select").value = id;
  if (!id) { closeImage(); return; }
  try {
    status("Reading the image…", false, -1);
    const info = await (await api(`/api/image/info?path=${encodeURIComponent(id)}`)).json();
    if (load !== V.load) return;
    play(false);
    Object.assign(V, { id, info, z: Math.floor(info.slices / 2), t: 0, preview: null, cand: [], slots: {}, lastDim: "", polys: {}, polysPrev: {}, sel: 0, selPt: null });
    V.multi.clear();
    const top = info.dtype === "uint8" ? 255 : 65535;
    V.ch = info.ranges.map(([lo, hi], c) => ({ on: true, lo, hi, auto: [lo, hi], top: Math.min(top, Math.max(hi * 3, lo + 10)), color: channelColor(c) }));
    const kept = store("fa.display");  // the contrast and on/off of the image shown before
    if (kept && kept.n === info.channels && kept.dtype === info.dtype && kept.scale === info.scale) {
      kept.ch.forEach((k, c) => { Object.assign(V.ch[c], { on: k.on, lo: k.lo, hi: k.hi }); V.ch[c].top = Math.min(top, Math.max(V.ch[c].top, k.hi + 1)); });
    }
    if (info.slices === 1) V.mode = "mip";
    $("z").max = info.slices - 1; $("z").value = V.z; $("t").max = info.frames - 1; $("t").value = 0;
    $("viewer-empty").hidden = true; showCell(null);
    const large = info.width * info.height >= LARGE;
    $("part-wrap").hidden = !large; $("part").checked = large;
    viewer.setImage(info.width, info.height); viewer.pixelSize = info.pixel_size; viewer.unit = info.unit;
    renderCellColor(); renderChannels(); syncMode(); syncLayers(); viewer.fit(); status("");
    previewInfo();
    if (analyze) runPreview();  // otherwise the user first checks the settings and presses Analyze this image,
    else await loadStored();  // unless the last batch run already analyzed it with these settings
  } catch (e) { status(e.message, true); }
}
function closeImage() {
  V.load++; Object.assign(V, { id: "", info: null, preview: null, cand: [], ch: [], slots: {}, polys: {}, polysPrev: {}, sel: 0, selPt: null });
  viewer.setImage(0, 0); $("viewer-empty").hidden = false; $("image-select").value = "";
  renderCellColor(); renderChannels(); syncMode(); syncLayers(); previewInfo(); showCell(null);
}
// The tiles shown before in this place stay until the new ones are there, so nothing flickers.
function slot(name, key) {
  const s = V.slots[name] || (V.slots[name] = {});
  if (s.key !== key) { s.prev = s.key; s.key = key; }
  return s.prev;
}
function saveDisplay() {
  if (V.info) store("fa.display", { n: V.info.channels, dtype: V.info.dtype, scale: V.info.scale, ch: V.ch.map((c) => ({ on: c.on, lo: c.lo, hi: c.hi })) });
}
// the result of the last batch run for this image and frame, when it used the current settings
async function loadStored() {
  const id = V.id, t = V.t;
  if (!id) return;
  try {
    const r = await (await api("/api/preview", { path: id, t, settings: S.settings, stored: true })).json();
    if (!r.stored || id !== V.id || t !== V.t || V.preview) return;
    const cand = await Promise.all(r.channels.map(async (ch) =>
      new Float32Array(await (await api(`/api/candidates?path=${encodeURIComponent(id)}&t=${t}&ch=${ch.channel}`)).arrayBuffer())));
    if (id !== V.id || t !== V.t || V.preview) return;
    V.preview = r; V.cand = cand;
    renderCellColor(); syncMode(); syncLayers(); previewInfo();
    status(`From the last run of Analyze all images (same settings) · ${r.cells} cells`);
  } catch { /* nothing stored: the image waits for Analyze this image */ }
}
function syncLayers() {
  saveDisplay();
  if (!V.info) { viewer.layers = []; viewer.markers = []; viewer.requestDraw(); return; }
  const det = DET[S.settings.foci_general.dimension];
  const z = V.mode === "slice" ? V.z : V.mode === "det" && det ? det[0] : "mip", base = `path=${encodeURIComponent(V.id)}&t=${V.t}&z=${z}`;
  const region = (extra) => (level, x, y, w, h) => `${extra}&${base}&level=${level}&x=${x}&y=${y}&w=${w}&h=${h}`;
  const layers = V.ch.map((ch, c) => {
    const key = `i|${V.id}|${V.t}|${c}|${z}`;
    return { type: "image", key, prevKey: slot(`i${c}`, key), url: region(`/api/image/region?c=${c}`), lo: ch.lo, hi: ch.hi, color: ch.color, on: ch.on };
  });
  viewer.markers = [];
  if (V.preview) {
    const k = V.preview.key, by = $("cell-color").value, smooth = $("smooth").checked; let key;
    if (by === "label") {  // every cell its own colour: where one cell ends and the next begins
      key = `l|cells|${z}|${k}`;
      layers.push({ type: "value", key, prevKey: slot("value", key), url: region("/api/overlay?kind=cells"), alpha: 0.5, cat: true, on: true });
    } else if (by) {  // under the outlines: every cell filled with the colour of one of its numbers
      key = `v|${by}|${z}|${k}`;
      layers.push({ type: "value", key, prevKey: slot("value", key), url: region(`/api/overlay?kind=value&col=${encodeURIComponent(by)}`), alpha: 0.6, range: V.vr, on: true });
    }
    key = `l|cells|${z}|${k}`;
    layers.push({ type: "labels", key, prevKey: slot("cells", key), url: region("/api/overlay?kind=cells"),
      color: [1, 1, 1, 0.85], fill: 0, on: $("show-cells").checked && !smooth });
    V.preview.channels.forEach((ch, i) => {
      const col = OUTLINE[i % OUTLINE.length], shown = V.fociOn[ch.channel] !== false;
      key = `l|foci${ch.channel}|${z}|${k}`;
      layers.push({ type: "labels", key, prevKey: slot(`foci${ch.channel}`, key), url: region(`/api/overlay?kind=foci&ch=${ch.channel}`),
        color: [...col, 1], fill: 0.35, on: $("show-foci").checked && shown && !smooth });
      viewer.markers.push({ pts: V.cand[i], color: hex(col), threshold: ch.live ?? ch.threshold, on: $("show-markers").checked && shown });
    });
    if ((V.preview.coloc || []).length) {  // the colocalization map, on top
      const on = $("show-coloc").checked;
      const line = on && !smooth;
      key = `l|coloc|${z}|${k}`;
      const prev = slot("coloc", key), url = region("/api/overlay?kind=coloc");
      layers.push({ type: "labels", key, prevKey: prev, url, color: [0, 0, 0, 1], width: 2.6, on: line });  // the black edge,
      layers.push({ type: "labels", key, prevKey: slot("coloc", key), url, color: [1, 1, 1, 1], width: 1, on: line });  // the white line on it
      key = `l|overlap|${z}|${k}`;
      layers.push({ type: "labels", key, prevKey: slot("overlap", key), url: region("/api/overlay?kind=overlap"), color: [0, 0, 0, 1], fill: 0.75, hatch: 3, on });
    }
  }
  viewer.box = V.preview && V.preview.part ? V.preview.part : null;
  viewer.layers = layers; viewer.zslice = V.mode === "slice" ? V.z : null;
  syncShapes();
}
// outlines as polygons, asked once per plane; until they are there the previous ones stay
function polys(kind, extra = "") {
  const z = shownZ(), name = kind + extra, k = `${name}|${z}|${V.preview.key}|${V.id}|${V.t}`;
  const have = V.polys[k];
  if (have === undefined) {
    V.polys[k] = null;
    const id = V.id;
    api(`/api/contours?kind=${kind}${extra}&path=${encodeURIComponent(V.id)}&t=${V.t}&z=${z}`).then((r) => r.arrayBuffer())
      .then((b) => { if (id === V.id) { V.polys[k] = parsePolys(b); syncShapes(); } })
      .catch((e) => {
        V.polys[k] = false;
        if (kind !== "cells" && $("smooth").checked) { $("smooth").checked = false; status(`${e.message}: outlines along the pixels`, true); syncLayers(); }
      });
  }
  if (have) V.polysPrev[name] = have;
  return have || V.polysPrev[name] || null;
}
function syncShapes() {
  const p = V.preview;
  if (!p || !V.info) { viewer.shapes = []; viewer.numbers = null; viewer.requestDraw(); return; }
  const smooth = $("smooth").checked, cellsOn = $("show-cells").checked, removed = new Set(p.removed || []);
  const cells = polys("cells"), sh = [];
  if (p.outside !== 0 && cellsOn) sh.push({ polys: polys("reach"), stroke: "rgba(215,215,215,.8)", width: 1, dash: [4, 3], on: true });
  if (smooth) {
    sh.push({ polys: cells, stroke: "rgba(255,255,255,.85)", width: 1.3, on: cellsOn });
    p.channels.forEach((ch, i) => {
      const col = OUTLINE[i % OUTLINE.length];
      sh.push({ polys: polys("foci", `&ch=${ch.channel}`), stroke: hex(col), width: 1.2, fill: rgba(col, 0.35),
        on: $("show-foci").checked && V.fociOn[ch.channel] !== false });
    });
    if ((p.coloc || []).length) sh.push({ polys: polys("coloc"), stroke: "#fff", width: 1, halo: 1.6, on: $("show-coloc").checked });
  }
  if (removed.size) sh.push({ polys: cells, only: removed, stroke: "#ff3b3b", width: 2, halo: 1, on: true });
  if (V.sel) sh.push({ polys: cells, only: new Set([V.sel]), stroke: "#ffd400", width: 2.6, halo: 1.6, on: true });
  if (V.multi.size) sh.push({ polys: cells, only: new Set(V.multi.keys()), stroke: "#ffd400", width: 2.6, halo: 1.6, on: true });
  viewer.shapes = sh;
  viewer.numbers = { on: $("show-numbers").checked, polys: cells, removed };
  viewer.requestDraw();
}
viewer.onview = () => {
  const i = V.info; if (!i) return;
  const cal = i.pixel_size ? ` · ${fmt(i.pixel_size)} ${i.unit}/px` : " · not calibrated";
  $("canvas").title = ""; $("image-select").title = `${i.width} × ${i.height} px · ${i.channels} channels · ${i.slices} slices${cal} · zoom ${Math.round(viewer.view.s * 100)} %`;
};

// For a large image the preview can be limited to what is visible (the part is fixed when the analysis starts).
const LARGE = 6e6;  // pixels: from here on the choice "Visible part only" is offered
function visiblePart() {
  if (!V.info || $("part-wrap").hidden || !$("part").checked) return null;
  const { s, tx, ty } = viewer.view, [w, h] = viewer.size();
  let x0 = -tx / s, y0 = -ty / s, x1 = (w - tx) / s, y1 = (h - ty) / s;
  const grow = (a, b) => { const d = Math.max(0, 1024 - (b - a)) / 2; return [a - d, b + d]; };  // zoomed in far: still enough image for the noise
  [x0, x1] = grow(x0, x1); [y0, y1] = grow(y0, y1);
  return [Math.floor(x0), Math.floor(y0), Math.ceil(x1), Math.ceil(y1)];
}

// ---- preview: the analysis of the shown image with the current settings
let previewBusy = false, previewAgain = false;
async function runPreview() {
  if (!V.id) { status("Choose an image first", true); return; }
  if (previewBusy) { previewAgain = true; return; }
  previewBusy = true; status("Analyzing the shown image…", false, -1, true);
  const id = V.id, t = V.t, start = performance.now();
  let best = 0;  // stages of several channels run side by side: the bar only moves forward
  // what the server is doing, asked a few times per second while it works
  const watch = setInterval(async () => {
    try {
      const a = await (await api("/api/activity")).json();
      if (!previewBusy || !a.busy) return;
      const secs = a.seconds >= 2 ? ` · ${Math.round(a.seconds)} s` : "";
      if (a.steps) best = Math.max(best, (a.step - 0.5) / a.steps);
      status(`${a.text}…${secs}`, false, a.steps || best ? best : -1, true);  // once it shows how far, it never falls back
    } catch { /* the answer of the analysis itself reports problems */ }
  }, 300);
  try {
    const r = await (await api("/api/preview", { path: id, t, settings: S.settings, box: visiblePart() })).json();
    const cand = await Promise.all(r.channels.map(async (ch) =>
      new Float32Array(await (await api(`/api/candidates?path=${encodeURIComponent(id)}&t=${t}&ch=${ch.channel}`)).arrayBuffer())));
    viewed[id] = JSON.stringify(S.settings);
    if (id === V.id && t === V.t) {
      V.preview = r; V.cand = cand;
      // foci found on a projection: show that projection (once, so a later choice of the user stays)
      if (V.info.slices > 1 && r.dimension !== "3D" && r.dimension !== V.lastDim) V.mode = DET[r.dimension] ? "det" : "mip";
      V.lastDim = r.dimension;
      renderCellColor(); syncMode(); syncLayers(); previewInfo(); selectCell(V.sel, V.selPt);
      const secs = (performance.now() - start) / 1000;
      status(`Ready · ${r.cells} ${r.cells === 1 ? "cell" : "cells"}, ${r.channels.map((c) => `${c.foci} foci in channel ${c.channel}`).join(", ")}`
        + ` · analyzed in ${secs < 10 ? secs.toFixed(1) : Math.round(secs)} s`);
    } else status("");
  } catch (e) {
    if (e.message === "stopped") { previewAgain = false; status("Stopped"); } else status(e.message, true);
  }
  clearInterval(watch);
  previewBusy = false;
  refreshAnalyzed();
  if (previewAgain) { previewAgain = false; runPreview(); }
}
function liveCount(pts, thr, removed) {  // foci above the threshold, not in a removed cell
  let n = countAbove(pts, thr);
  if (removed.size) for (let i = n - 1; i >= 0; i--) if (removed.has(pts[i * 6 + 4])) n--;
  return n;
}
function countAbove(pts, thr) {  // heights are sorted from high to low
  let lo = 0, hi = pts.length / 6;
  while (lo < hi) { const m = (lo + hi) >> 1; if (pts[m * 6 + 3] >= thr) lo = m + 1; else hi = m; }
  return lo;
}
function previewInfo() {
  const p = V.preview, box = $("preview-info");
  if (!p) { box.replaceChildren(V.id ? "Not analyzed yet. Check the settings, then press Analyze this image (A)." : ""); return; }
  const lines = [p.stored ? el("div", { class: "sub", style: "margin-left:0", textContent: "From the last run of Analyze all images:" }) : null,
    p.part ? el("div", { class: "sub", style: "margin-left:0", textContent: "Of the part inside the dashed frame only:" }) : null, el("div", {}, el("b", { textContent: `${p.cells} ${p.cells === 1 ? "cell" : "cells"}` }))];
  p.channels.forEach((ch, i) => {
    const live = ch.live !== undefined, n = live ? liveCount(V.cand[i], ch.live, new Set(p.removed || [])) : ch.foci;
    const thr = ch.height_unit === "× noise" ? `noise ${fmt(ch.noise)}, threshold ${fmt((ch.live ?? ch.threshold) * ch.noise)}`
      : ch.height_unit === "filtered" ? `threshold ${fmt(ch.threshold)} (noise ${fmt(ch.noise)})` : `probability ≥ ${fmt(ch.live ?? ch.threshold)}`;
    lines.push(el("div", {}, el("span", { class: "dot", style: `background:${hex(OUTLINE[i % OUTLINE.length])}` }),
      el("b", { textContent: `Channel ${ch.channel}: ${live ? "≈ " : ""}${n} foci` }),
      el("div", { class: "sub", textContent: `${fmt(n / Math.max(1, p.cells))} per cell · ${thr}` })));
  });
  for (const c of p.coloc || []) lines.push(el("div", {}, `Channel ${c.a} foci on a channel ${c.b} focus: `,
    el("b", { textContent: `${c.with} of ${c.n} (${Math.round(100 * c.with / c.n)} %)` }),
    c.chance !== null ? el("div", { class: "sub", textContent: `by chance about ${Math.round(100 * c.chance)} %` }) : null));
  if ((p.coloc || []).length) lines.push(el("div", { class: "sub", textContent: "Map (key O): white outline with a black edge = colocalizing focus; hatched = where foci of two channels overlap." }));
  else if (p.channels.length > 1 && !S.settings.colocalization.enabled) lines.push(el("div", { class: "sub", textContent: "Colocalization is switched off (see the settings on the left)." }));
  if ((p.removed || []).length) lines.push(el("div", { class: "sub", style: "margin-left:0", textContent: `${p.removed.length} ${p.removed.length === 1 ? "cell" : "cells"} removed by hand (red outline), not counted` }));
  for (const r of p.remarks) lines.push(el("div", { textContent: r }));
  box.replaceChildren(...lines.filter(Boolean));
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
// the choice "Colour cells by" and its legend, from the numbers of the analyzed image
// the numbers of a cell in groups: the cell itself, the foci per channel, colocalization, the rest
function colourGroups(names) {
  const groups = new Map(), add = (g, n, text) => { if (!groups.has(g)) groups.set(g, []); groups.get(g).push([n, text]); };
  const FOCI = { "": "number of foci", total_intensity: "total intensity", total_size: "total size", mean_size: "mean size", mean_intensity: "mean intensity" };
  for (const n of names) {
    let m;
    if (n === "area") add("Cell", n, "area");
    else if ((m = n.match(/^mean_ch(\d+)$/))) add("Cell", n, `mean intensity, channel ${m[1]}`);
    else if ((m = n.match(/^foci_ch(\d+)_?(.*)$/))) add(`Foci in channel ${m[1]}`, n, FOCI[m[2]] || m[2].replaceAll("_", " "));
    else if (/^(coloc|fraction|manders|pearson)_/.test(n)) add("Colocalization", n, n.replace(/_ch(\d+)/g, " ch$1").replaceAll("_", " "));
    else add("Other", n, n.replaceAll("_", " "));
  }
  return [...groups];
}
function renderCellColor() {
  const p = V.preview, sel = $("cell-color"), keep = sel.value;
  const names = p && p.table.cell ? Object.keys(p.table).filter((n) => !["cell", "x", "y"].includes(n) && p.table[n].some((v) => typeof v === "number")) : [];
  $("cell-color-wrap").hidden = !names.length; $("coloc-wrap").hidden = !(p && (p.coloc || []).length);
  sel.replaceChildren(el("option", { value: "", textContent: "nothing" }), el("option", { value: "label", textContent: "label (a colour per cell)" }),
    ...colourGroups(names).map(([group, items]) => el("optgroup", { label: group }, items.map(([n, text]) => el("option", { value: n, textContent: text })))));
  sel.value = names.includes(keep) || (keep === "label" && names.length) ? keep : "";
  cellLegend();
}
function cellLegend() {
  const by = $("cell-color").value, v = by && by !== "label" && V.preview && V.preview.table[by] ? V.preview.table[by].filter((x) => typeof x === "number") : [];
  $("cell-legend").hidden = !v.length;
  if (!v.length) return;
  const lo = Math.min(...v), hi = Math.max(...v), at = (f) => lo + f * (hi - lo);
  $("legend-lo").textContent = fmt(at(V.vr[0])); $("legend-hi").textContent = fmt(at(V.vr[1]));
  $("legend-a").value = Math.round(1000 * V.vr[0]); $("legend-b").value = Math.round(1000 * V.vr[1]);
  $("legend-fill").style.left = `${100 * V.vr[0]}%`; $("legend-fill").style.right = `${100 - 100 * V.vr[1]}%`;
}
$("cell-color").onchange = () => { V.vr = [0, 1]; cellLegend(); syncLayers(); };
$("legend-a").oninput = (e) => { V.vr = [Math.min(+e.target.value / 1000, V.vr[1] - 0.01), V.vr[1]]; cellLegend(); syncLayers(); };
$("legend-b").oninput = (e) => { V.vr = [V.vr[0], Math.max(+e.target.value / 1000, V.vr[0] + 0.01)]; cellLegend(); syncLayers(); };

function showCell(rows) {
  $("cell-info").replaceChildren(rows ? el("table", {}, rows) : V.preview ? "Click a cell in the image for its numbers." : "");
}
viewer.onclick = (x, y, e) => {
  const p = V.preview, ctrl = !!e && (e.ctrlKey || e.metaKey);
  if (!ctrl) V.multi.clear();
  if (!p || !p.table.cell) { selectCell(0, null); return; }
  const cells = polys("cells");  // the cell under the click; else the nearest centre
  let id = 0;
  if (cells) id = (cells.find((poly) => inside(poly, x, y)) || {}).id || 0;
  else {
    const k = p.binning / p.pixel_size; let bestD = Infinity;
    p.table.cell.forEach((c, i) => {
      const d = Math.hypot(p.table.x[i] * k + p.border + p.offset[0] - x, p.table.y[i] * k + p.border + p.offset[1] - y);
      if (d < bestD && d < 1.2 * Math.sqrt(p.table.area[i] / Math.PI) * k) { id = c; bestD = d; }
    });
  }
  if (ctrl) { if (id) pickCells([[id, [x, y]]], true); return; }
  selectCell(id, id ? [x, y] : null);
};
// several cells at once: Ctrl + click adds or drops one, Shift + drag a box adds those whose centre is in it
function cellCentres() {
  const p = V.preview, k = p.binning / p.pixel_size;
  return p.table.cell.map((c, i) => [c, [p.table.x[i] * k + p.border + p.offset[0], p.table.y[i] * k + p.border + p.offset[1]]]);
}
function pickCells(list, toggle) {
  if (!V.multi.size && V.sel && V.selPt) V.multi.set(V.sel, V.selPt);  // the cell selected before counts too
  for (const [id, pt] of list) { if (toggle && V.multi.has(id)) V.multi.delete(id); else V.multi.set(id, pt); }
  V.sel = 0; V.selPt = null;
  showSelection();
}
viewer.onbox = (x0, y0, x1, y1) => {
  if (!V.preview || !V.preview.table.cell) return;
  const cells = polys("cells");
  const hit = cellCentres().filter(([, [x, y]]) => x >= x0 && x <= x1 && y >= y0 && y <= y1)
    // a point that lies in the cell itself (the centre of a bent cell can lie outside it)
    .map(([id, pt]) => { const poly = cells && cells.find((q) => q.id === id); return [id, poly && !inside(poly, ...pt) ? pointIn(poly, pt) : pt]; });
  if (hit.length) pickCells(hit, false); else status("No cell has its centre in the box", true);
};
function pointIn(poly, [cx, cy]) {  // a point just inside the outline, on the way from a corner to the centre
  for (let i = 0; i < poly.pts.length; i += 2) {
    const x = poly.pts[i], y = poly.pts[i + 1], d = Math.hypot(cx - x, cy - y) || 1;
    for (const step of [0.5, 1, 2]) { const px = x + step * (cx - x) / d, py = y + step * (cy - y) / d; if (inside(poly, px, py)) return [px, py]; }
  }
  return [cx, cy];
}
function showSelection() {
  const n = V.multi.size, removed = new Set((V.preview && V.preview.removed) || []);
  if (!n) { selectCell(0, null); return; }
  const back = [...V.multi.keys()].every((id) => removed.has(id));
  $("cell-actions").hidden = false;
  $("toggle-cell").firstChild.textContent = back ? `Put ${n} ${n === 1 ? "cell" : "cells"} back ` : `Remove ${n} ${n === 1 ? "cell" : "cells"} from the analysis `;
  $("cell-info").replaceChildren(el("div", { textContent: `${n} ${n === 1 ? "cell" : "cells"} selected: ${[...V.multi.keys()].sort((a, b) => a - b).join(", ")}` }),
    el("div", { class: "muted", textContent: "Ctrl + click: add or drop a cell. Shift + drag: a box. Esc: clear." }));
  syncShapes();
}
function selectCell(id, pt) {
  const p = V.preview;
  V.sel = p ? id : 0; V.selPt = V.sel ? pt : null;
  const removed = !!p && (p.removed || []).includes(V.sel), row = p && p.table.cell ? p.table.cell.indexOf(V.sel) : -1;
  $("cell-actions").hidden = !V.sel;
  $("toggle-cell").firstChild.textContent = removed ? "Put back in the analysis " : "Remove from the analysis ";
  if (!V.sel) showCell(null);
  else if (removed) $("cell-info").replaceChildren(el("div", { textContent: `Cell ${V.sel} is removed from the analysis.` }));
  else if (row >= 0) showCell(Object.entries(p.table).filter(([name, col]) => !["x", "y"].includes(name) && col[row] !== null && col[row] !== "")
    .map(([name, col]) => el("tr", {}, el("td", { textContent: name.replaceAll("_", " ") }), el("td", { textContent: fmt(col[row]) }))));
  else $("cell-info").replaceChildren(el("div", { textContent: `Cell ${V.sel} is not in the tables (for example on the image edge).` }));
  syncShapes();
}
async function toggleCell() {
  if (V.multi.size && V.preview) {
    const removed = new Set(V.preview.removed || []), ids = [...V.multi.keys()], back = ids.every((id) => removed.has(id));
    const z = V.mode === "slice" ? V.z : -1;
    try {
      const r = await (await api("/api/cells/set", { path: V.id, t: V.t, remove: !back, points: [...V.multi.values()].map(([x, y]) => [x, y, z]) })).json();
      V.preview = r; V.multi.clear(); renderCellColor(); syncLayers(); previewInfo(); selectCell(0, null);
      status(back ? `${r.changed} cells are back in the analysis`
        : `${r.changed} cells removed from the analysis; remembered in the results folder (Analyze all images again to update its tables)`);
    } catch (e) { status(e.message, true); }
    return;
  }
  if (!V.sel || !V.selPt || !V.preview) return;
  const cell = V.sel;
  try {
    const r = await (await api("/api/cells/toggle", { path: V.id, t: V.t, x: V.selPt[0], y: V.selPt[1], z: V.mode === "slice" ? V.z : -1 })).json();
    V.preview = r; renderCellColor(); syncLayers(); previewInfo(); selectCell(cell, V.selPt);
    status(r.toggled === "removed" ? `Cell ${cell} removed from the analysis; remembered in the results folder (Analyze all images again to update its tables)`
      : `Cell ${cell} is back in the analysis`);
  } catch (e) { status(e.message, true); }
}
$("toggle-cell").onclick = toggleCell;

function syncMode() {
  const i = V.info, stack = !!i && i.slices > 1, det = S.settings.foci_general && DET[S.settings.foci_general.dimension];
  if (V.mode === "det" && !det) V.mode = "mip";
  const b = $("mode").querySelector("[data-mode=det]"); b.hidden = !det; b.textContent = det ? det[1] : "";
  for (const m of $("mode").children) m.classList.toggle("on", m.dataset.mode === V.mode);
  // the Z slider keeps its place in a projection (only dimmed), so the bar and the image do not jump
  $("mode-wrap").hidden = !stack; $("z-wrap").hidden = !stack;
  $("z-wrap").classList.toggle("off", V.mode !== "slice"); $("z").disabled = V.mode !== "slice";
  $("t-wrap").hidden = !i || i.frames === 1;
  $("z-label").textContent = i ? `${V.z + 1} / ${i.slices}` : ""; $("t-label").textContent = i ? `${V.t + 1} / ${i.frames}` : "";
}
function setMode(mode) { if (!V.info || V.info.slices === 1) return; V.mode = mode; syncMode(); syncLayers(); }
function nextMode() {
  const order = ["slice", "mip", ...(DET[S.settings.foci_general.dimension] ? ["det"] : [])];
  setMode(order[(order.indexOf(V.mode) + 1) % order.length]);
}
function stepZ(d) {  // the slice keys: they also leave a projection, otherwise nothing would be seen
  if (!V.info || V.info.slices === 1) return;
  if (V.mode !== "slice") { V.mode = "slice"; setZ(V.z); } else setZ(V.z + d);
}
function setZ(z) { if (!V.info) return; V.z = Math.min(V.info.slices - 1, Math.max(0, z)); $("z").value = V.z; syncMode(); syncLayers(); }
function setT(t) {
  if (!V.info) return;
  t = Math.min(V.info.frames - 1, Math.max(0, t)); if (t === V.t) return;
  const had = !!V.preview;
  V.t = t; $("t").value = t; V.preview = null; V.cand = []; V.sel = 0; V.multi.clear(); renderCellColor(); syncMode(); syncLayers(); previewInfo(); selectCell(0, null);
  if (had && $("live").checked) runPreview(); else loadStored();
}
$("mode").addEventListener("click", (e) => { if (e.target.dataset.mode) setMode(e.target.dataset.mode); });
$("z").addEventListener("input", (e) => setZ(+e.target.value));
$("t").addEventListener("change", (e) => setT(+e.target.value));
viewer.onstep = (axis, d) => axis === "z" ? stepZ(d) : setT(V.t + d);
// play through the time frames, 4 per second, round and round
let playing = null;
function play(on) {
  clearTimeout(playing); playing = null;
  $("t-play").innerHTML = on ? "&#10074;&#10074;" : "&#9654;"; $("t-play").setAttribute("aria-label", on ? "Pause" : "Play");
  if (!on || !V.info || V.info.frames < 2) return;
  const tick = () => { if (!playing) return; setT((V.t + 1) % V.info.frames); playing = setTimeout(tick, 250); };
  playing = setTimeout(tick, 250);
}
$("t-play").onclick = () => play(!playing);
$("fit").onclick = () => viewer.fit();
$("image-select").onchange = (e) => openImage(e.target.value);
$("analyze-now").onclick = () => runPreview();
for (const id of ["show-cells", "show-foci", "show-markers", "show-coloc"]) $(id).onchange = syncLayers;
for (const id of ["show-numbers", "smooth"]) {
  $(id).checked = store(`fa.${id}`) === true;
  $(id).onchange = (e) => { store(`fa.${id}`, e.target.checked); syncLayers(); };
}
$("live").checked = store("fa.live") !== false;
$("live").onchange = (e) => store("fa.live", e.target.checked);

function renderChannels() {
  $("channels").replaceChildren(...V.ch.map((ch, c) => {
    const values = el("span", { class: "values" }), fill = el("div", { class: "fill" });
    const show = () => {
      values.textContent = `${ch.lo} – ${ch.hi}`;
      fill.style.left = `${100 * ch.lo / ch.top}%`; fill.style.right = `${100 - 100 * Math.min(ch.hi, ch.top) / ch.top}%`;
    };
    const lo = el("input", { type: "range", min: 0, max: ch.top, value: ch.lo, title: "Darkest value shown", "aria-label": `Minimum of channel ${c + 1}`,
      oninput: (e) => { ch.lo = Math.min(+e.target.value, ch.hi - 1); e.target.value = ch.lo; show(); syncLayers(); } });
    const hi = el("input", { type: "range", min: 0, max: ch.top, value: ch.hi, title: "Brightest value shown", "aria-label": `Maximum of channel ${c + 1}`,
      oninput: (e) => { ch.hi = Math.max(+e.target.value, ch.lo + 1); e.target.value = ch.hi; show(); syncLayers(); } });
    show();
    return el("div", { class: "channel" },
      el("div", { class: "channel-head" },
        el("input", { type: "checkbox", id: `ch${c}`, checked: ch.on, title: `Show channel ${c + 1} (key ${c + 1})`, onchange: (e) => { ch.on = e.target.checked; syncLayers(); } }),
        el("input", { type: "color", value: hex(ch.color), title: "Colour of this channel", "aria-label": `Colour of channel ${c + 1}`,
          oninput: (e) => { ch.color = unhex(e.target.value); store("fa.colors", { ...(store("fa.colors") || {}), [c]: e.target.value }); syncLayers(); } }),
        el("label", { htmlFor: `ch${c}`, textContent: `Ch ${c + 1}` }), c < 9 ? el("kbd", { textContent: c + 1 }) : null, values,
        el("button", { textContent: "Auto", title: "Contrast from the image (0.1 to 99.9 percentile)",
          onclick: () => { [ch.lo, ch.hi] = ch.auto; lo.value = ch.lo; hi.value = ch.hi; show(); syncLayers(); } })),
      el("div", { class: "dual" }, el("div", { class: "track" }), fill, lo, hi));
  }));
}
// per foci channel: show or hide its foci and markers
function renderOverlays() {
  $("foci-toggles").replaceChildren(...((S.settings.foci || []).length > 1 ? S.settings.foci : []).map((f, i) => el("label", {},
    el("input", { type: "checkbox", checked: V.fociOn[f.channel] !== false, onchange: (e) => { V.fociOn[f.channel] = e.target.checked; syncLayers(); } }),
    " ", el("span", { class: "dot", style: `background:${hex(OUTLINE[i % OUTLINE.length])}` }), `channel ${f.channel}`)));
}

document.addEventListener("keydown", (e) => {
  if ($("page-analyze").hidden || e.ctrlKey || e.metaKey || e.altKey) return;
  if (e.key === "Escape" && V.multi.size) { V.multi.clear(); selectCell(0, null); return; }
  const tag = e.target.tagName, type = e.target.type;
  if (tag === "SELECT" || tag === "TEXTAREA" || (tag === "INPUT" && !["checkbox", "range", "color"].includes(type))) return;
  const toggle = (id) => { $(id).checked = !$(id).checked; $(id).onchange({ target: $(id) }); };
  const k = { BracketLeft: "[", BracketRight: "]", Comma: ",", Period: "." }[e.code] || e.key.toLowerCase();
  const frames = V.info && V.info.frames > 1;
  if (k >= "1" && k <= "9" && V.ch[+k - 1]) { V.ch[+k - 1].on = !V.ch[+k - 1].on; renderChannels(); syncLayers(); }
  else if (k === "c") toggle("show-cells"); else if (k === "f") toggle("show-foci"); else if (k === "m") toggle("show-markers");
  else if (k === "o") toggle("show-coloc"); else if (k === "l") toggle("live"); else if (k === "p") nextMode(); else if (k === "a") runPreview();
  else if (k === "n") toggle("show-numbers"); else if (k === "s") toggle("smooth"); else if (k === "x") toggleCell();
  else if (k === "," || (k === "[" && !frames)) stepZ(-1); else if (k === "." || (k === "]" && !frames)) stepZ(1);
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
  if (b.running) status(`Analyzing all images: ${b.done} of ${b.total}${b.current ? ` · ${b.current}` : ""}`, false, b.total ? b.done / b.total : -1, true);
  $("progress-text").textContent = b.running ? `${b.done} of ${b.total} images${b.current ? ` · ${b.current}` : ""}`
    : b.errors && b.errors.length ? b.errors.join("\n") : "";
}
async function poll() {
  clearTimeout(pollTimer);
  try {
    const b = await (await api("/api/batch")).json(); showProgress(b);
    if (b.running) pollTimer = setTimeout(poll, 600);
    else if (b.stop) status("Stopped");
    else if (b.total) {
      status(`Analyzed ${b.done - (b.errors || []).length} of ${b.total} images`, !!(b.errors || []).length);
      refreshAnalyzed();
      if (!(b.errors || []).length) showPage("results");
    }
  } catch (e) { status(e.message, true); }
}
async function run(images) {
  try { showProgress(await (await api("/api/batch/start", { settings: S.settings, images, labels: $("save-labels").checked, images_overlay: $("save-images").checked })).json()); poll(); }
  catch (e) { status(e.message, true); }
}
$("save-labels").checked = store("fa.labels") === true;
$("save-labels").onchange = (e) => store("fa.labels", e.target.checked);
$("save-images").checked = store("fa.images") === true;
$("save-images").onchange = (e) => store("fa.images", e.target.checked);

// ---------------------------------------------------------------- settings files and saving pictures
// saving from the Save images card: the button shows a green tick for a moment, and a line under the buttons
// names the file, with a button that opens its folder
async function saved(what, request, button) {
  status(`Saving ${what}…`, false, -1);
  const label = button ? button.textContent : "";
  if (button) button.disabled = true;
  try {
    const r = await (await request()).json();
    status(r.cancelled ? "" : `Saved ${what}: ${r.path}`);
    if (button && !r.cancelled) {
      button.classList.add("saved"); button.textContent = "✓ Saved";
      setTimeout(() => { button.classList.remove("saved"); button.textContent = label; }, 1500);
      const file = String(r.path).split(/[\\/]/), sub = file.length > 1 ? file.at(-2) : "";
      $("save-note").hidden = false;
      $("save-note").replaceChildren(`Saved ${file.at(-1)}`, el("button", { textContent: "Show in folder", title: r.path,
        onclick: async () => { try { await api("/api/results/open", { sub }); } catch (e) { status(e.message, true); } } }));
    }
    return r;
  } catch (e) { status(e.message, true); return null; } finally { if (button) button.disabled = false; }
}
$("load-settings").onclick = async () => {
  try {
    const r = await (await api("/api/settings/load", {})).json();
    if (r.cancelled) return;
    S.settings = r.settings; renderSettings(); updateChannelColors();
    status(`Settings loaded from ${r.path}`);
    if (V.id && V.preview && $("live").checked) runPreview();
  } catch (e) { status(e.message, true); }
};
$("save-settings").onclick = () => saved("the settings", () => api("/api/settings/save", { settings: S.settings }));
function shownZ() {
  const det = DET[S.settings.foci_general.dimension];
  return V.mode === "slice" ? V.z : V.mode === "det" && det ? det[0] : "mip";
}
$("save-view").onclick = () => {
  if (!V.info) return status("Choose an image first", true);
  const name = `${S.images.find((im) => im.id === V.id)?.name || "view"}_view`;
  return saved("the view", () => api("/api/export/view", { png: viewer.snapshot(), name }), $("save-view"));
};
$("save-rgb").onclick = () => {
  if (!V.info) return status("Choose an image first", true);
  const p = V.preview;
  const overlays = p ? { cells: $("show-cells").checked, coloc: $("show-coloc").checked && (p.coloc || []).length > 0,
    foci: $("show-foci").checked ? p.channels.map((c) => c.channel).filter((c) => V.fociOn[c] !== false) : [] } : {};
  const display = V.ch.map((c) => ({ on: c.on, lo: c.lo, hi: c.hi, colour: c.color.map((v) => Math.round(255 * v)) }));
  return saved("the full image", () => api("/api/export/rgb", { path: V.id, t: V.t, z: shownZ(), display, overlays }), $("save-rgb"));
};
$("save-imagej").onclick = () => {
  if (!V.info) return status("Choose an image first", true);
  return saved("the image for Fiji", () => api("/api/export/imagej", { path: V.id, colours: V.ch.map((c) => c.color), ranges: V.ch.map((c) => [c.lo, c.hi]) }), $("save-imagej"));
};
$("run").onclick = () => run(null);
$("run-one").onclick = () => { if (V.id) run([V.id]); else status("Choose an image first", true); };

$("cancel").onclick = async () => {
  status("Stopping…", false, -1);
  try { await api("/api/cancel", {}); } catch (e) { status(e.message, true); }
};

// ---------------------------------------------------------------- settings page
async function loadProgram(change) {
  status("Asking the networks which versions they have…", false, -1);
  try {
    const p = await (await api("/api/program", change)).json();
    const inputs = {};
    S.installed = p.installed; S.gpu = p.gpu_ready; renderSettings();
    $("opt-gpu").checked = p.gpu !== false; $("opt-workers").value = p.workers || 0;
    const version = (n) => n.version ? n.version : !p.can_install ? "not installed"
      : el("button", { textContent: p.install.running && p.install.network === n.name ? "Installing…" : "Install", disabled: p.install.running,
        title: "Download and install this network into the folder of Foci Analyzer (a few GB; takes minutes)", onclick: () => installNetwork(n.name, n.label) });
    const device = (n) => n.error ? n.error : !n.version ? "–" : n.name === "stardist" ? "processor (TensorFlow)"
      : n.device === "cuda" ? `graphics card: ${n.gpu}${n.gpu_memory ? `, ${n.gpu_memory} GB` : ""}` : n.device === "cpu" ? "processor (no usable graphics card)" : "–";
    const general = el("input", { type: "text", value: p.python, placeholder: p.default_python || "no networks installed", "aria-label": "Python of the networks" });
    const apply = () => loadProgram({ python: general.value, pythons: Object.fromEntries(Object.entries(inputs).map(([k, i]) => [k, i.value])) });
    $("program").replaceChildren(
      el("table", {}, el("thead", {}, el("tr", {}, ["Network", "Version", "Runs on", "Python used", "Own Python for this network"].map((h) => el("th", { textContent: h })))),
        el("tbody", {}, p.networks.map((n) => el("tr", {},
          el("td", { textContent: n.label || n.name }), el("td", {}, version(n)),
          el("td", { class: n.error ? "bad" : "", textContent: device(n) }), el("td", { class: "path", textContent: n.python, title: n.python }),
          el("td", {}, inputs[n.name] = el("input", { type: "text", value: p.pythons[n.name] || "", placeholder: "the same as the others", "aria-label": `Python for ${n.name}` })))))),
      el("pre", { id: "program-log", hidden: !(p.install.running || p.install.ok === false) || !p.install.lines.length, textContent: p.install.lines.join("\n") }),
      el("div", { class: "row" }, el("label", { textContent: "Python of the networks" }), helpIcon("The python.exe of the environment that holds the networks, for example C:/PythonProjects/Cellpose4/Scripts/python.exe. Empty = the environment installed with Foci Analyzer. A conda environment works too (…/envs/name/python.exe). A network can have its own in the table above."),
        el("div", { style: "flex:1" }, general), el("button", { class: "primary", textContent: "Apply", onclick: apply })));
    general.style.width = "100%";
    $("program-about").replaceChildren(`Foci Analyzer ${p.version} · Python ${p.app_python} · the window keeps its state in ${p.home}`);
    status("");
  } catch (e) { status(e.message, true); }
}

async function installNetwork(name, label) {
  try { await api("/api/program/install", { network: name }); } catch (e) { status(e.message, true); return; }
  const poll = async () => {
    try {
      const st = await (await api("/api/program/install")).json();
      const log = $("program-log");
      if (log) { log.hidden = false; log.textContent = st.lines.join("\n"); log.scrollTop = log.scrollHeight; }
      if (st.running) { status(`Installing ${label}… ${st.lines.at(-1) || ""}`, false, -1); setTimeout(poll, 1500); return; }
      status(st.ok ? `${label} is installed` : `Installing ${label} failed: see the messages on the Settings page`, !st.ok);
      loadProgram();
    } catch (e) { status(e.message, true); }
  };
  loadProgram().then(poll);
}
$("opt-apply").onclick = () => loadProgram({ gpu: $("opt-gpu").checked, workers: parseInt($("opt-workers").value, 10) || 0 });

$("quit").onclick = async () => {
  if (S.batch.running && !window.confirm("An analysis is running. Stop Foci Analyzer anyway?")) return;
  try { await api("/api/quit", {}); } catch { /* it is stopping */ }
  document.body.replaceChildren(el("p", { style: "margin:60px auto;font-size:16px", textContent: "Foci Analyzer has stopped. You can close this tab." }));
};

// ---------------------------------------------------------------- updates (updates.py)
let U = null, updateLater = false;
function renderUpdate(u) {
  U = u;
  $("upd-auto").checked = u.check;
  const when = u.last_check ? new Date(u.last_check * 1000).toLocaleString() : "not yet";
  const what = u.newer ? `Version ${u.latest.version} is available.` : u.latest ? "This is the newest version." : u.last_check ? "No release was found on GitHub." : "";
  $("upd-info").replaceChildren(`You have Foci Analyzer ${u.current}. Last check: ${when}. ${what} `,
    u.newer ? el("button", { textContent: "Show the new version…", onclick: () => showUpdate(true) }) : "");
}
function showUpdate(asked) {
  if (!U || !U.newer || (!asked && (!U.show || updateLater))) return;
  const d = $("upd-dialog");
  const choose = async (choice) => {
    d.hidden = true;
    if (choice === "later") { updateLater = true; return; }
    try { renderUpdate(await (await api("/api/update/choice", { choice })).json()); } catch (e) { status(e.message, true); }
  };
  const go = U.can_update ? el("button", { class: "primary", textContent: "Update now", onclick: startUpdate })
    : el("a", { href: U.latest.page, target: "_blank", rel: "noopener", textContent: "Open the release page" });
  d.replaceChildren(el("div", { class: "dialog-box" },
    el("h2", { textContent: `Foci Analyzer ${U.latest.version} is available` }),
    el("p", { class: "muted", textContent: `You have ${U.current}.${U.latest.date ? ` Released ${U.latest.date}.` : ""} ` + (U.can_update
      ? "Updating takes a few minutes; Foci Analyzer then starts again. Your settings, experiment, results and networks stay."
      : "This copy was not installed with the setup (a developer's copy): update it with git pull.") }),
    el("div", { class: "notes", textContent: (U.latest.notes || "No release notes.").replace(/^#+\s*v?[\d.]+[^\n]*\n+/, "").replace(/^#+\s*/gm, "").replace(/\*\*/g, "") }),
    el("div", { class: "row" }, go,
      el("button", { textContent: "Later", title: "Ask again the next time Foci Analyzer starts", onclick: () => choose("later") }),
      el("button", { textContent: "Not for a week", onclick: () => choose("week") }),
      el("button", { textContent: "Skip this version", title: "Only a version newer than this one is offered again", onclick: () => choose("skip") }))));
  d.hidden = false;
}
async function loadUpdate() {
  try { renderUpdate(await (await api("/api/update")).json()); showUpdate(false); } catch { /* not important */ }
}
async function startUpdate() {
  $("upd-dialog").hidden = true;
  if (S.batch.running) { status("An analysis of all images is running: wait for it or stop it, then update.", true); return; }
  try { await api("/api/update/start", {}); } catch (e) { status(e.message, true); return; }
  const poll = async () => {
    let u;
    try { u = await (await api("/api/update")).json(); } catch { waitForUpdater(); return; }
    if (u.state === "failed") { status(`The update failed: ${u.error}`, true); return; }
    if (u.state === "restarting") { status("Foci Analyzer stops for the update…", false, -1); setTimeout(waitForUpdater, 1000); return; }
    status(`Downloading Foci Analyzer ${U.latest.version}…`, false, u.total ? u.done / u.total : -1);
    setTimeout(poll, 500);
  };
  poll();
}
function waitForUpdater() {  /* the page of the updater, or already the new window: anything but this program */
  fetch("/", { cache: "no-store" }).then((r) => r.text()).then((t) => {
    if (t.indexOf("(updating)") >= 0 || t.indexOf(window.TOKEN) < 0) location.reload(); else setTimeout(waitForUpdater, 1000);
  }).catch(() => setTimeout(waitForUpdater, 1000));
}
$("upd-auto").onchange = async () => { try { renderUpdate(await (await api("/api/update/choice", { check: $("upd-auto").checked })).json()); } catch (e) { status(e.message, true); } };
$("upd-check").onclick = async () => {
  status("Looking for a newer version on GitHub…", false, -1);
  try {
    renderUpdate(await (await api("/api/update/check", {})).json());
    status(U.newer ? `Foci Analyzer ${U.latest.version} is available` : "This is the newest version");
    if (U.newer) showUpdate(true);
  } catch (e) { status(e.message, true); }
};

// ---------------------------------------------------------------- start
window.FA = { viewer, V, S };  // for looking at the state from the browser console
initResults({
  openCell: async (id, frame, x, y, width) => {  // from a clicked point on the Results page
    showPage("analyze");
    if (V.id !== id) await openImage(id);
    if (frame) setT(frame - 1);
    if (x !== undefined) viewer.centerOn(x, y, width);
    if (!V.preview) runPreview();  // after centring: of a large image only the part around the cell is analyzed
  },
});
(async () => {
  try {
    Object.assign(S, await (await api("/api/state")).json());
    if (!S.conditions.length) S.conditions.push({ name: "Condition 1", folders: [], pattern: "*", recursive: false });
    renderConditions(); renderImages(); renderSettings(); syncMode(); showProgress(S.batch);
    if (S.batch.running) poll();
  } catch (e) { status(e.message, true); }
  loadUpdate(); setTimeout(loadUpdate, 20000);  /* the check at start runs in the background */
})();
