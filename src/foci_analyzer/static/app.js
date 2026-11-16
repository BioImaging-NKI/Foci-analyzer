// Foci Analyzer window: pages, the experiment, the settings form (drawn from the schema) and a simple viewer.
// The viewer here draws whole planes on a 2D canvas; the WebGL tile viewer of plan.md section 6 replaces it.

const $ = (id) => document.getElementById(id);
function el(tag, props = {}, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (k === "class") e.className = v; else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else if (k in e) e[k] = v; else e.setAttribute(k, v);
  }
  for (const kid of kids.flat()) if (kid != null) e.append(kid);
  return e;
}
const clamp = (v, a, b) => Math.min(b, Math.max(a, v));

async function api(path, body) {
  const r = await fetch(path, {
    method: body ? "POST" : "GET",
    headers: { "X-Token": window.TOKEN, ...(body ? { "Content-Type": "application/json" } : {}) },
    body: body ? JSON.stringify(body) : undefined,
  });
  if (!r.ok) { let m = r.statusText; try { m = (await r.json()).error; } catch { /* not JSON */ } throw new Error(m); }
  return r;
}
function status(text, error = false) { const s = $("status"); s.textContent = text; s.classList.toggle("error", error); }

const S = { schema: [], settings: {}, conditions: [], images: [], remarks: [] };

// ---------------------------------------------------------------- pages
function showPage(name) {
  for (const b of $("tabs").children) b.classList.toggle("on", b.dataset.page === name);
  for (const p of document.querySelectorAll(".page")) p.hidden = p.id !== `page-${name}`;
  if (name === "analyze") { resizeCanvas(); draw(); }
}
$("tabs").addEventListener("click", (e) => { if (e.target.dataset.page) showPage(e.target.dataset.page); });

// ---------------------------------------------------------------- experiment
function renderConditions() {
  const box = $("conditions"); box.replaceChildren();
  S.conditions.forEach((c, i) => {
    const folders = el("div");
    c.folders.forEach((f, j) => folders.append(el("div", { class: "folder" },
      el("input", { type: "text", value: f, placeholder: "Folder", "aria-label": "Folder", oninput: (e) => { c.folders[j] = e.target.value; } }),
      el("button", { class: "link", textContent: "Remove", onclick: () => { c.folders.splice(j, 1); renderConditions(); } }))));
    box.append(el("div", { class: "condition" },
      el("div", { class: "row" },
        el("input", { type: "text", class: "name", value: c.name, placeholder: "Condition name", "aria-label": "Condition name", oninput: (e) => { c.name = e.target.value; } }),
        el("label", {}, "Files named ", el("input", { type: "text", class: "pattern", value: c.pattern, title: "Pattern for the file name, * = anything", oninput: (e) => { c.pattern = e.target.value; } })),
        el("label", {}, el("input", { type: "checkbox", checked: c.recursive, onchange: (e) => { c.recursive = e.target.checked; } }), " with sub-folders"),
        el("button", { class: "link", textContent: "Remove condition", onclick: () => { S.conditions.splice(i, 1); renderConditions(); } })),
      folders,
      el("button", { textContent: "Add folder", onclick: () => addFolder(c) })));
  });
}
async function addFolder(c) {
  let folder = "";
  try { folder = (await (await api("/api/pick-folder", {})).json()).folder; if (!folder) return; }
  catch (e) { status(e.message, true); }  // no dialog: an empty box to type or paste the folder in
  c.folders.push(folder); renderConditions();
}
function renderImages() {
  const body = $("images").tBodies[0]; body.replaceChildren();
  for (const im of S.images) body.append(el("tr", { onclick: () => { showPage("analyze"); openImage(im.path); } },
    el("td", { textContent: im.name }), el("td", { textContent: im.condition }), el("td", { class: "path", textContent: im.folder, title: im.folder })));
  $("image-count").textContent = S.images.length ? `(${S.images.length})` : "";
  $("remarks").replaceChildren(...S.remarks.map((r) => el("li", { textContent: r })));
  const sel = $("image-select"); sel.replaceChildren(el("option", { value: "", textContent: "Choose an image…" }));
  for (const im of S.images) sel.append(el("option", { value: im.path, textContent: `${im.condition} · ${im.name}` }));
  sel.value = S.images.some((im) => im.path === V.path) ? V.path : "";
}
$("add-condition").onclick = () => { S.conditions.push({ name: `Condition ${S.conditions.length + 1}`, folders: [], pattern: "*", recursive: false }); renderConditions(); };
$("scan").onclick = async () => {
  try {
    Object.assign(S, await (await api("/api/experiment", { conditions: S.conditions })).json());
    renderConditions(); renderImages(); status(`${S.images.length} images`);
  } catch (e) { status(e.message, true); }
};

// ---------------------------------------------------------------- settings form
let saveTimer = 0;
function saveSettings() {
  clearTimeout(saveTimer);
  saveTimer = setTimeout(async () => {
    try { S.settings = await (await api("/api/settings", { settings: S.settings })).json(); status("Settings saved"); updateChannels(); }
    catch (e) { status(e.message, true); }
  }, 300);
}
function control(f, group) {
  const set = (v) => { group[f.key] = v; saveSettings(); };
  const id = `f${Math.random().toString(36).slice(2)}`, v = group[f.key];
  let input;
  if (f.type === "bool") input = el("input", { id, type: "checkbox", checked: v, onchange: (e) => set(e.target.checked) });
  else if (f.type === "choice") { input = el("select", { id, onchange: (e) => set(e.target.value) }, f.choices.map((c) => el("option", { value: c, textContent: c }))); input.value = v; }
  else if (f.type === "str") input = el("input", { id, type: "text", value: v, onchange: (e) => set(e.target.value) });
  else input = el("input", { id, type: "number", value: v, min: f.min ?? "", max: f.max ?? "", step: f.step ?? (f.type === "int" ? 1 : "any"),
    onchange: (e) => { const n = f.type === "int" ? parseInt(e.target.value, 10) : parseFloat(e.target.value); if (Number.isFinite(n)) set(n); else e.target.value = group[f.key]; } });
  return el("div", { class: `field${f.type === "bool" ? " check" : ""}`, title: f.help || "" },
    el("label", { htmlFor: id, textContent: f.label }),
    el("div", { class: "control" }, input, f.unit ? el("span", { class: "unit", textContent: f.unit }) : null));
}
function fieldsOf(section, group) {
  const basic = section.fields.filter((f) => !f.advanced), adv = section.fields.filter((f) => f.advanced);
  return [...basic.map((f) => control(f, group)),
    adv.length ? el("details", { class: "advanced" }, el("summary", { textContent: "Advanced" }), adv.map((f) => control(f, group))) : null];
}
function renderSettings() {
  const box = $("settings"); box.replaceChildren();
  for (const sec of S.schema) {
    const div = el("div", { class: "section" });
    if (!sec.repeat) { div.append(el("h3", { textContent: sec.title }), ...fieldsOf(sec, S.settings[sec.key]).filter(Boolean)); box.append(div); continue; }
    const list = S.settings[sec.key];
    div.append(el("h3", {}, `${sec.title}s`, el("button", { textContent: "Add", onclick: () => {
      const g = Object.fromEntries(sec.fields.map((f) => [f.key, f.default]));
      if ("channel" in g) g.channel = Math.max(...list.map((x) => x.channel)) + 1;
      list.push(g); renderSettings(); saveSettings(); } })));
    list.forEach((g, i) => div.append(el("div", { class: "group" },
      el("div", { class: "group-head" }, `${sec.title} ${i + 1}`,
        list.length > 1 ? el("button", { class: "link", textContent: "Remove", onclick: () => { list.splice(i, 1); renderSettings(); saveSettings(); } }) : null),
      ...fieldsOf(sec, g).filter(Boolean))));
    box.append(div);
  }
}

// ---------------------------------------------------------------- viewer (simple; see plan.md section 6)
const canvas = $("canvas"), wrap = $("canvas-wrap");
const V = { path: "", info: null, mode: "mip", z: 0, ch: [], planes: new Map(), view: { s: 1, tx: 0, ty: 0 }, off: null, load: 0 };
const NUCLEI = [0, 0.5, 1], CELLS = [1, 0.5, 0], OTHER = [[0, 1, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]];

function channelColor(c) {
  if (c + 1 === S.settings.image.nuclei_channel) return NUCLEI;
  if (c + 1 === S.settings.image.cell_channel) return CELLS;
  const foci = S.settings.foci.map((f) => f.channel - 1), i = foci.indexOf(c);
  return OTHER[(i >= 0 ? i : foci.length + c) % OTHER.length];
}
function updateChannels() { V.ch.forEach((ch, c) => { ch.color = channelColor(c); }); renderChannels(); compose(); }

async function openImage(path) {
  const load = ++V.load;
  $("image-select").value = path;
  if (!path) return;
  try {
    status("Reading image…");
    const info = await (await api(`/api/image/info?path=${encodeURIComponent(path)}`)).json();
    if (load !== V.load) return;
    Object.assign(V, { path, info, z: Math.floor(info.slices / 2), planes: new Map() });
    V.ch = info.ranges.map(([lo, hi], c) => ({ on: true, lo, hi, auto: [lo, hi], color: channelColor(c) }));
    if (info.slices === 1) V.mode = "mip";
    $("mode").hidden = info.slices === 1;
    $("z").max = info.slices - 1; $("z").value = V.z;
    $("viewer-empty").hidden = true;
    renderChannels(); syncMode(); await refresh(); fit(); status("");
  } catch (e) { status(e.message, true); }
}
async function plane(c) {
  const z = V.mode === "mip" ? "mip" : V.z, key = `${c}|${z}`;
  if (!V.planes.has(key)) {
    const r = await api(`/api/image/region?path=${encodeURIComponent(V.path)}&t=0&c=${c}&z=${z}&level=0`);
    V.planes.set(key, new Uint16Array(await r.arrayBuffer()));
    while (V.planes.size > 80) V.planes.delete(V.planes.keys().next().value);
  }
  return V.planes.get(key);
}
async function refresh() {
  if (!V.info) return;
  const load = V.load, want = `${V.mode}|${V.z}`;
  try {
    const data = await Promise.all(V.ch.map((ch, c) => (ch.on ? plane(c) : null)));
    if (load !== V.load || want !== `${V.mode}|${V.z}`) return;  // a newer request is on its way
    V.data = data; compose();
  } catch (e) { status(e.message, true); }
}
function compose() {
  if (!V.info || !V.data) return;
  const { width: w, height: h } = V.info, n = w * h;
  if (!V.off || V.off.width !== w || V.off.height !== h) V.off = new OffscreenCanvas(w, h);
  const ctx = V.off.getContext("2d"), img = ctx.createImageData(w, h), out = img.data;
  const acc = new Float32Array(n * 3);
  V.ch.forEach((ch, c) => {
    const d = V.data[c]; if (!ch.on || !d) return;
    const k = 255 / Math.max(1, ch.hi - ch.lo), [r, g, b] = ch.color;
    for (let p = 0, q = 0; p < n; p++, q += 3) {
      const v = clamp((d[p] - ch.lo) * k, 0, 255);
      acc[q] += v * r; acc[q + 1] += v * g; acc[q + 2] += v * b;
    }
  });
  for (let p = 0, q = 0, o = 0; p < n; p++, q += 3, o += 4) { out[o] = acc[q]; out[o + 1] = acc[q + 1]; out[o + 2] = acc[q + 2]; out[o + 3] = 255; }
  ctx.putImageData(img, 0, 0); draw();
}
function resizeCanvas() {
  const r = wrap.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
  const W = Math.round(r.width * dpr), H = Math.round(r.height * dpr);
  if (canvas.width !== W || canvas.height !== H) { canvas.width = W; canvas.height = H; }
}
function draw() {
  const ctx = canvas.getContext("2d"), dpr = window.devicePixelRatio || 1;
  ctx.setTransform(1, 0, 0, 1, 0, 0); ctx.clearRect(0, 0, canvas.width, canvas.height);
  if (!V.off || !V.info) return;
  const { s, tx, ty } = V.view;
  ctx.imageSmoothingEnabled = s < 1; ctx.imageSmoothingQuality = "high";  // real pixels when zoomed in
  ctx.setTransform(s * dpr, 0, 0, s * dpr, tx * dpr, ty * dpr); ctx.drawImage(V.off, 0, 0);
  const i = V.info, cal = i.pixel_size ? ` · ${i.pixel_size.toFixed(3)} ${i.unit}/px` : " · not calibrated";
  $("view-info").textContent = `${i.width} × ${i.height} px · ${i.slices} slice${i.slices > 1 ? "s" : ""}${cal} · zoom ${Math.round(s * 100)} %`;
}
function fit() {
  if (!V.info) return;
  resizeCanvas();
  const r = wrap.getBoundingClientRect(), s = Math.min(r.width / V.info.width, r.height / V.info.height);
  V.view = { s, tx: (r.width - V.info.width * s) / 2, ty: (r.height - V.info.height * s) / 2 }; draw();
}
// zoom glides to its target around the cursor, so a wheel notch does not jump
let glide = null;
function zoomAt(px, py, factor) {
  if (!V.info) return;
  const r = wrap.getBoundingClientRect(), smin = Math.min(r.width / V.info.width, r.height / V.info.height) / 2;
  const target = clamp((glide ? glide.s : V.view.s) * factor, smin, 40);
  const start = !glide; glide = { s: target, px, py, t: performance.now() };
  if (start) requestAnimationFrame(glideStep);
}
function glideStep(now) {
  if (!glide) return;
  const k = 1 - Math.exp(-(now - glide.t) / 35); glide.t = now;
  let s = V.view.s * Math.pow(glide.s / V.view.s, k);
  const done = Math.abs(Math.log(glide.s / s)) < 0.002; if (done) s = glide.s;
  const ix = (glide.px - V.view.tx) / V.view.s, iy = (glide.py - V.view.ty) / V.view.s;
  V.view = { s, tx: glide.px - ix * s, ty: glide.py - iy * s }; draw();
  if (done) glide = null; else requestAnimationFrame(glideStep);
}
wrap.addEventListener("wheel", (e) => {
  e.preventDefault(); const r = wrap.getBoundingClientRect();
  zoomAt(e.clientX - r.left, e.clientY - r.top, Math.pow(1.33, -e.deltaY / 100));
}, { passive: false });
let drag = null;
canvas.addEventListener("pointerdown", (e) => { drag = { x: e.clientX, y: e.clientY }; canvas.setPointerCapture(e.pointerId); canvas.classList.add("dragging"); });
canvas.addEventListener("pointermove", (e) => { if (!drag) return; V.view.tx += e.clientX - drag.x; V.view.ty += e.clientY - drag.y; drag = { x: e.clientX, y: e.clientY }; draw(); });
canvas.addEventListener("pointerup", () => { drag = null; canvas.classList.remove("dragging"); });
new ResizeObserver(() => { resizeCanvas(); draw(); }).observe(wrap);

function syncMode() {
  for (const b of $("mode").children) b.classList.toggle("on", b.dataset.mode === V.mode);
  $("z-wrap").hidden = V.mode !== "slice" || !V.info || V.info.slices === 1;
  $("z-label").textContent = V.info ? `${V.z + 1} / ${V.info.slices}` : "";
}
$("mode").addEventListener("click", (e) => { if (!e.target.dataset.mode) return; V.mode = e.target.dataset.mode; syncMode(); refresh(); });
$("z").addEventListener("input", (e) => { V.z = +e.target.value; syncMode(); refresh(); });
$("fit").onclick = fit;
$("image-select").onchange = (e) => openImage(e.target.value);

function renderChannels() {
  $("channels").replaceChildren(...V.ch.map((ch, c) => {
    const top = Math.max(ch.auto[1] * 3, ch.auto[0] + 10);
    return el("div", { class: "channel" },
      el("input", { type: "checkbox", id: `ch${c}`, checked: ch.on, onchange: (e) => { ch.on = e.target.checked; refresh(); } }),
      el("span", { class: "swatch", style: `background: rgb(${ch.color.map((v) => v * 255).join(",")})` }),
      el("label", { htmlFor: `ch${c}`, textContent: `Channel ${c + 1}` }),
      el("input", { type: "range", min: ch.auto[0] + 1, max: top, value: ch.hi, title: "Brightness (double click: automatic)", "aria-label": `Brightness channel ${c + 1}`,
        oninput: (e) => { ch.hi = +e.target.value; compose(); }, ondblclick: (e) => { ch.hi = ch.auto[1]; e.target.value = ch.hi; compose(); } }));
  }));
}

// ---------------------------------------------------------------- start
(async () => {
  try {
    Object.assign(S, await (await api("/api/state")).json());
    if (!S.conditions.length) S.conditions.push({ name: "Condition 1", folders: [], pattern: "*", recursive: false });
    renderConditions(); renderImages(); renderSettings(); syncMode();
  } catch (e) { status(e.message, true); }
})();
