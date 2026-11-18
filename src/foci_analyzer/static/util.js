// Small helpers shared by the pages: elements, requests, the status line, context help and the theme.

export const $ = (id) => document.getElementById(id);
export function el(tag, props = {}, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (k === "class") e.className = v; else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else if (k === "help") { if (v) e.dataset.help = v; }
    else if (k in e) e[k] = v; else e.setAttribute(k, v);
  }
  for (const kid of kids.flat()) if (kid != null) e.append(kid);
  return e;
}
export const clamp = (v, a, b) => Math.min(b, Math.max(a, v));

export async function api(path, body) {
  const r = await fetch(path, {
    method: body ? "POST" : "GET",
    headers: { "X-Token": window.TOKEN, ...(body ? { "Content-Type": "application/json" } : {}) },
    body: body ? JSON.stringify(body) : undefined,
  });
  if (!r.ok) { let m = r.statusText; try { m = (await r.json()).error; } catch { /* not JSON */ } throw new Error(m); }
  return r;
}
// The status bar at the bottom. `progress`: undefined = not busy, -1 = busy for an unknown time, 0..1 = how far.
// The bar only grows while one task runs; when it ends it fills up and then empties at once (no sliding back).
let barRun = 0;
function barJump(bar, width) { bar.style.transition = "none"; bar.style.width = width; void bar.offsetWidth; bar.style.transition = ""; }
export function status(text, error = false, progress, cancellable = false) {
  const s = $("status"), busy = progress !== undefined, bar = $("busy-bar"), wasBusy = !$("busy").classList.contains("idle"), wasUnknown = $("busy").classList.contains("unknown");
  s.textContent = text || (busy ? "Working…" : "Ready"); s.title = s.textContent;
  s.classList.toggle("error", error); s.classList.toggle("busy", busy && !error);
  $("busy").classList.toggle("idle", !busy);  /* always there, so the text never moves */ $("cancel").hidden = !(busy && cancellable); $("busy").classList.toggle("unknown", progress === -1);
  const run = ++barRun, shown = parseFloat(bar.style.width) || 0;
  if (busy && (!wasBusy || (wasUnknown && progress >= 0))) barJump(bar, "0%");  // a new task, or the first known step: from empty
  if (busy && progress >= 0) bar.style.width = `${Math.round(100 * Math.min(1, Math.max(progress, wasBusy ? shown / 100 : 0)))}%`;
  else if (busy) bar.style.width = "";
  else if (wasBusy && shown > 0 && !error) {
    bar.style.width = "100%";
    setTimeout(() => { if (run === barRun) barJump(bar, ""); }, 400);
  } else barJump(bar, "");
}

// a number as people read it: 3 significant digits, no exponent for ordinary sizes
export function fmt(v) {
  if (v === null || v === undefined || Number.isNaN(v)) return "–";
  if (typeof v !== "number") return String(v);
  if (Number.isInteger(v)) return String(v);
  const a = Math.abs(v);
  if (a !== 0 && (a < 0.001 || a >= 1e6)) return v.toExponential(2);
  return String(+v.toPrecision(3));
}
export function pvalue(p) { return p === null || p === undefined ? "–" : p < 0.0001 ? "< 0.0001" : String(+p.toPrecision(2)); }

// ---- context help: any element with class "help" and a data-help text
export function helpIcon(text, link) {
  if (!text) return el("span");
  const h = el("span", { class: "help", tabIndex: 0, role: "button", "aria-label": "Help", help: text });
  if (link) h.dataset.link = link;
  return h;
}
const pop = $("help-pop");
function showHelp(target) {
  pop.textContent = target.dataset.help; pop.hidden = false;
  if (target.dataset.link) pop.append(" ", el("a", { href: target.dataset.link, target: "_blank", rel: "noopener", textContent: "More about these models ↗" }));
  const r = target.getBoundingClientRect(), w = pop.offsetWidth, h = pop.offsetHeight;
  pop.style.left = `${clamp(r.left, 8, window.innerWidth - w - 8)}px`;
  pop.style.top = `${r.bottom + h + 12 > window.innerHeight ? Math.max(8, r.top - h - 6) : r.bottom + 6}px`;
}
document.addEventListener("click", (e) => {
  const t = e.target.closest ? e.target.closest(".help") : null;
  if (t && t.dataset.help) { e.preventDefault(); e.stopPropagation(); showHelp(t); } else pop.hidden = true;
});
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape") pop.hidden = true;
  if ((e.key === "Enter" || e.key === " ") && e.target.classList && e.target.classList.contains("help")) { e.preventDefault(); showHelp(e.target); }
});

// ---- theme: auto (as the computer), light or dark
const THEMES = ["auto", "light", "dark"];
export function theme() { try { return localStorage.getItem("fa.theme") || "auto"; } catch { return "auto"; } }
export function isDark() { const t = theme(); return t === "dark" || (t === "auto" && window.matchMedia("(prefers-color-scheme: dark)").matches); }
function applyTheme(t) {
  if (t === "auto") delete document.documentElement.dataset.theme; else document.documentElement.dataset.theme = t;
  $("theme").textContent = `Theme: ${t}`;
  document.dispatchEvent(new Event("themechange"));
}
$("theme").onclick = () => {
  const t = THEMES[(THEMES.indexOf(theme()) + 1) % 3];
  try { localStorage.setItem("fa.theme", t); } catch { /* private window */ }
  applyTheme(t);
};
applyTheme(theme());

export function store(key, value) {
  try {
    if (value === undefined) return JSON.parse(localStorage.getItem(key) || "null");
    localStorage.setItem(key, JSON.stringify(value));
  } catch { /* private window: nothing is remembered */ }
  return null;
}
