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
export function status(text, error = false) { const s = $("status"); s.textContent = text; s.classList.toggle("error", error); }

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
export function helpIcon(text) { return text ? el("span", { class: "help", tabIndex: 0, role: "button", "aria-label": "Help", help: text }) : el("span"); }
const pop = $("help-pop");
function showHelp(target) {
  pop.textContent = target.dataset.help; pop.hidden = false;
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
