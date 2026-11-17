// The image viewer: WebGL2, tiles from a pyramid, channels composed in the shader (plan.md, section 6).
// Layers are images (16-bit tiles, contrast and colour applied in the shader, so they change without a
// request) and label images (drawn as outlines). Markers and the scale bar are drawn on a 2D canvas on top.

const TILE = 512, MAX_TILES = 320, MAX_REQUESTS = 6;

const VS = `#version 300 es
in vec2 a; uniform vec4 rect; uniform vec4 uvr; out vec2 uv;
void main() { uv = uvr.xy + a * uvr.zw; gl_Position = vec4(rect.xy + a * rect.zw, 0.0, 1.0); }`;
const FS_IMAGE = `#version 300 es
precision highp float; precision highp usampler2D;
uniform usampler2D tex; uniform vec2 range; uniform vec3 color; in vec2 uv; out vec4 o;
void main() {
  float v = float(texture(tex, uv).r);
  o = vec4(color * clamp((v - range.x) / max(range.y - range.x, 1.0), 0.0, 1.0), 1.0);
}`;
const FS_LABELS = `#version 300 es
precision highp float; precision highp usampler2D;
uniform usampler2D tex; uniform vec2 d; uniform vec4 color; uniform float fill; uniform float hatch; in vec2 uv; out vec4 o;
void main() {
  uint c = texture(tex, uv).r;
  if (c == 0u) discard;
  bool edge = texture(tex, uv + vec2(d.x, 0.0)).r != c || texture(tex, uv - vec2(d.x, 0.0)).r != c
           || texture(tex, uv + vec2(0.0, d.y)).r != c || texture(tex, uv - vec2(0.0, d.y)).r != c;
  if (!edge && fill <= 0.0) discard;
  if (!edge && hatch > 0.0) {  // black and white stripes: visible on any colour
    float s = mod(floor((gl_FragCoord.x + gl_FragCoord.y) / hatch), 2.0);
    o = vec4(vec3(s), fill);
    return;
  }
  o = edge ? color : vec4(color.rgb, color.a * fill);  // the inside lightly filled, the outline solid
}`;

// cells coloured by a measure: values 1..255 through a colour map (dark blue, teal, yellow)
const FS_VALUE = `#version 300 es
precision highp float; precision highp usampler2D;
uniform usampler2D tex; uniform float alpha; in vec2 uv; out vec4 o;
void main() {
  uint v = texture(tex, uv).r;
  if (v == 0u) discard;
  float t = float(v - 1u) / 254.0;
  vec3 a = vec3(0.23, 0.10, 0.55), b = vec3(0.13, 0.62, 0.55), c = vec3(0.99, 0.91, 0.14);
  o = vec4(t < 0.5 ? mix(a, b, t * 2.0) : mix(b, c, t * 2.0 - 1.0), alpha);
}`;

function program(gl, fs) {
  const p = gl.createProgram();
  for (const [type, src] of [[gl.VERTEX_SHADER, VS], [gl.FRAGMENT_SHADER, fs]]) {
    const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    gl.attachShader(p, s);
  }
  gl.bindAttribLocation(p, 0, "a"); gl.linkProgram(p);
  const u = {};
  for (let i = 0; i < gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS); i++) {
    const name = gl.getActiveUniform(p, i).name; u[name] = gl.getUniformLocation(p, name);
  }
  return { p, u };
}

export class Viewer {
  // fetchTile(url) → Promise<Response> with the headers X-Width and X-Height and 16-bit pixels as body
  constructor(canvas, top, fetchTile) {
    this.canvas = canvas; this.top = top; this.fetchTile = fetchTile;
    this.W = 0; this.H = 0; this.view = { s: 1, tx: 0, ty: 0 };
    this.layers = []; this.markers = []; this.zslice = null; this.box = null; this.pixelSize = null; this.unit = "µm";
    this.onclick = null; this.onview = null;
    this.tiles = new Map(); this.queue = []; this.busy = 0; this.wanted = new Set(); this.frame = 0; this.dirty = false;
    const gl = this.gl = canvas.getContext("webgl2", { alpha: false, antialias: false });
    if (!gl) throw new Error("this browser has no WebGL2");
    this.pImage = program(gl, FS_IMAGE); this.pLabels = program(gl, FS_LABELS); this.pValue = program(gl, FS_VALUE);
    const buf = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1); gl.enable(gl.BLEND);
    this._events();
    new ResizeObserver(() => this.draw()).observe(canvas.parentElement);
  }

  setImage(width, height) { this.W = width; this.H = height; this.maxLevel = Math.min(8, Math.max(0, Math.ceil(Math.log2(Math.max(width, height) / TILE)))); }
  requestDraw() {
    if (this.dirty) return;
    this.dirty = true;
    const go = () => { if (this.dirty) { this.dirty = false; this.draw(); } };
    requestAnimationFrame(go); setTimeout(go, 150);  // the timer: a browser may hold back animation frames
  }

  // ---- view
  size() { const r = this.canvas.parentElement.getBoundingClientRect(); return [r.width, r.height]; }
  fit() {
    if (!this.W) return;
    const [w, h] = this.size(), s = Math.min(w / this.W, h / this.H);
    this.view = { s, tx: (w - this.W * s) / 2, ty: (h - this.H * s) / 2 }; this.requestDraw();
  }
  centerOn(x, y, widthPx) {  // show `widthPx` image pixels around (x, y)
    const [w, h] = this.size(), s = Math.min(40, Math.min(w, h) / widthPx);
    this.view = { s, tx: w / 2 - x * s, ty: h / 2 - y * s }; this.requestDraw();
  }
  zoomAt(px, py, factor) {  // glides to its target around the cursor, so a wheel notch does not jump
    if (!this.W) return;
    const [w, h] = this.size(), smin = Math.min(w / this.W, h / this.H) / 2;
    const target = Math.min(40, Math.max(smin, (this.glide ? this.glide.s : this.view.s) * factor));
    const start = !this.glide; this.glide = { s: target, px, py, t: null };
    if (start) requestAnimationFrame((now) => this._glide(now));
  }
  _glide(now) {
    const g = this.glide; if (!g) return;
    const dt = g.t === null ? 16 : Math.max(0, now - g.t); g.t = now;  // never a negative step: that zooms the wrong way
    const v = this.view; let s = v.s * Math.pow(g.s / v.s, 1 - Math.exp(-dt / 35));
    const done = Math.abs(Math.log(g.s / s)) < 0.002; if (done) s = g.s;
    const ix = (g.px - v.tx) / v.s, iy = (g.py - v.ty) / v.s;
    this.view = { s, tx: g.px - ix * s, ty: g.py - iy * s }; this.draw();
    if (done) this.glide = null; else requestAnimationFrame((t) => this._glide(t));
  }
  _events() {
    const top = this.top; let drag = null;
    top.addEventListener("wheel", (e) => {
      e.preventDefault(); const r = top.getBoundingClientRect();
      this.zoomAt(e.clientX - r.left, e.clientY - r.top, Math.pow(1.33, -e.deltaY / 100));
    }, { passive: false });
    top.addEventListener("pointerdown", (e) => { drag = { x: e.clientX, y: e.clientY, moved: 0 }; top.setPointerCapture(e.pointerId); top.classList.add("dragging"); });
    top.addEventListener("pointermove", (e) => {
      if (!drag) return;
      const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
      drag.moved += Math.abs(dx) + Math.abs(dy); drag.x = e.clientX; drag.y = e.clientY;
      this.view.tx += dx; this.view.ty += dy; this.requestDraw();
    });
    top.addEventListener("pointerup", (e) => {
      top.classList.remove("dragging");
      if (drag && drag.moved < 4 && this.onclick && this.W) {
        const r = top.getBoundingClientRect();
        this.onclick((e.clientX - r.left - this.view.tx) / this.view.s, (e.clientY - r.top - this.view.ty) / this.view.s);
      }
      drag = null;
    });
  }

  // ---- tiles
  _tile(layer, key, L, tx, ty) {
    const id = `${key}|${L}|${tx}|${ty}`; this.wanted.add(id);
    let t = this.tiles.get(id);
    if (!t && layer) { t = { ready: false, used: this.frame }; this.tiles.set(id, t); this.queue.push({ id, layer, L, tx, ty }); this._pump(); }
    if (t) t.used = this.frame;
    return t && t.ready ? t : null;
  }
  _pump() {
    while (this.busy < MAX_REQUESTS && this.queue.length) {
      const job = this.queue.pop();  // newest first: what is on screen now
      if (!this.wanted.has(job.id)) { this.tiles.delete(job.id); continue; }
      this.busy++;
      this.fetchTile(job.layer.url(job.L, job.tx * TILE, job.ty * TILE, TILE, TILE)).then(async (r) => {
        const w = +r.headers.get("X-Width"), h = +r.headers.get("X-Height"), data = new Uint16Array(await r.arrayBuffer());
        const gl = this.gl, tex = gl.createTexture(), t = this.tiles.get(job.id);
        if (!t) return;
        gl.bindTexture(gl.TEXTURE_2D, tex);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.R16UI, w, h, 0, gl.RED_INTEGER, gl.UNSIGNED_SHORT, data);
        for (const p of [gl.TEXTURE_MIN_FILTER, gl.TEXTURE_MAG_FILTER]) gl.texParameteri(gl.TEXTURE_2D, p, gl.NEAREST);
        for (const p of [gl.TEXTURE_WRAP_S, gl.TEXTURE_WRAP_T]) gl.texParameteri(gl.TEXTURE_2D, p, gl.CLAMP_TO_EDGE);
        Object.assign(t, { tex, w, h, ready: true }); this.requestDraw();
      }).catch(() => { this.tiles.delete(job.id); }).finally(() => { this.busy--; this._pump(); });
    }
  }
  _trim() {
    if (this.tiles.size <= MAX_TILES) return;
    const old = [...this.tiles].filter(([, t]) => t.used < this.frame && t.ready).sort((a, b) => a[1].used - b[1].used);
    for (const [id, t] of old.slice(0, this.tiles.size - MAX_TILES)) { this.gl.deleteTexture(t.tex); this.tiles.delete(id); }
  }

  // ---- drawing
  draw() {
    const gl = this.gl, dpr = window.devicePixelRatio || 1, [cw, ch] = this.size();
    const W = Math.max(1, Math.round(cw * dpr)), H = Math.max(1, Math.round(ch * dpr));
    for (const c of [this.canvas, this.top]) if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
    gl.viewport(0, 0, W, H); gl.clearColor(0.043, 0.051, 0.071, 1); gl.clear(gl.COLOR_BUFFER_BIT);
    this.frame++; this.wanted = new Set();
    if (this.W) {
      const { s, tx, ty } = this.view;
      const L = Math.min(this.maxLevel, Math.max(0, Math.floor(Math.log2(1 / (s * dpr)))));
      const f = 1 << L, lw = Math.ceil(this.W / f), lh = Math.ceil(this.H / f);
      const x0 = Math.max(0, Math.floor(-tx / s / f / TILE)), y0 = Math.max(0, Math.floor(-ty / s / f / TILE));
      const x1 = Math.min(Math.ceil(lw / TILE) - 1, Math.floor((cw - tx) / s / f / TILE));
      const y1 = Math.min(Math.ceil(lh / TILE) - 1, Math.floor((ch - ty) / s / f / TILE));
      for (const layer of this.layers) {
        if (!layer.on) continue;
        const image = layer.type === "image", value = layer.type === "value", pr = image ? this.pImage : value ? this.pValue : this.pLabels;
        gl.useProgram(pr.p); gl.uniform1i(pr.u.tex, 0);
        if (image) { gl.blendFunc(gl.ONE, gl.ONE); gl.uniform2f(pr.u.range, layer.lo, layer.hi); gl.uniform3fv(pr.u.color, layer.color); }
        else if (value) { gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); gl.uniform1f(pr.u.alpha, layer.alpha ?? 0.6); }
        else { gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA); gl.uniform4fv(pr.u.color, layer.color); gl.uniform1f(pr.u.fill, layer.fill || 0); gl.uniform1f(pr.u.hatch, (layer.hatch || 0) * dpr); }
        this._tile(layer, layer.key, this.maxLevel, 0, 0);  // the coarsest level is always kept as a stand-in
        for (let ty_ = y0; ty_ <= y1; ty_++) for (let tx_ = x0; tx_ <= x1; tx_++) {
          const tw = Math.min(TILE, lw - tx_ * TILE), th = Math.min(TILE, lh - ty_ * TILE);
          let t = this._tile(layer, layer.key, L, tx_, ty_), uv = [0, 0, 1, 1], texel = f;
          if (!t && layer.prevKey) t = this._tile(null, layer.prevKey, L, tx_, ty_);  // the previous result until the new one is there
          for (let l = L + 1; !t && l <= this.maxLevel; l++) {  // or a coarser tile that is loaded
            const k = 1 << (l - L), ax = Math.floor(tx_ / k), ay = Math.floor(ty_ / k);
            for (const key of [layer.key, layer.prevKey]) {
              const a = key && this._tile(null, key, l, ax, ay);
              if (a) { t = a; texel = f * k; uv = [(tx_ * TILE / k - ax * TILE) / a.w, (ty_ * TILE / k - ay * TILE) / a.h, tw / k / a.w, th / k / a.h]; break; }
            }
          }
          if (!t) continue;
          const ix = tx_ * TILE * f, iy = ty_ * TILE * f, iw = Math.min(tw * f, this.W - ix), ih = Math.min(th * f, this.H - iy);
          gl.uniform4f(pr.u.rect, (tx + ix * s) / cw * 2 - 1, 1 - (ty + iy * s) / ch * 2, iw * s / cw * 2, -ih * s / ch * 2);
          gl.uniform4f(pr.u.uvr, ...uv);
          if (!image && !value) { const px = Math.min((layer.width || 1) * 1.5 / (texel * s * dpr), 1); /* a thin line inside the edge pixel when zoomed in */ gl.uniform2f(pr.u.d, px / t.w, px / t.h); }
          gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, t.tex); gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
        }
      }
    }
    this._trim(); this._top(dpr, cw, ch);
    if (this.onview) this.onview();
  }
  // a PNG of what the viewer shows: the image layers and the markers and scale bar on top
  snapshot() {
    this.draw();  // the WebGL picture is only readable in the same task as it was drawn
    const out = document.createElement("canvas");
    out.width = this.canvas.width; out.height = this.canvas.height;
    const ctx = out.getContext("2d");
    ctx.drawImage(this.canvas, 0, 0); ctx.drawImage(this.top, 0, 0);
    return out.toDataURL("image/png");
  }
  _top(dpr, cw, ch) {
    const ctx = this.top.getContext("2d"); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, cw, ch);
    if (!this.W) return;
    // dots that shrink a little when zoomed out, so they do not cover the image
    const { s, tx, ty } = this.view, r = Math.max(1, Math.min(2.6, 1.3 * Math.sqrt(s)));
    for (const m of this.markers) {
      if (!m.on || !m.pts) continue;
      ctx.fillStyle = m.color; ctx.beginPath(); const rings = new Path2D();
      const p = m.pts, n = p.length / 6; let shown = 0;
      for (let i = 0; i < n; i++) {
        const o = i * 6;
        if (p[o + 3] < m.threshold) break;  // sorted from high to low
        if (this.zslice !== null && p[o + 2] >= 0 && Math.abs(p[o + 2] - this.zslice) > 1) continue;
        const x = tx + (p[o] + 0.5) * s, y = ty + (p[o + 1] + 0.5) * s;
        if (x < -4 || y < -4 || x > cw + 4 || y > ch + 4) continue;
        if (++shown > 60000) break;
        ctx.moveTo(x + r, y); ctx.arc(x, y, r, 0, 6.2832);
        if (p[o + 5] === 2 && s >= 1.5) { /* only when zoomed in: at a distance the rings would hide the image */ rings.moveTo(x + r + 1.3, y); rings.arc(x, y, r + 1.3, 0, 6.2832); }  // colocalizes
      }
      ctx.fill(); ctx.strokeStyle = "#000"; ctx.lineWidth = 2.4; ctx.stroke(rings); ctx.strokeStyle = "#fff"; ctx.lineWidth = 1; ctx.stroke(rings);
    }
    if (this.box) {  // the analyzed part of a large image
      const [bx, by, bw, bh] = this.box;
      ctx.strokeStyle = "rgba(255,255,255,.7)"; ctx.lineWidth = 1; ctx.setLineDash([6, 4]);
      ctx.strokeRect(tx + bx * s, ty + by * s, bw * s, bh * s); ctx.setLineDash([]);
    }
    if (this.pixelSize) {  // scale bar
      const um = 110 / s * this.pixelSize, e = Math.pow(10, Math.floor(Math.log10(um)));
      const len = [5, 2, 1].map((k) => k * e).find((v) => v <= um) || e, px = len / this.pixelSize * s;
      ctx.fillStyle = "#fff"; ctx.fillRect(14, ch - 18, px, 3);
      ctx.font = "12px system-ui, sans-serif"; ctx.shadowColor = "#000"; ctx.shadowBlur = 3;
      ctx.fillText(`${len >= 1 ? len : len.toPrecision(1)} ${this.unit}`, 14, ch - 24); ctx.shadowBlur = 0;
    }
  }
}
