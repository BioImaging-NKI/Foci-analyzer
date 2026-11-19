"""The figure of the multi-scale section of the manual (manual/images/multiscale_compare.png): watershed at a small
and at a large foci size against multi-scale, on two crops of the max projection of a 3D image with foci of very
different sizes. The crops are chosen automatically: one where the small size splits a large focus, one where the
large size misses small foci.

    <Foci Analyzer python> tools/manual/multiscale_figure.py IMAGE [--nuclei 3] [--channel 2]

The image used: TestData/Hi-res 3D/ZT_CO_8Y_AP_63x_1.tif (not in git; the user allowed its use in the manual).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from foci_analyzer import config, io, pipeline
from foci_analyzer.stages import StageCache

OUT = Path(__file__).resolve().parents[2] / "manual" / "images" / "multiscale_compare.png"
SMALL, LARGE, RANGE = 0.3, 1.0, (0.3, 2.5)
H, ZOOM = 300, 2  # crop size (px of the image) and enlargement


def analyze(data, info, nuclei, channel, **foci):
    s = config.Settings()
    s.image.nuclei_channel = nuclei
    s.segmentation.method, s.segmentation.exclude_on_edges = "threshold", False
    s.segmentation.min_diameter, s.segmentation.max_diameter = 6, 40
    s.foci_general.dimension = "max projection"
    s.foci = [config.FociChannelSettings(channel=channel, **foci)]
    f = pipeline.analyze_frame(StageCache(), "figure", data[0], info, s).foci
    return f.x.to_numpy() / info.pixel_size, f.y.to_numpy() / info.pixel_size, f


def crops(ws_small, ws_large, ms, shape, px):
    """Top-left corners of two crops that do not overlap: the most splits of large foci by the small size, and
    the most small foci that the large size misses (both counted against multi-scale)."""
    mx, my, mf = ms
    size = mf.detected_size.to_numpy()
    big = size >= 1.0
    # per multi-scale focus: how many extra markers the small size puts inside a large one, and whether the
    # large size misses a small one
    split = np.array([max(0, int((np.hypot(ws_small[0] - x, ws_small[1] - y) < d / px / 2).sum()) - 1) if b else 0
                      for x, y, d, b in zip(mx, my, size, big, strict=True)])
    missed = np.array([(not b) and np.min(np.hypot(ws_large[0] - x, ws_large[1] - y), initial=99) > 4
                       for x, y, b in zip(mx, my, big, strict=True)], int)
    windows = []
    for y0 in range(0, shape[0] - H, H // 4):
        for x0 in range(0, shape[1] - H, H // 4):
            inside = (mx >= x0 + 20) & (mx < x0 + H - 20) & (my >= y0 + 20) & (my < y0 + H - 20)
            windows.append((int(split[inside].sum()), int(missed[inside].sum()), (y0, x0)))
    a = max(windows, key=lambda w: (w[0], w[1]))[2]
    apart = [w for w in windows if abs(w[2][0] - a[0]) >= H or abs(w[2][1] - a[1]) >= H]
    b = max(apart, key=lambda w: (w[1], w[0]))[2]
    return [a, b]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("image")
    p.add_argument("--nuclei", type=int, default=3)
    p.add_argument("--channel", type=int, default=2)
    a = p.parse_args()
    data, info = io.read_image(a.image)
    px = info.pixel_size
    ws_small = analyze(data, info, a.nuclei, a.channel, size=SMALL)
    ws_large = analyze(data, info, a.nuclei, a.channel, size=LARGE)
    ms = analyze(data, info, a.nuclei, a.channel, method="multi-scale", size=RANGE[0], size_max=RANGE[1])
    mip = data[0, a.channel - 1].max(0).astype(np.float32)
    corners = crops(ws_small, ws_large, ms, mip.shape, px)
    names = [f"watershed, foci size {SMALL} µm", f"watershed, foci size {LARGE} µm",
             f"multi-scale, {RANGE[0]}–{RANGE[1]} µm"]
    try:  # a font with µ and the en dash (the built-in one has neither)
        font = ImageFont.truetype("segoeui.ttf", 22)
    except OSError:
        font = ImageFont.truetype("DejaVuSans.ttf", 20)
    rows = []
    for y0, x0 in corners:
        crop = mip[y0:y0 + H, x0:x0 + H]
        lo, hi = np.percentile(crop, [1, 99.8])
        gray = (np.clip((crop - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
        panels = []
        for name, (xs, ys, f) in zip(names, (ws_small, ws_large, ms), strict=True):
            im = Image.fromarray(gray).convert("RGB").resize((ZOOM * H, ZOOM * H), Image.NEAREST)
            d = ImageDraw.Draw(im)
            sizes = f.detected_size.to_numpy() if "detected_size" in f else None
            for i, (x, y) in enumerate(zip(xs - x0, ys - y0, strict=True)):
                if 0 <= x < H and 0 <= y < H:
                    cx, cy = ZOOM * x, ZOOM * y
                    if sizes is not None:  # a circle as wide as the size at which the focus was found
                        r = max(4, ZOOM * sizes[i] / px / 2)
                        d.ellipse([cx - r, cy - r, cx + r, cy + r], outline=(255, 70, 190), width=3)
                    else:
                        d.ellipse([cx - 5, cy - 5, cx + 5, cy + 5], fill=(255, 205, 50), outline=(0, 0, 0))
            d.rectangle([0, 0, ZOOM * H, 38], fill=(0, 0, 0))
            d.text((10, 3), name, fill=(255, 255, 255), font=font)
            bar = ZOOM * 5 / px  # 5 µm
            d.line([14, ZOOM * H - 18, 14 + bar, ZOOM * H - 18], fill=(255, 255, 255), width=5)
            panels.append(np.pad(np.asarray(im), ((4, 4), (4, 4), (0, 0)), constant_values=255))
        rows.append(np.concatenate(panels, 1))
    Image.fromarray(np.concatenate(rows, 0)).save(OUT, optimize=True)
    print(OUT, "crops at", corners)


if __name__ == "__main__":
    main()
