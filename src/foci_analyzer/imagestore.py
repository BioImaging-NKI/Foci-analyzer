"""Open images for the viewer: planes (a slice or the max projection) at pyramid levels, kept in memory."""

from __future__ import annotations

import threading
from collections import OrderedDict

import numpy as np

from .io import ImageInfo, read_image


def downsample(plane: np.ndarray, level: int) -> np.ndarray:
    """Average 2^level × 2^level pixels; a partial block at the right or bottom edge is averaged as it is."""
    f = 1 << level
    if f == 1:
        return plane
    h, w = plane.shape
    H, W = -(-h // f), -(-w // f)
    padded = np.full((H * f, W * f), np.nan, dtype=np.float32)
    padded[:h, :w] = plane
    return np.nanmean(padded.reshape(H, f, W, f), axis=(1, 3))


def to_uint16(plane: np.ndarray) -> np.ndarray:
    if plane.dtype == np.uint16:
        return plane
    return np.clip(np.rint(plane), 0, 65535).astype(np.uint16)


class OpenImage:
    def __init__(self, data: np.ndarray, info: ImageInfo):
        self.data, self.info = data, info
        self._planes: dict[tuple, np.ndarray] = {}
        self._lock = threading.Lock()

    def plane(self, t: int, c: int, z: int | None, level: int = 0) -> np.ndarray:
        """Channel c (0-based) of frame t as 16-bit: slice z, or the max projection when z is None."""
        key = (t, c, z, level)
        with self._lock:
            if key not in self._planes:
                stack = self.data[t, c]
                full = stack.max(axis=0) if z is None else stack[z]
                self._planes[key] = to_uint16(downsample(full, level))
            return self._planes[key]

    def display_range(self, c: int) -> tuple[int, int]:
        """A contrast range for channel c: 0.1 and 99.9 percentile of the max projection of the first frame."""
        lo, hi = np.percentile(self.plane(0, c, None, 0), [0.1, 99.9])
        return int(lo), int(max(hi, lo + 1))


class ImageStore:
    """The last few opened images (reading a file again is slower than keeping it)."""

    def __init__(self, keep: int = 3):
        self.keep = keep
        self._open: OrderedDict[str, OpenImage] = OrderedDict()
        self._lock = threading.Lock()

    def get(self, path: str) -> OpenImage:
        with self._lock:
            if path in self._open:
                self._open.move_to_end(path)
                return self._open[path]
        image = OpenImage(*read_image(path))
        with self._lock:
            self._open[path] = image
            while len(self._open) > self.keep:
                self._open.popitem(last=False)
        return image
