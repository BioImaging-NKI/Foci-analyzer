"""The analysis as a chain of cached stages (plan.md, section 5).

A stage result is stored under a key made of the stage name, its own settings and the keys of its inputs.
Changing a setting therefore recomputes that stage and what follows, and nothing before it: this is what
makes the preview fast.
"""

from __future__ import annotations

import hashlib
import json
import threading
import time
from collections import OrderedDict
from collections.abc import Callable
from typing import Any


def stage_key(name: str, params: dict, *inputs: str) -> str:
    text = json.dumps([name, params, list(inputs)], sort_keys=True, default=str)
    return hashlib.sha1(text.encode()).hexdigest()[:16]


class Canceled(ValueError):
    """The user stopped the analysis."""


class StageCache:
    def __init__(self, keep: int = 64):
        self.keep = keep
        self._items: OrderedDict[str, Any] = OrderedDict()
        self._lock = threading.Lock()
        self.computed: list[str] = []  # names of the stages that really ran (for tests)
        self.seconds: dict[str, float] = {}  # time spent per stage name
        self.canceled = False  # set to stop the analysis before its next stage
        self.now: tuple | None = None  # the stage being computed: name, note, step, of steps, start time

    def run(self, name: str, params: dict, inputs: tuple[str, ...], compute: Callable[[], Any],
            note: str = "", step: tuple[int, int] = (0, 0)) -> tuple[str, Any]:
        """The key and the result of a stage; `compute` is only called when the result is not there yet.

        `note` (what is being done, for the user) and `step` (number, of how many) are shown in the status
        bar while the stage is computed.
        """
        key = stage_key(name, params, *inputs)
        with self._lock:
            if key in self._items:
                self._items.move_to_end(key)
                return key, self._items[key]
        if self.canceled:
            raise Canceled("stopped")
        start = time.perf_counter()
        self.now = (name, note or name, step[0], step[1], time.time())
        try:
            value = compute()
        finally:
            self.now = None
        self.seconds[name] = self.seconds.get(name, 0.0) + time.perf_counter() - start
        with self._lock:
            self._items[key] = value
            self.computed.append(name)
            while len(self._items) > self.keep:
                self._items.popitem(last=False)
        return key, value
