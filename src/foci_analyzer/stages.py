"""The analysis as a chain of cached stages (plan.md, section 5).

A stage result is stored under a key made of the stage name, its own settings and the keys of its inputs.
Changing a setting therefore recomputes that stage and what follows, and nothing before it: this is what
makes the preview fast.
"""

from __future__ import annotations

import hashlib
import json
import threading
from collections import OrderedDict
from collections.abc import Callable
from typing import Any


def stage_key(name: str, params: dict, *inputs: str) -> str:
    text = json.dumps([name, params, list(inputs)], sort_keys=True, default=str)
    return hashlib.sha1(text.encode()).hexdigest()[:16]


class StageCache:
    def __init__(self, keep: int = 64):
        self.keep = keep
        self._items: OrderedDict[str, Any] = OrderedDict()
        self._lock = threading.Lock()
        self.computed: list[str] = []  # names of the stages that really ran (for tests and timing)

    def run(self, name: str, params: dict, inputs: tuple[str, ...], compute: Callable[[], Any]) -> tuple[str, Any]:
        """The key and the result of a stage; `compute` is only called when the result is not there yet."""
        key = stage_key(name, params, *inputs)
        with self._lock:
            if key in self._items:
                self._items.move_to_end(key)
                return key, self._items[key]
        value = compute()
        with self._lock:
            self._items[key] = value
            self.computed.append(name)
            while len(self._items) > self.keep:
                self._items.popitem(last=False)
        return key, value
