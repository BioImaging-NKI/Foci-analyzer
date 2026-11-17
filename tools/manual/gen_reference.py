"""The settings reference of the manual (manual/reference/settings.md), made from the program's own labels and
context help (config.py), so it cannot drift from the window. tests/test_manual.py fails when it is stale.

    python tools/manual/gen_reference.py          (with the Python of Foci Analyzer)
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from foci_analyzer import config

OUT = ROOT / "manual" / "reference" / "settings.md"


def _value(v) -> str:
    if isinstance(v, bool):
        return "on" if v else "off"
    if v == "":
        return "(empty)"
    return f"{v:g}" if isinstance(v, float) else str(v)


def _shown_when(f: dict, labels: dict) -> str:
    if not f["when"]:
        return ""
    parts = [f"*{labels.get(k, k)}* is {' or '.join(f'`{v}`' for v in values)}" for k, values in f["when"].items()]
    return " Shown when " + " and ".join(parts) + "."


def render() -> str:
    lines = ["<!-- Made by tools/manual/gen_reference.py from config.py: do not edit by hand. -->", "",
             "# All settings", "",
             ("Every analysis setting of the panel on the left of the **Analyze** page, with its default and the "
             "text of its **?**. Settings marked *advanced* are under **Advanced** in their section. Settings "
             "files (**Save settings…**, and `settings.toml` in every results folder) use the names in "
             "`code`, under the section name in brackets."), ""]
    for sec in config.schema():
        labels = {f["key"]: f["label"] for f in sec["fields"]}
        title = f"{sec['title']} (one per foci channel)" if sec["repeat"] else sec["title"]
        lines += [f"## {title}", "", f"In a settings file: `[{'[' + sec['key'] + ']' if sec['repeat'] else sec['key']}]`.", ""]
        lines += ["| Setting | Default | Description |", "|---|---|---|"]
        for f in sec["fields"]:
            name = f"**{f['label']}**<br>`{f['key']}`" + ("<br>*advanced*" if f["advanced"] else "")
            default = _value(f["default"]) + (f" {f['unit']}" if f["unit"] and f["default"] != "" else "")
            text = f["help"] or ""
            if f["choices"]:
                text += " Choices: " + ", ".join(f"`{c}`" for c in f["choices"]) + "."
            if f["min"] is not None or f["max"] is not None:
                lo = "" if f["min"] is None else _value(f["min"])
                hi = "" if f["max"] is None else _value(f["max"])
                text += f" Range: {lo} to {hi}." if lo and hi else (f" At least {lo}." if lo else f" At most {hi}.")
            if f.get("catalog"):
                text += " The list of models is on [Models](models.md)."
            text += _shown_when(f, labels)
            lines.append(f"| {name} | {default} | {text.strip().replace('|', '/')} |")
        lines.append("")
    return "\n".join(lines)


def render_models() -> str:
    from foci_analyzer import models
    lines = ["<!-- Made by tools/manual/gen_reference.py from models.py: do not edit by hand. -->", "",
             "# Models", "",
             ("The pretrained models that the **Model** lists offer, per network. A model needs its network to be "
             "installed (see [The Settings page](../guide/settings-page.md)); models of a network that is not "
             "installed are greyed out. *Own model…*, at the bottom of every list, takes the folder or file of a "
             "model you trained yourself."), ""]
    names = {"stardist": "StarDist", "cellpose": "Cellpose", "instanseg": "InstanSeg", "spotiflow": "Spotiflow"}
    env = {"cellpose": "Cellpose 3", "cellpose4": "Cellpose 4"}
    for net, items in models.CATALOG.items():
        lines += [f"## {names[net]}", "", f"More about these models: <{models.DOCS[net]}>", "",
                  "| Name | In the list | Runs in |", "|---|---|---|"]
        for name, label, runs in items:
            lines.append(f"| `{name or '(default)'}` | {label} | {env.get(runs, names.get(runs, runs))} |")
        lines.append("")
    return "\n".join(lines)


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(render(), encoding="utf-8")
    (OUT.parent / "models.md").write_text(render_models(), encoding="utf-8")
    print(f"wrote {OUT} and models.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
