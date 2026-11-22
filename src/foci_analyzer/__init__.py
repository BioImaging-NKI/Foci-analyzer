"""Foci Analyzer: foci in nuclei and cells, with a browser window."""

__version__ = "2.0.0b4"


def build() -> str:
    """When and from which commit this copy was built (the setup writes BUILD.txt); a development copy says
    its commit. Empty when unknown."""
    from pathlib import Path

    here = Path(__file__).parent
    try:
        return (here / "BUILD.txt").read_text(encoding="utf-8").strip()
    except OSError:
        pass
    head = here.parents[1] / ".git" / "HEAD"
    try:
        ref = head.read_text().strip()
        sha = (head.parent / ref[5:]).read_text().strip() if ref.startswith("ref: ") else ref
        return f"development copy, commit {sha[:7]}"
    except OSError:
        return ""
