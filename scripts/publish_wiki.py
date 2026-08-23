#!/usr/bin/env python3
"""Publish docs/*.md to the GitHub wiki (nicoaira/ginflow.wiki.git).

The wiki git repo does not exist until someone creates the first page in
the GitHub UI. After that, this script can clone, rewrite, and push.
"""
from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
WIKI_REMOTE = "https://github.com/nicoaira/ginflow.wiki.git"

MAPPING = {
    "README.md": "Home.md",
    "getting-started.md": "Getting-started.md",
    "how-it-works.md": "How-it-works.md",
    "usage.md": "Usage.md",
    "indexes.md": "Indexes.md",
    "clustering-and-alignment.md": "Clustering-and-alignment.md",
    "statistics.md": "E-values.md",
    "output.md": "Output.md",
    "plotting.md": "Plotting-and-report.md",
    "parameters.md": "Parameters.md",
    "profiles.md": "Profiles-and-hardware.md",
    "development.md": "Development.md",
    "faq.md": "FAQ.md",
}

LINK_MAP = [
    ("clustering-and-alignment.md", "Clustering-and-alignment"),
    ("getting-started.md", "Getting-started"),
    ("how-it-works.md", "How-it-works"),
    ("parameters.md", "Parameters"),
    ("development.md", "Development"),
    ("statistics.md", "E-values"),
    ("plotting.md", "Plotting-and-report"),
    ("profiles.md", "Profiles-and-hardware"),
    ("indexes.md", "Indexes"),
    ("output.md", "Output"),
    ("usage.md", "Usage"),
    ("faq.md", "FAQ"),
    ("README.md", "Home"),
]

RAW_IMG = "https://raw.githubusercontent.com/nicoaira/ginflow/main/docs/"
BLOB = "https://github.com/nicoaira/ginflow/blob/main/"

SIDEBAR = """**GINflow**

* [Home](Home)
* [Getting started](Getting-started)
* [How it works](How-it-works)
* [Usage](Usage)

**Using the pipeline**

* [Window indexes](Indexes)
* [Clustering and alignment](Clustering-and-alignment)
* [E-values](E-values)
* [Output](Output)
* [Plotting and report](Plotting-and-report)
* [Parameters](Parameters)
* [Profiles and hardware](Profiles-and-hardware)

**Reference**

* [Development](Development)
* [FAQ](FAQ)
"""

FOOTER = (
    "GINflow · [source repository](https://github.com/nicoaira/ginflow) · "
    "[docs/ in git](https://github.com/nicoaira/ginflow/tree/main/docs)\n"
)


def rewrite(text: str) -> str:
    text = re.sub(
        r"\]\((images/[^)]+)\)",
        lambda match: f"]({RAW_IMG}{match.group(1)})",
        text,
    )
    text = re.sub(
        r'src="(images/[^"]+)"',
        lambda match: f'src="{RAW_IMG}{match.group(1)}"',
        text,
    )
    text = re.sub(
        r"\]\(\.\./([^)]+)\)",
        lambda match: f"]({BLOB}{match.group(1)})",
        text,
    )
    for fname, wiki in LINK_MAP:
        text = text.replace(f"]({fname})", f"]({wiki})")
        text = text.replace(f"]({fname}#", f"]({wiki}#")
    return text


def write_pages(dest: Path) -> None:
    for src_name, wiki_name in MAPPING.items():
        (dest / wiki_name).write_text(rewrite((DOCS / src_name).read_text()))
    (dest / "_Sidebar.md").write_text(SIDEBAR)
    (dest / "_Footer.md").write_text(FOOTER)


def git(cwd: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["git", *args],
        cwd=cwd,
        check=True,
        text=True,
        capture_output=True,
    )


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="ginflow-wiki-") as tmp:
        dest = Path(tmp) / "wiki"
        clone = subprocess.run(
            ["git", "clone", WIKI_REMOTE, str(dest)],
            text=True,
            capture_output=True,
        )
        if clone.returncode != 0:
            dest.mkdir()
            git(dest, "init", "-b", "master")
            git(dest, "remote", "add", "origin", WIKI_REMOTE)
            print(
                "Wiki remote is not initialized yet.\n"
                "Open https://github.com/nicoaira/ginflow/wiki and click "
                "'Create the first page', then re-run this script.",
                file=sys.stderr,
            )
            print(clone.stderr, file=sys.stderr)
            return 1

        for path in dest.glob("*.md"):
            path.unlink()
        write_pages(dest)
        git(dest, "add", "-A")
        status = git(dest, "status", "--porcelain")
        if not status.stdout.strip():
            print("Wiki already up to date.")
            return 0
        git(dest, "commit", "-m", "Sync wiki from docs/")
        git(dest, "push", "origin", "HEAD")
        print("Pushed", WIKI_REMOTE)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(exc.stderr or str(exc))
        sys.exit(exc.returncode)
