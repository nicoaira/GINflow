#!/usr/bin/env python3
"""Plot query/target 2Ds from an alignments TSV with RNArtistCore."""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


GRAY = "#adb5bd"
GRAY_LETTER = "#6c757d"

KTS = """rnartist {{
    ss {{
        bn {{
            seq   = "{seq}"
            value = "{dot}"
            name  = "{name}"
        }}
    }}
    theme {{
        details {{ value = 4 }}
{theme}    }}
    svg {{
        path = "{out_dir}"
    }}
}}
"""


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    return cleaned or "rna"


def _color_block(value: str, type_name: str, location: tuple[int, int] | None = None) -> str:
    lines = [
        "        color {",
        f'            type = "{type_name}"',
        f'            value = "{value}"',
    ]
    if location is not None:
        start, end = location
        lines.append(f"            location {{ {start} to {end} }}")
    lines.append("        }")
    return "\n".join(lines)


def theme_block(start: int, end: int, colour: str) -> str:
    """Gray the whole 2D, then paint the aligned span (0-based half-open)."""
    parts = [
        _color_block(GRAY, "N"),
        _color_block(GRAY_LETTER, "n"),
        _color_block(GRAY, "phosphodiester_bond"),
        _color_block(GRAY, "secondary_interaction"),
    ]
    if end > start:
        loc = (start + 1, end)
        parts.extend([
            _color_block(colour, "N", loc),
            _color_block("#212529", "n", loc),
            _color_block(colour, "phosphodiester_bond", loc),
            _color_block(colour, "secondary_interaction", loc),
        ])
    return "\n".join(parts) + "\n"


def draw_one(seq: str, dot: str, name: str, start: int, end: int, colour: str, work: Path) -> Path | None:
    seq = seq.replace("T", "U").replace("t", "u")
    if not seq or len(seq) != len(dot):
        print(f"warning: skip {name}: sequence/structure length mismatch", file=sys.stderr)
        return None
    script = work / f"{name}.kts"
    script.write_text(KTS.format(
        seq=seq, dot=dot, name=name, out_dir=str(work),
        theme=theme_block(start, end, colour),
    ))
    proc = subprocess.run(["rnartistcore", str(script)], cwd=work, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"warning: rnartistcore failed for {name}: {proc.stderr.strip()}", file=sys.stderr)
        return None
    svg = work / f"{name}.svg"
    if not svg.is_file():
        matches = list(work.glob(f"**/{name}.svg"))
        if not matches:
            print(f"warning: no SVG produced for {name}", file=sys.stderr)
            return None
        svg = matches[0]
    return svg


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--highlight-colour", default="#24B064")
    parser.add_argument("--max-pairs", type=int, default=25)
    parser.add_argument("--cpus", type=int, default=1)
    return parser.parse_args(argv)


def pair_stem(row: dict[str, str], index: int) -> str:
    qid = row.get("query_id", f"query_{index}")
    tid = row.get("target_id", f"target_{index}")
    cluster = row.get("cluster_id", str(index))
    return f"{safe_name(cluster)}_{safe_name(qid)}__{safe_name(tid)}"


def collect_jobs(rows: list[dict[str, str]], max_pairs: int, colour: str) -> list[tuple]:
    jobs: list[tuple] = []
    for index, row in enumerate(rows[:max_pairs]):
        prefix = pair_stem(row, index)
        for kind, seq_key, ss_key, start_key, end_key in (
            ("query", "query_sequence", "query_structure", "query_start", "query_end"),
            ("target", "target_sequence", "target_structure", "target_start", "target_end"),
        ):
            jobs.append((
                row.get(seq_key, ""),
                row.get(ss_key, ""),
                f"{prefix}_{kind}",
                int(float(row[start_key])),
                int(float(row[end_key])),
                colour,
            ))
    return jobs


def draw_job(job: tuple) -> tuple[str, bytes] | None:
    seq, dot, name, start, end, colour = job
    with tempfile.TemporaryDirectory(prefix=f"rnartist_{name}_") as tmp:
        svg = draw_one(seq, dot, name, start, end, colour, Path(tmp))
        if svg is None:
            return None
        return svg.name, svg.read_bytes()


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    with args.alignments.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    jobs = collect_jobs(rows, args.max_pairs, args.highlight_colour)
    workers = max(1, args.cpus)
    drawn = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(draw_job, jobs):
            if result is None:
                continue
            name, payload = result
            (args.outdir / name).write_bytes(payload)
            drawn += 1
    print(f"wrote {drawn} SVGs to {args.outdir} ({workers} workers)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
