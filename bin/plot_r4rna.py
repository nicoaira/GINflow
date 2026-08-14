#!/usr/bin/env python3
"""Plot query/target arc diagrams from an alignments TSV with R4RNA."""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


R_SCRIPT = r'''
library(R4RNA)
seq <- "{seq}"
ss <- "{ss}"
len <- nchar(seq)
lo <- as.integer("{lo}")
hi <- as.integer("{hi}")
hit_col <- "{colour}"
gray_col <- "{gray}"
svg("{svg}", width = max(4, len / 50), height = 4)
helix <- viennaToHelix(ss)
helix$col <- gray_col
if (hi >= lo) {{
    both <- which(helix$i >= lo & helix$i <= hi & helix$j >= lo & helix$j <= hi)
    helix$col[both] <- hit_col
}}
plotHelix(helix, line = FALSE, arrow = TRUE, lwd = 1.5)
lines(c(0.5, len + 0.5), c(0, 0), col = gray_col, lwd = 3)
if (hi >= lo) {{
    lines(c(lo - 0.5, hi + 0.5), c(0, 0), col = hit_col, lwd = 5)
}}
dev.off()
'''


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    return cleaned or "rna"


def r_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


GRAY = "#B0B0B0"


def draw_one(seq: str, dot: str, name: str, start: int, end: int, colour: str, work: Path) -> Path | None:
    seq = seq.replace("T", "U").replace("t", "u")
    if not seq or len(seq) != len(dot):
        print(f"warning: skip {name}: sequence/structure length mismatch", file=sys.stderr)
        return None
    svg = work / f"{name}.svg"
    script = work / f"{name}.R"
    lo = start + 1 if end > start else 0
    hi = end if end > start else -1
    script.write_text(R_SCRIPT.format(
        seq=r_escape(seq),
        ss=r_escape(dot),
        svg=r_escape(str(svg)),
        lo=lo,
        hi=hi,
        colour=r_escape(colour),
        gray=GRAY,
    ))
    proc = subprocess.run(["Rscript", str(script)], cwd=work, capture_output=True, text=True)
    if proc.returncode != 0 or not svg.is_file():
        print(f"warning: R4RNA failed for {name}: {proc.stderr.strip()}", file=sys.stderr)
        return None
    return svg


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--highlight-colour", default="#CC0000")
    parser.add_argument("--max-plots", type=int, default=50)
    parser.add_argument("--cpus", type=int, default=1)
    return parser.parse_args(argv)


def collect_jobs(rows: list[dict[str, str]], max_plots: int, colour: str) -> list[tuple]:
    jobs: list[tuple] = []
    for index, row in enumerate(rows):
        if len(jobs) >= max_plots:
            break
        qid = row.get("query_id", f"query_{index}")
        tid = row.get("target_id", f"target_{index}")
        prefix = f"{index:04d}_{safe_name(qid)}__{safe_name(tid)}"
        for kind, seq_key, ss_key, start_key, end_key in (
            ("query", "query_sequence", "query_structure", "query_start", "query_end"),
            ("target", "target_sequence", "target_structure", "target_start", "target_end"),
        ):
            if len(jobs) >= max_plots:
                break
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
    with tempfile.TemporaryDirectory(prefix=f"r4rna_{name}_") as tmp:
        svg = draw_one(seq, dot, name, start, end, colour, Path(tmp))
        if svg is None:
            return None
        return svg.name, svg.read_bytes()


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    with args.alignments.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    jobs = collect_jobs(rows, args.max_plots, args.highlight_colour)
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
