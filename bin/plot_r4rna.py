#!/usr/bin/env python3
"""Plot a query/target alignment as one mirrored arc diagram."""
from __future__ import annotations

import argparse
import csv
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


INK = "#212529"
MUTE = "#6c757d"
TRACK = "#e9ecef"
QUERY_ONLY = "#396E35"
TARGET_ONLY = "#3F2B29"
MISMATCH = "#ECDC86"
GAP = "#adb5bd"
RESIDUE = "#212529"
PAIR = "#24B064"

OPENS = {"(": ")", "<": ">", "[": "]", "{": "}", "A": "a", "B": "b", "C": "c", "D": "d"}
CLOSES = {end: start for start, end in OPENS.items()}

MARGIN_LEFT = 18
MARGIN_RIGHT = 18
MARGIN_TOP = 34
MARGIN_BOTTOM = 36
ARC_MAX = 118.0
AXIS_GAP = 24.0
BACKBONE_RESIDUE = 2.6
BACKBONE_GAP = 1.15
RIBBON = 8.4


def safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_")
    return cleaned or "rna"


def xml_escape(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def parse_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def normalize_base(char: str) -> str:
    char = char.upper()
    return "U" if char == "T" else char


def vienna_pairs(dot: str) -> list[tuple[int, int]]:
    stacks: dict[str, list[int]] = {close: [] for close in CLOSES}
    pairs: list[tuple[int, int]] = []
    for index, char in enumerate(dot):
        if char in OPENS:
            stacks[OPENS[char]].append(index)
        elif char in CLOSES:
            stack = stacks[char]
            if stack:
                pairs.append((stack.pop(), index))
    return pairs


def position_map(aligned: str, start: int) -> dict[int, int]:
    mapping: dict[int, int] = {}
    pos = start
    for column, char in enumerate(aligned):
        if char != "-":
            mapping[pos] = column
            pos += 1
    return mapping


def aligned_pairs(dot: str, aligned: str, start: int) -> list[tuple[int, int]]:
    mapping = position_map(aligned, start)
    pairs = []
    for left, right in vienna_pairs(dot):
        if left in mapping and right in mapping:
            column_i = mapping[left]
            column_j = mapping[right]
            pairs.append((min(column_i, column_j), max(column_i, column_j)))
    return pairs


def identity_states(query_aln: str, target_aln: str) -> list[str]:
    states = []
    for query, target in zip(query_aln, target_aln):
        if query == "-" or target == "-":
            states.append("gap")
        elif normalize_base(query) == normalize_base(target):
            states.append("match")
        else:
            states.append("mismatch")
    return states


def runs(values: list[str]) -> list[tuple[int, int, str]]:
    out: list[tuple[int, int, str]] = []
    index = 0
    while index < len(values):
        stop = index + 1
        while stop < len(values) and values[stop] == values[index]:
            stop += 1
        out.append((index, stop, values[index]))
        index = stop
    return out


def resolved_alignment(row: dict[str, str]) -> tuple[str, str] | None:
    query_aln = row.get("query_aligned", "").strip()
    target_aln = row.get("target_aligned", "").strip()
    if query_aln and target_aln:
        if len(query_aln) != len(target_aln):
            return None
        return query_aln.replace("T", "U").replace("t", "u"), target_aln.replace("T", "U").replace("t", "u")
    query_seq = row.get("query_sequence", "")
    target_seq = row.get("target_sequence", "")
    query_start = parse_int(row.get("query_start", ""))
    query_end = parse_int(row.get("query_end", ""))
    target_start = parse_int(row.get("target_start", ""))
    target_end = parse_int(row.get("target_end", ""))
    query_slice = query_seq[query_start:query_end]
    target_slice = target_seq[target_start:target_end]
    if query_slice and len(query_slice) == len(target_slice):
        return query_slice.replace("T", "U").replace("t", "u"), target_slice.replace("T", "U").replace("t", "u")
    return None


def column_width(count: int) -> float:
    if count <= 0:
        return 8.0
    return max(3.2, min(8.5, 1000.0 / count))


def pair_key(pair: tuple[int, int]) -> str:
    return f"{pair[0]}:{pair[1]}"


def colour_pairs(
    query_pairs: list[tuple[int, int]],
    target_pairs: list[tuple[int, int]],
    highlight: str,
) -> tuple[list[tuple[int, int, str]], list[tuple[int, int, str]]]:
    shared = {pair_key(pair) for pair in query_pairs} & {pair_key(pair) for pair in target_pairs}
    query_coloured = [
        (left, right, highlight if pair_key((left, right)) in shared else QUERY_ONLY)
        for left, right in query_pairs
    ]
    target_coloured = [
        (left, right, highlight if pair_key((left, right)) in shared else TARGET_ONLY)
        for left, right in target_pairs
    ]
    return query_coloured, target_coloured


def arc_path(x1: float, x2: float, y: float, flip: bool) -> str:
    radius_x = abs(x2 - x1) / 2.0
    if radius_x <= 0:
        return ""
    radius_y = min(radius_x, ARC_MAX)
    sweep = 0 if flip else 1
    start, end = (x1, x2) if x1 <= x2 else (x2, x1)
    return (
        f"M {start:.2f},{y:.2f} "
        f"A {radius_x:.2f},{radius_y:.2f} 0 0 {sweep} {end:.2f},{y:.2f}"
    )


def swatch(x: float, y: float, colour: str, kind: str) -> str:
    if kind == "line":
        return (
            f'<line x1="{x:.1f}" y1="{y:.1f}" x2="{x + 12:.1f}" y2="{y:.1f}" '
            f'stroke="{colour}" stroke-width="1.8" stroke-linecap="round"/>'
        )
    return (
        f'<rect x="{x:.1f}" y="{y - 4:.1f}" width="12" height="8" rx="1.5" fill="{colour}"/>'
    )


def render_svg(
    *,
    query_id: str,
    target_id: str,
    query_aln: str,
    target_aln: str,
    query_pairs: list[tuple[int, int]],
    target_pairs: list[tuple[int, int]],
    query_start: int,
    query_end: int,
    target_start: int,
    target_end: int,
    highlight: str,
    marker_id: str,
) -> str:
    count = len(query_aln)
    col_w = column_width(count)
    plot_left = MARGIN_LEFT
    plot_width = count * col_w
    width = max(640.0, plot_left + plot_width + MARGIN_RIGHT)
    mid_y = MARGIN_TOP + ARC_MAX + AXIS_GAP / 2.0
    y_query = mid_y - AXIS_GAP / 2.0
    y_target = mid_y + AXIS_GAP / 2.0
    height = MARGIN_TOP + ARC_MAX + AXIS_GAP + ARC_MAX + MARGIN_BOTTOM
    match_col = highlight or PAIR

    def col_center(index: int) -> float:
        return plot_left + (index + 0.5) * col_w

    def col_edge(index: int) -> float:
        return plot_left + index * col_w

    query_coloured, target_coloured = colour_pairs(query_pairs, target_pairs, match_col)
    query_coloured.sort(key=lambda item: item[1] - item[0], reverse=True)
    target_coloured.sort(key=lambda item: item[1] - item[0], reverse=True)

    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.1f} {height:.1f}" '
        f'width="{width:.1f}" height="{height:.1f}" '
        f'font-family="ui-sans-serif, system-ui, sans-serif" font-size="11" fill="{INK}" '
        f'role="img" aria-label="Aligned R4RNA arc plot of {xml_escape(query_id)} versus {xml_escape(target_id)}">',
        f"<title>{xml_escape(query_id)} vs {xml_escape(target_id)}</title>",
        f'<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<defs><marker id="{xml_escape(marker_id)}" viewBox="0 0 10 10" refX="8" refY="5" '
        f'markerWidth="5.5" markerHeight="5.5" orient="auto">'
        f'<path d="M 0 1.2 L 10 5 L 0 8.8 z" fill="{RESIDUE}"/></marker></defs>',
    ]

    def draw_arcs(items: list[tuple[int, int, str]], y: float, flip: bool) -> None:
        for left, right, colour in items:
            path = arc_path(col_center(left), col_center(right), y, flip)
            if not path:
                continue
            parts.append(
                f'<path d="{path}" fill="none" stroke="{colour}" stroke-width="1.35" '
                f'stroke-linecap="round" stroke-linejoin="round"/>'
            )

    draw_arcs(query_coloured, y_query, flip=False)
    draw_arcs(target_coloured, y_target, flip=True)

    track_x = col_edge(0)
    track_w = plot_width
    track_y = mid_y - RIBBON / 2.0 - 2.2
    parts.append(
        f'<rect x="{track_x:.2f}" y="{track_y:.2f}" width="{track_w:.2f}" '
        f'height="{RIBBON + 4.4:.2f}" rx="3.2" fill="{TRACK}"/>'
    )

    states = identity_states(query_aln, target_aln)
    colours = {"match": match_col, "mismatch": MISMATCH, "gap": GAP}
    for start, stop, state in runs(states):
        parts.append(
            f'<rect x="{col_edge(start):.2f}" y="{mid_y - RIBBON / 2.0:.2f}" '
            f'width="{(stop - start) * col_w:.2f}" height="{RIBBON:.2f}" '
            f'rx="1.2" fill="{colours[state]}"/>'
        )

    def draw_backbone(aligned: str, y: float) -> None:
        chars = list(aligned)
        flags = ["gap" if char == "-" else "residue" for char in chars]
        for start, stop, kind in runs(flags):
            parts.append(
                f'<line x1="{col_edge(start):.2f}" y1="{y:.2f}" x2="{col_edge(stop):.2f}" y2="{y:.2f}" '
                f'stroke="{GAP if kind == "gap" else RESIDUE}" '
                f'stroke-width="{BACKBONE_GAP if kind == "gap" else BACKBONE_RESIDUE}" '
                f'stroke-linecap="butt"/>'
            )
        tip = col_edge(len(chars)) + 1.2
        parts.append(
            f'<line x1="{col_edge(len(chars)) - 0.2:.2f}" y1="{y:.2f}" x2="{tip:.2f}" y2="{y:.2f}" '
            f'stroke="{RESIDUE}" stroke-width="{BACKBONE_RESIDUE}" '
            f'stroke-linecap="round" marker-end="url(#{xml_escape(marker_id)})"/>'
        )

    draw_backbone(query_aln, y_query)
    draw_backbone(target_aln, y_target)

    parts.append(
        f'<text x="{plot_left:.1f}" y="16" font-size="12" font-weight="600">{xml_escape(query_id)}</text>'
    )
    parts.append(
        f'<text x="{plot_left:.1f}" y="28" fill="{MUTE}" font-size="10">'
        f"query {query_start}–{query_end}</text>"
    )
    parts.append(
        f'<text x="{plot_left:.1f}" y="{height - 16:.1f}" font-size="12" font-weight="600">'
        f"{xml_escape(target_id)}</text>"
    )
    parts.append(
        f'<text x="{plot_left:.1f}" y="{height - 4:.1f}" fill="{MUTE}" font-size="10">'
        f"target {target_start}–{target_end}</text>"
    )

    legend = (
        ("shared pair", match_col, "line"),
        ("query only", QUERY_ONLY, "line"),
        ("target only", TARGET_ONLY, "line"),
        ("same base", match_col, "box"),
        ("different", MISMATCH, "box"),
        ("gap", GAP, "box"),
    )
    cursor = width - MARGIN_RIGHT
    for label, colour, kind in reversed(legend):
        label_w = 6.1 * len(label)
        cursor -= label_w + 22
        parts.append(swatch(cursor, 14, colour, kind))
        parts.append(
            f'<text x="{cursor + 16:.1f}" y="17.5" fill="{MUTE}" font-size="10">{label}</text>'
        )

    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def pair_stem(row: dict[str, str], index: int) -> str:
    query_id = row.get("query_id", f"query_{index}")
    target_id = row.get("target_id", f"target_{index}")
    cluster = row.get("cluster_id", str(index))
    return f"{safe_name(cluster)}_{safe_name(query_id)}__{safe_name(target_id)}"


def draw_row(row: dict[str, str], index: int, highlight: str) -> tuple[str, str] | None:
    aligned = resolved_alignment(row)
    query_id = row.get("query_id", f"query_{index}")
    target_id = row.get("target_id", f"target_{index}")
    if aligned is None:
        print(f"warning: skip {query_id} vs {target_id}: missing gapped alignment", file=sys.stderr)
        return None
    query_aln, target_aln = aligned
    query_start = parse_int(row.get("query_start", ""))
    query_end = parse_int(row.get("query_end", ""))
    target_start = parse_int(row.get("target_start", ""))
    target_end = parse_int(row.get("target_end", ""))
    query_pairs = aligned_pairs(row.get("query_structure", ""), query_aln, query_start)
    target_pairs = aligned_pairs(row.get("target_structure", ""), target_aln, target_start)
    stem = pair_stem(row, index)
    svg = render_svg(
        query_id=query_id,
        target_id=target_id,
        query_aln=query_aln,
        target_aln=target_aln,
        query_pairs=query_pairs,
        target_pairs=target_pairs,
        query_start=query_start,
        query_end=query_end,
        target_start=target_start,
        target_end=target_end,
        highlight=highlight,
        marker_id=f"r4_{stem}_arrow",
    )
    return f"{stem}_alignment.svg", svg


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--highlight-colour", default="#24B064")
    parser.add_argument("--max-pairs", type=int, default=25)
    parser.add_argument("--cpus", type=int, default=1)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    with args.alignments.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    jobs = [(row, index, args.highlight_colour) for index, row in enumerate(rows[: args.max_pairs])]
    workers = max(1, min(args.cpus, len(jobs) or 1))
    drawn = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(lambda item: draw_row(*item), jobs):
            if result is None:
                continue
            name, payload = result
            (args.outdir / name).write_text(payload)
            drawn += 1
    print(f"wrote {drawn} SVGs to {args.outdir} ({workers} workers)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
