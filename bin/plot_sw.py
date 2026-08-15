#!/usr/bin/env python3
"""Plot crop cosine and substitution-score matrices with the SW traceback."""
from __future__ import annotations

import argparse
import base64
import csv
import json
import re
import struct
import sys
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from ginfinity_sw import ScoringParameters, align, similarity_matrix, transform_scores

SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")
MAX_HEATMAP = 480
PANEL = 420
MARGIN_LEFT = 58
MARGIN_RIGHT = 78
MARGIN_TOP = 42
MARGIN_BOTTOM = 52
VIRIDIS = (
    (0.00, (68, 1, 84)),
    (0.25, (59, 82, 139)),
    (0.50, (33, 145, 140)),
    (0.75, (94, 201, 98)),
    (1.00, (253, 231, 37)),
)
RDBU = (
    (0.00, (5, 48, 97)),
    (0.25, (67, 147, 195)),
    (0.50, (247, 247, 247)),
    (0.75, (214, 96, 77)),
    (1.00, (103, 0, 31)),
)


def safe_name(value: str) -> str:
    cleaned = SAFE_NAME.sub("_", value).strip("_")
    return cleaned or "rna"


def xml_escape(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def load_parameters(path: Path) -> ScoringParameters:
    payload = load_json(path)
    return ScoringParameters(**payload.get("scoring_parameters", payload))


def load_clusters(path: Path) -> dict[str, dict]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    clusters = {}
    for row in rows:
        clusters[row["cluster_id"]] = row
    return clusters


def load_npz_shards(paths: list[Path]) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    for path in paths:
        with np.load(path) as archive:
            for key in archive.files:
                if key in arrays:
                    raise ValueError(f"duplicate embedding id {key!r} in {path}")
                arrays[key] = np.asarray(archive[key])
    return arrays


def resolve_targets(database: Path) -> dict[str, np.ndarray]:
    embeddings = database / "embeddings.npz"
    if not embeddings.is_file():
        raise ValueError(
            f"{database} is missing embeddings.npz; "
            "rebuild the database with this pipeline version"
        )
    return load_npz_shards([embeddings])


def crop_bounds(start: int, end: int, length: int, pad: int) -> tuple[int, int]:
    return max(0, start - pad), min(length, end + pad)


def colormap_lut(stops: tuple, size: int = 256) -> np.ndarray:
    lut = np.empty((size, 3), dtype=np.uint8)
    xs = np.linspace(0.0, 1.0, size)
    knots = np.array([stop[0] for stop in stops], dtype=np.float64)
    colors = np.array([stop[1] for stop in stops], dtype=np.float64)
    for channel in range(3):
        lut[:, channel] = np.clip(np.interp(xs, knots, colors[:, channel]), 0, 255)
    return lut


def colorize(matrix: np.ndarray, vmin: float, vmax: float, stops: tuple) -> np.ndarray:
    span = vmax - vmin
    if not np.isfinite(span) or span <= 0:
        norm = np.full(matrix.shape, 0.5, dtype=np.float64)
    else:
        norm = np.clip((matrix - vmin) / span, 0.0, 1.0)
    lut = colormap_lut(stops)
    index = np.rint(norm * (len(lut) - 1)).astype(np.intp)
    return lut[index]


def encode_png(rgb: np.ndarray) -> bytes:
    height, width, _ = rgb.shape
    raw = b"".join(b"\x00" + rgb[row].tobytes() for row in range(height))

    def chunk(tag: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(tag + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)

    return (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )


def downsample(matrix: np.ndarray) -> tuple[np.ndarray, float, float]:
    rows, cols = matrix.shape
    scale_i = max(1, int(np.ceil(rows / MAX_HEATMAP)))
    scale_j = max(1, int(np.ceil(cols / MAX_HEATMAP)))
    if scale_i == 1 and scale_j == 1:
        return matrix, 1.0, 1.0
    new_i = (rows + scale_i - 1) // scale_i
    new_j = (cols + scale_j - 1) // scale_j
    padded = np.full((new_i * scale_i, new_j * scale_j), np.nan, dtype=np.float64)
    padded[:rows, :cols] = matrix
    blocks = padded.reshape(new_i, scale_i, new_j, scale_j)
    reduced = np.nanmean(blocks, axis=(1, 3))
    return reduced, scale_i, scale_j


def ticks(start: int, length: int, count: int = 5) -> list[tuple[float, str]]:
    if length <= 1:
        return [(0.5, str(start))]
    values = np.linspace(0, length - 1, min(count, length))
    marks = []
    seen = set()
    for value in values:
        index = int(round(value))
        if index in seen:
            continue
        seen.add(index)
        marks.append((index + 0.5, str(start + index)))
    return marks


def path_points(
    columns: tuple[tuple[int, int], ...],
    scale_i: float,
    scale_j: float,
) -> list[tuple[float, float]]:
    points = []
    last_q = 0.0
    last_t = 0.0
    for query, target in columns:
        if query >= 0:
            last_q = (query + 0.5) / scale_i
        if target >= 0:
            last_t = (target + 0.5) / scale_j
        if query < 0 and target < 0:
            continue
        points.append((last_t, last_q))
    return points


def colorbar_png(stops: tuple, height: int = 256) -> bytes:
    column = np.linspace(1.0, 0.0, height, dtype=np.float64).reshape(height, 1)
    return encode_png(colorize(column, 0.0, 1.0, stops))


def matrix_svg(
    matrix: np.ndarray,
    *,
    title: str,
    subtitle: str,
    query_offset: int,
    target_offset: int,
    vmin: float,
    vmax: float,
    stops: tuple,
    columns: tuple[tuple[int, int], ...] | None,
    colour: str,
    cbar_low: str,
    cbar_high: str,
    marker_id: str,
) -> str:
    reduced, scale_i, scale_j = downsample(matrix)
    rows, cols = reduced.shape
    rgb = colorize(reduced, vmin, vmax, stops)
    payload = base64.standard_b64encode(encode_png(rgb)).decode("ascii")
    cbar = base64.standard_b64encode(colorbar_png(stops)).decode("ascii")
    width = MARGIN_LEFT + PANEL + MARGIN_RIGHT
    height = MARGIN_TOP + PANEL + MARGIN_BOTTOM
    x_ticks = ticks(target_offset, matrix.shape[1])
    y_ticks = ticks(query_offset, matrix.shape[0])
    tick_svg = []
    for pos, label in x_ticks:
        x = MARGIN_LEFT + PANEL * pos / cols
        y = MARGIN_TOP + PANEL
        tick_svg.append(
            f'<line x1="{x:.1f}" y1="{y}" x2="{x:.1f}" y2="{y + 5}" stroke="#6c757d"/>'
            f'<text x="{x:.1f}" y="{y + 18}" text-anchor="middle">{xml_escape(label)}</text>'
        )
    for pos, label in y_ticks:
        x = MARGIN_LEFT
        y = MARGIN_TOP + PANEL * pos / rows
        tick_svg.append(
            f'<line x1="{x - 5}" y1="{y:.1f}" x2="{x}" y2="{y:.1f}" stroke="#6c757d"/>'
            f'<text x="{x - 8}" y="{y + 3:.1f}" text-anchor="end">{xml_escape(label)}</text>'
        )
    path_svg = ""
    if columns:
        points = path_points(columns, scale_i, scale_j)
        if points:
            coords = " ".join(
                f"{MARGIN_LEFT + PANEL * x / cols:.2f},"
                f"{MARGIN_TOP + PANEL * y / rows:.2f}"
                for x, y in points
            )
            start = points[0]
            sx = MARGIN_LEFT + PANEL * start[0] / cols
            sy = MARGIN_TOP + PANEL * start[1] / rows
            path_svg = (
                f'<polyline points="{coords}" fill="none" stroke="#ffffff" '
                f'stroke-width="3.2" stroke-linejoin="round" stroke-linecap="round"/>'
                f'<polyline points="{coords}" fill="none" stroke="{xml_escape(colour)}" '
                f'stroke-width="1.8" stroke-linejoin="round" stroke-linecap="round" '
                f'marker-end="url(#{xml_escape(marker_id)})"/>'
                f'<circle cx="{sx:.2f}" cy="{sy:.2f}" r="3.2" fill="{xml_escape(colour)}" '
                f'stroke="#ffffff" stroke-width="1"/>'
            )
    return f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}" font-family="ui-sans-serif, system-ui, sans-serif" font-size="11" fill="#212529">
  <defs>
    <marker id="{xml_escape(marker_id)}" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="{xml_escape(colour)}"/>
    </marker>
  </defs>
  <rect width="100%" height="100%" fill="#ffffff"/>
  <text x="{MARGIN_LEFT}" y="16" font-size="13" font-weight="600">{xml_escape(title)}</text>
  <text x="{MARGIN_LEFT}" y="30" fill="#6c757d">{xml_escape(subtitle)}</text>
  <image href="data:image/png;base64,{payload}" x="{MARGIN_LEFT}" y="{MARGIN_TOP}" width="{PANEL}" height="{PANEL}" preserveAspectRatio="none"/>
  <rect x="{MARGIN_LEFT}" y="{MARGIN_TOP}" width="{PANEL}" height="{PANEL}" fill="none" stroke="#dee2e6"/>
  {path_svg}
  {''.join(tick_svg)}
  <text x="{MARGIN_LEFT + PANEL / 2}" y="{height - 8}" text-anchor="middle">target residue</text>
  <text x="14" y="{MARGIN_TOP + PANEL / 2}" text-anchor="middle" transform="rotate(-90 14 {MARGIN_TOP + PANEL / 2})">query residue</text>
  <image href="data:image/png;base64,{cbar}" x="{MARGIN_LEFT + PANEL + 14}" y="{MARGIN_TOP}" width="12" height="{PANEL}" preserveAspectRatio="none"/>
  <rect x="{MARGIN_LEFT + PANEL + 14}" y="{MARGIN_TOP}" width="12" height="{PANEL}" fill="none" stroke="#dee2e6"/>
  <text x="{MARGIN_LEFT + PANEL + 30}" y="{MARGIN_TOP + 8}" font-size="10">{xml_escape(cbar_high)}</text>
  <text x="{MARGIN_LEFT + PANEL + 30}" y="{MARGIN_TOP + PANEL}" font-size="10">{xml_escape(cbar_low)}</text>
</svg>
"""


def cluster_crop(
    row: dict[str, str],
    cluster: dict[str, str] | None,
    query_len: int,
    target_len: int,
    pad: int,
) -> tuple[int, int, int, int]:
    source = cluster or row
    q0, q1 = crop_bounds(int(float(source["query_start"])), int(float(source["query_end"])), query_len, pad)
    t0, t1 = crop_bounds(int(float(source["target_start"])), int(float(source["target_end"])), target_len, pad)
    return q0, q1, t0, t1


def plot_hit(
    row: dict[str, str],
    index: int,
    cluster: dict[str, str] | None,
    query_emb: np.ndarray,
    target_emb: np.ndarray,
    params: ScoringParameters,
    pad: int,
    max_cells: int,
    colour: str,
) -> list[tuple[str, str]]:
    q0, q1, t0, t1 = cluster_crop(row, cluster, query_emb.shape[0], target_emb.shape[0], pad)
    query = query_emb[q0:q1]
    target = target_emb[t0:t1]
    if query.size == 0 or target.size == 0:
        raise ValueError("empty crop")
    cells = query.shape[0] * target.shape[0]
    if cells > max_cells:
        raise ValueError(f"crop {query.shape[0]}x{target.shape[0]} exceeds --max-cells {max_cells}")
    cosine = similarity_matrix(query, target)
    scores = transform_scores(cosine, params)
    result = align(query, target, params=params, max_cells=max_cells)
    prefix = (
        f"{safe_name(str(row.get('cluster_id', index)))}_"
        f"{safe_name(row['query_id'])}__{safe_name(row['target_id'])}"
    )
    crop_note = f"crop q[{q0}:{q1}] × t[{t0}:{t1}]"
    cosine_svg = matrix_svg(
        cosine,
        title=f"cosine  {row['query_id']} → {row['target_id']}",
        subtitle=crop_note,
        query_offset=q0,
        target_offset=t0,
        vmin=float(np.min(cosine)),
        vmax=float(np.max(cosine)),
        stops=VIRIDIS,
        columns=None,
        colour=colour,
        cbar_low=f"{np.min(cosine):.2f}",
        cbar_high=f"{np.max(cosine):.2f}",
        marker_id=f"arrow-{prefix}-cosine",
    )
    score_svg = matrix_svg(
        scores,
        title=f"SW scores  {row['query_id']} → {row['target_id']}",
        subtitle=f"{crop_note}  score={result.score:.2f}",
        query_offset=q0,
        target_offset=t0,
        vmin=float(np.min(scores)),
        vmax=float(np.max(scores)),
        stops=RDBU,
        columns=result.columns,
        colour=colour,
        cbar_low=f"{np.min(scores):.2f}",
        cbar_high=f"{np.max(scores):.2f}",
        marker_id=f"arrow-{prefix}-scores",
    )
    return [
        (f"{prefix}_similarity.svg", cosine_svg),
        (f"{prefix}_scores.svg", score_svg),
    ]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignments", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--query-embeddings", type=Path, nargs="+", required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--pad", type=int, default=32)
    parser.add_argument("--max-cells", type=int, default=16_777_216)
    parser.add_argument("--max-pairs", type=int, default=25)
    parser.add_argument("--highlight-colour", default="#24B064")
    parser.add_argument("--cpus", type=int, default=1)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.pad < 0:
        print("error: --pad must be >= 0", file=sys.stderr)
        return 2
    if args.max_pairs < 1:
        print("error: --max-pairs must be >= 1", file=sys.stderr)
        return 2
    try:
        with args.alignments.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        clusters = load_clusters(args.clusters)
        params = load_parameters(args.parameters)
        query_emb = load_npz_shards(args.query_embeddings)
        target_emb = resolve_targets(args.database)
    except (OSError, ValueError, TypeError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    jobs = []
    for index, row in enumerate(rows[:args.max_pairs]):
        jobs.append((row, index, clusters.get(row.get("cluster_id", ""))))

    args.outdir.mkdir(parents=True, exist_ok=True)
    drawn = 0

    def run(job: tuple) -> list[tuple[str, str]]:
        row, index, cluster = job
        query_id = row["query_id"]
        target_id = row["target_id"]
        try:
            if query_id not in query_emb:
                raise ValueError(f"query {query_id!r} missing from query embeddings")
            if target_id not in target_emb:
                raise ValueError(f"target {target_id!r} missing from database embeddings")
            return plot_hit(
                row,
                index,
                cluster,
                query_emb[query_id],
                target_emb[target_id],
                params,
                args.pad,
                args.max_cells,
                args.highlight_colour,
            )
        except (ValueError, KeyError, TypeError) as exc:
            print(f"warning: skip {query_id} vs {target_id}: {exc}", file=sys.stderr)
            return []

    workers = max(1, min(args.cpus, len(jobs) or 1))
    if jobs:
        align(np.eye(2), np.eye(2), params=params, max_cells=4)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for result in pool.map(run, jobs):
            for name, body in result:
                (args.outdir / name).write_text(body)
                drawn += 1
    print(f"wrote {drawn} SVGs to {args.outdir} ({workers} workers)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
