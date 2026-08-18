#!/usr/bin/env python3
"""Render validated recall@100-versus-QPS benchmark plots and a Markdown table.

Each point is one measured backend/configuration.  Lines are intentionally not
drawn between points: parameter sweeps are observations, not a continuous
performance function.  Recall is the horizontal axis and QPS is the
logarithmic vertical axis because the useful QPS range normally spans orders
of magnitude.
"""
from __future__ import annotations

import argparse
import html
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from results_common import read_records
from validate_results import validate_records


RECALL_THRESHOLD = 0.80
BACKEND_STYLE = {
    "faiss": {"color": "#2F6690", "marker": "circle"},
    "scann": {"color": "#B65D00", "marker": "square"},
    "ngt": {"color": "#68743A", "marker": "diamond"},
    "cuvs": {"color": "#9B3F72", "marker": "triangle"},
}


def _slug(value: Any) -> str:
    result = re.sub(r"[^a-z0-9]+", "-", str(value).lower()).strip("-")
    return result or "unknown"


def _format_number(value: float | int | None, *, precision: int = 2) -> str:
    if value is None:
        return "—"
    value = float(value)
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.{precision}f}M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:.{precision}f}k"
    if abs(value) >= 100:
        return f"{value:.0f}"
    if abs(value) >= 10:
        return f"{value:.1f}"
    return f"{value:.{precision}f}"


def _format_bytes(value: float | int | None) -> str:
    if value is None:
        return "—"
    value = float(value)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    for unit in units:
        if abs(value) < 1024 or unit == units[-1]:
            return f"{value:.2f} {unit}" if unit != "B" else f"{value:.0f} B"
        value /= 1024
    raise AssertionError("unreachable")


def _raw_database_bytes(row: dict[str, Any]) -> int:
    """Return the normalized float32 database payload represented by a row."""

    return int(row["dataset_window_count"]) * int(row["dimension"]) * 4


def _escape(value: Any) -> str:
    return html.escape(str(value), quote=True)


def _markdown_text(value: Any, *, line_break: str = "<br>") -> str:
    """Escape untrusted result metadata for a Markdown table or heading."""

    text = str(value).replace("\r\n", "\n").replace("\r", "\n")
    text = html.escape(text, quote=False).replace("\\", "\\\\")
    for character in ("`", "*", "_", "[", "]", "|"):
        text = text.replace(character, f"\\{character}")
    return text.replace("\n", line_break)


def _recall_ticks(threshold: float) -> tuple[list[float], int]:
    """Return readable tick values and an exact-enough display precision."""

    span = 1.0 - threshold
    target_step = span / 4.0
    exponent = math.floor(math.log10(target_step))
    scale = 10.0**exponent
    step = next(multiplier * scale for multiplier in (1.0, 2.0, 2.5, 5.0, 10.0) if multiplier * scale >= target_step)

    def decimal_places(value: float) -> int:
        for places in range(11):
            scaled = value * (10**places)
            if math.isclose(scaled, round(scaled), rel_tol=0.0, abs_tol=1e-9):
                return places
        return 10

    precision = max(2, decimal_places(step), decimal_places(threshold))
    epsilon = step * 1e-8
    ticks = [threshold]
    multiple = math.floor((threshold + epsilon) / step) + 1
    while (value := multiple * step) < 1.0 - epsilon:
        ticks.append(round(value, 12))
        multiple += 1
    if not math.isclose(ticks[-1], 1.0, rel_tol=0.0, abs_tol=epsilon):
        ticks.append(1.0)
    return ticks, precision


def _duplicate_parameter_labels(rows: Iterable[dict[str, Any]]) -> set[tuple[str, str]]:
    """Identify labels that need a run ID to stay distinguishable in an SVG."""

    counts: dict[tuple[str, str], int] = defaultdict(int)
    for row in rows:
        counts[(str(row["backend"]), str(row["parameter_label"]))] += 1
    return {key for key, count in counts.items() if count > 1}


def _svg_marker(kind: str, x: float, y: float, color: str, size: float = 6.0) -> str:
    if kind == "square":
        return f'<rect x="{x - size:.2f}" y="{y - size:.2f}" width="{2 * size:.2f}" height="{2 * size:.2f}" fill="{color}" stroke="#20252B" stroke-width="1" />'
    if kind == "diamond":
        points = f"{x:.2f},{y - size:.2f} {x + size:.2f},{y:.2f} {x:.2f},{y + size:.2f} {x - size:.2f},{y:.2f}"
        return f'<polygon points="{points}" fill="{color}" stroke="#20252B" stroke-width="1" />'
    if kind == "triangle":
        points = f"{x:.2f},{y - size:.2f} {x + size:.2f},{y + size:.2f} {x - size:.2f},{y + size:.2f}"
        return f'<polygon points="{points}" fill="{color}" stroke="#20252B" stroke-width="1" />'
    return f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{size:.2f}" fill="{color}" stroke="#20252B" stroke-width="1" />'


def _qps_ticks(minimum: float, maximum: float) -> tuple[float, float, list[float]]:
    """Choose a small set of readable log-scale ticks covering a QPS range."""

    if minimum <= 0 or maximum <= 0:
        raise ValueError("QPS must be positive for a logarithmic plot")
    if math.isclose(minimum, maximum):
        minimum /= 1.8
        maximum *= 1.8
    lower_log = math.log10(minimum)
    upper_log = math.log10(maximum)
    padding = max(0.08, (upper_log - lower_log) * 0.08)
    lower_log -= padding
    upper_log += padding
    ticks: list[float] = []
    for exponent in range(math.floor(lower_log) - 1, math.ceil(upper_log) + 2):
        for multiplier in (1, 2, 5):
            value = multiplier * (10**exponent)
            if 10**lower_log <= value <= 10**upper_log:
                ticks.append(float(value))
    return lower_log, upper_log, sorted(set(ticks))


def _label_positions(points: list[tuple[float, float, str]]) -> dict[str, float]:
    """Prevent direct labels from overlapping without changing the data position."""

    if not points:
        return {}
    ordered = sorted(points, key=lambda item: (item[1], item[2]))
    positions: dict[str, float] = {}
    previous = -math.inf
    for _, y, key in ordered:
        positions[key] = max(y, previous + 14)
        previous = positions[key]
    ceiling = max(y for _, y, _ in points) + 46
    overflow = max(0.0, max(positions.values()) - ceiling)
    if overflow:
        for key in positions:
            positions[key] -= overflow
    return positions


def render_scope_svg(
    scope: list[dict[str, Any]],
    output: Path,
    source_note: str,
    *,
    recall_threshold: float = RECALL_THRESHOLD,
    show_all: bool = False,
) -> None:
    """Write one publication-ready SVG for a comparable dataset/hardware scope.

    The default is the requested high-recall view. ``show_all`` writes a
    second full-range view so low-recall measurements remain visible without
    changing the strict thresholded plot.
    """

    if not scope:
        raise ValueError("cannot render an empty scope")
    if not 0.0 <= recall_threshold < 1.0:
        raise ValueError("recall_threshold must be in [0, 1)")
    if not show_all and any(float(row["recall_median"]) <= recall_threshold for row in scope):
        raise ValueError("scope contains a point at or below the strict recall threshold")
    ordered_scope = sorted(scope, key=lambda item: (item["backend"], item["parameter_label"], item["run_id"]))
    first = ordered_scope[0]
    width, height = 1360, 820
    left, top, plot_width, plot_height = 120, 154, 850, 520
    right = left + plot_width
    bottom = top + plot_height
    qps_min = min(float(row["qps_min"]) for row in ordered_scope)
    qps_max = max(float(row["qps_max"]) for row in ordered_scope)
    log_min, log_max, y_ticks = _qps_ticks(qps_min, qps_max)
    x_floor = 0.0 if show_all else recall_threshold
    x_ticks, recall_precision = _recall_ticks(x_floor)

    def x_position(value: float) -> float:
        return left + ((value - x_floor) / (1.0 - x_floor)) * plot_width

    def y_position(value: float) -> float:
        return bottom - ((math.log10(value) - log_min) / (log_max - log_min)) * plot_height

    title = f"Recall@100 vs. queries/s — {first['dataset_id']}"
    subtitle = (
        f"{first['dataset_window_count']:,} indexed windows · {first['dimension']:,} dimensions · "
        f"hardware: {first['hardware_id']} · metric: cosine · query batch: {first.get('query_batch_size', 'unknown')}"
    )
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-labelledby="title description">',
        f"<title id=\"title\">{_escape(title)}</title>",
        '<desc id="description">Measured recall@100 versus queries per second. Points are median repeats and vertical whiskers show the measured QPS range. No interpolation is drawn.</desc>',
        '<rect width="100%" height="100%" fill="#FFFFFF" />',
        '<style>text { font-family: Inter, Arial, sans-serif; fill: #20252B; } .mono { font-family: "SFMono-Regular", Consolas, monospace; } .grid { stroke: #D9DEE5; stroke-width: 1; } .axis { stroke: #2D333B; stroke-width: 1.4; } .note { fill: #4E5965; font-size: 14px; }</style>',
        f'<text x="{left}" y="50" font-size="28" font-weight="700">{_escape(title)}</text>',
        f'<text x="{left}" y="78" font-size="16" fill="#4E5965">{_escape(subtitle)}</text>',
        (
            f'<text x="{left}" y="106" font-size="14" fill="#4E5965">All validated points are shown. Point = median of timed repeats; vertical whisker = measured min–max QPS.</text>'
            if show_all
            else f'<text x="{left}" y="106" font-size="14" fill="#4E5965">Only recall@100 &gt; {recall_threshold:.{recall_precision}f} is shown. Point = median of timed repeats; vertical whisker = measured min–max QPS.</text>'
        ),
    ]
    for tick in x_ticks:
        x = x_position(tick)
        lines.extend(
            [
                f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{bottom}" />',
                f'<text class="mono" x="{x:.2f}" y="{bottom + 26}" font-size="13" text-anchor="middle">{tick:.{recall_precision}f}</text>',
            ]
        )
    for tick in y_ticks:
        y = y_position(tick)
        lines.extend(
            [
                f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{right}" y2="{y:.2f}" />',
                f'<text class="mono" x="{left - 14}" y="{y + 5:.2f}" font-size="13" text-anchor="end">{_escape(_format_number(tick))}</text>',
            ]
        )
    lines.extend(
        [
            f'<line class="axis" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" />',
            f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" />',
            f'<text x="{left + plot_width / 2:.2f}" y="{bottom + 62}" font-size="16" font-weight="600" text-anchor="middle">Recall@100</text>',
            f'<text x="32" y="{top + plot_height / 2:.2f}" font-size="16" font-weight="600" text-anchor="middle" transform="rotate(-90 32 {top + plot_height / 2:.2f})">Queries per second (log scale)</text>',
            f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#B7C0CA" stroke-width="1" />',
        ]
    )
    legend_backends = [backend for backend in BACKEND_STYLE if any(str(row["backend"]) == backend for row in ordered_scope)]
    legend_x = left
    legend_y = 132
    lines.append(f'<text x="{legend_x}" y="{legend_y}" font-size="14" font-weight="700">Backend:</text>')
    cursor = legend_x + 76
    for backend in legend_backends:
        style = BACKEND_STYLE[backend]
        lines.append(_svg_marker(style["marker"], cursor, legend_y - 5, style["color"], size=5.0))
        lines.append(f'<text x="{cursor + 11}" y="{legend_y}" font-size="13" fill="{style["color"]}">{_escape(backend.upper())}</text>')
        cursor += 112
    label_points: list[tuple[float, float, str]] = []
    plotted: list[tuple[dict[str, Any], float, float, str]] = []
    duplicate_labels = _duplicate_parameter_labels(ordered_scope)
    for index, row in enumerate(ordered_scope):
        x = x_position(float(row["recall_median"]))
        y = y_position(float(row["qps_median"]))
        key = str(index)
        label_points.append((x, y, key))
        plotted.append((row, x, y, key))
    labels = _label_positions(label_points)
    for row, x, y, key in plotted:
        style = BACKEND_STYLE[row["backend"]]
        low_y = y_position(float(row["qps_min"]))
        high_y = y_position(float(row["qps_max"]))
        label_y = labels[key]
        label_x = right + 27
        parameter_text = str(row["parameter_label"])
        if (str(row["backend"]), parameter_text) in duplicate_labels:
            parameter_text = f"{parameter_text} ({row['run_id']})"
        lines.extend(
            [
                f'<line x1="{x:.2f}" y1="{low_y:.2f}" x2="{x:.2f}" y2="{high_y:.2f}" stroke="{style["color"]}" stroke-width="2" />',
                f'<line x1="{x - 4:.2f}" y1="{low_y:.2f}" x2="{x + 4:.2f}" y2="{low_y:.2f}" stroke="{style["color"]}" stroke-width="2" />',
                f'<line x1="{x - 4:.2f}" y1="{high_y:.2f}" x2="{x + 4:.2f}" y2="{high_y:.2f}" stroke="{style["color"]}" stroke-width="2" />',
                _svg_marker(style["marker"], x, y, style["color"]),
                f'<line x1="{x + 8:.2f}" y1="{y:.2f}" x2="{label_x - 6:.2f}" y2="{label_y - 4:.2f}" stroke="#AAB3BD" stroke-width="1" />',
                f'<text x="{label_x}" y="{label_y:.2f}" font-size="13"><tspan font-weight="700" fill="{style["color"]}">{_escape(row["backend"])}</tspan><tspan> / {_escape(parameter_text)}</tspan></text>',
            ]
        )
    lines.append(f'<text x="{left}" y="{height - 72}" class="note">Source: {_escape(source_note)}</text>')
    lines.append(f'<text x="{left}" y="{height - 48}" class="note">Comparison identity: query set {_escape(str(first["query_ids_sha256"])[:12])}… · exact ground truth {_escape(str(first["ground_truth_ids_sha256"])[:12])}…</text>')
    lines.append('</svg>')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _markdown_table(rows: Iterable[dict[str, Any]], scope_title: str) -> list[str]:
    lines = [
        f"## {scope_title}",
        "",
        "| Backend | Configuration | Run ID | Raw database | Recall@100, median (min–max) | QPS, median (min–max) | Timed repeats | Index size | Build time | Build peak RSS | Build peak VRAM | Search peak RSS | Search peak VRAM |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in sorted(rows, key=lambda item: (item["backend"], item["parameter_label"], item["run_id"])):
        measurement = row.get("measurement") or {}
        recall = f"{row['recall_median']:.4f} ({row['recall_min']:.4f}–{row['recall_max']:.4f})"
        qps = f"{_format_number(row['qps_median'])} ({_format_number(row['qps_min'])}–{_format_number(row['qps_max'])})"
        lines.append(
            "| "
            + " | ".join(
                (
                    _markdown_text(row["backend"]),
                    _markdown_text(row["parameter_label"]),
                    _markdown_text(row["run_id"]),
                    _format_bytes(_raw_database_bytes(row)),
                    recall,
                    qps,
                    str(row["repeat_count"]),
                    _format_bytes(row.get("index_bytes")),
                    f"{_format_number(row.get('build_seconds'))} s",
                    _format_bytes(measurement.get("build_peak_rss_bytes")),
                    _format_bytes(measurement.get("build_peak_vram_bytes")),
                    _format_bytes(row.get("peak_rss_bytes")),
                    _format_bytes(row.get("peak_vram_bytes")),
                )
            )
            + " |"
        )
    return lines


def render_markdown_summary(
    report: dict[str, Any],
    output: Path,
    plot_names: dict[str, dict[str, str]],
    threshold: float,
) -> None:
    """Write a compact results table without claiming a winner from synthetic data."""

    aggregates = report["aggregates"]
    by_scope: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in aggregates:
        by_scope[row["scope_id"]].append(row)
    lines = [
        "# GINflow vector-index benchmark results",
        "",
        f"Validation generated: {_markdown_text(report['generated_at'], line_break=' ')}",
        "",
        "Each row aggregates timed repeats for one backend/configuration. QPS is the median; the parenthesized range is the measured min–max. Curves are not interpolated.",
        "",
        "All successful configurations appear in the tables and full-range SVG plots. "
        f"A second focused SVG per scope contains only configurations with median recall@100 > {threshold:.2f}.",
        "",
        f"Validation: {report['valid_successful_records']} successful repeat records across {report['successful_configurations']} configurations; {len(report['unavailable'])} skipped/error rows retained below.",
        "",
    ]
    for scope_id in sorted(by_scope):
        rows = by_scope[scope_id]
        first = rows[0]
        title = f"{_markdown_text(first['dataset_id'], line_break=' ')} — {_markdown_text(first['hardware_id'], line_break=' ')}"
        lines.extend(_markdown_table(rows, title))
        if scope_id in plot_names:
            names = plot_names[scope_id]
            lines.extend(["", f"All measured points: [{names['all']}]({names['all']})"])
            if "high_recall" in names:
                lines.append(f"Focused recall@100 > {threshold:.2f}: [{names['high_recall']}]({names['high_recall']})")
        excluded = [row for row in rows if row["recall_median"] <= threshold]
        if excluded:
            lines.extend(
                [
                    "",
                    "Below or equal to the focused recall threshold (retained in the table and full-range plot): "
                    + ", ".join(
                        f"{_markdown_text(row['backend'])}/{_markdown_text(row['parameter_label'])} ({row['recall_median']:.4f})"
                        for row in excluded
                    )
                    + ".",
                ]
            )
        lines.append("")
    if report["unavailable"]:
        lines.extend(
            [
                "## Skipped or errored configurations",
                "",
                "| Status | Backend | Dataset | Configuration | Reason |",
                "|---|---|---|---|---|",
            ]
        )
        for row in report["unavailable"]:
            lines.append(
                "| " + " | ".join(
                    _markdown_text(row.get(key) or "—")
                    for key in ("status", "backend", "dataset_id", "parameter_label", "error")
                ) + " |"
            )
        lines.append("")
    lines.extend(
        [
            "## Reproducibility guardrails",
            "",
            "- Rows in one plot share dataset window count, embedding dimension, cosine metric, k=100, warm-up/timed query counts, fixed query batch size, query-set SHA-256, exact-ground-truth SHA-256, cache provenance, runner revision, and hardware identifier.",
            "- The table retains skipped/error runs so unavailable backends are not silently treated as zero throughput or omitted from coverage review.",
            "- See `benchmarks/results-template.md` for the durable report table shape when no real records have been collected yet.",
            "",
        ]
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="result files or directories to scan recursively")
    parser.add_argument("--output-dir", type=Path, required=True, help="directory for SVGs, validation.json, and benchmark_results.md")
    parser.add_argument("--min-repeats", type=int, default=2, help="minimum valid timed repeats per configuration (default: 2)")
    parser.add_argument("--recall-threshold", type=float, default=RECALL_THRESHOLD, help="strict plotting threshold (default: 0.80)")
    parser.add_argument("--source-note", default="validated GINflow benchmark result JSON", help="source note written in SVG footers")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.min_repeats < 1:
        raise SystemExit("--min-repeats must be at least 1")
    if not 0 <= args.recall_threshold < 1:
        raise SystemExit("--recall-threshold must be in [0, 1)")
    records, read_issues, ignored = read_records(args.inputs)
    report = validate_records(records, read_issues=read_issues, min_repeats=args.min_repeats)
    report["inputs"] = [str(Path(path)) for path in args.inputs]
    report["ignored_files"] = ignored
    args.output_dir.mkdir(parents=True, exist_ok=True)
    validation_path = args.output_dir / "validation.json"
    validation_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not report["valid"]:
        print(f"validation failed; inspect {validation_path}", file=sys.stderr)
        for issue in report["errors"]:
            print(f"ERROR {issue['code']}: {issue['message']}", file=sys.stderr)
        return 2

    by_scope: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in report["aggregates"]:
        by_scope[row["scope_id"]].append(row)
    plot_names: dict[str, dict[str, str]] = {}
    for scope_id in sorted(by_scope):
        scope = by_scope[scope_id]
        first = scope[0]
        stem = f"{_slug(first['dataset_id'])}_{_slug(first['hardware_id'])}_{scope_id}"
        all_filename = f"recall_qps_all_{stem}.svg"
        render_scope_svg(
            scope,
            args.output_dir / all_filename,
            args.source_note,
            recall_threshold=args.recall_threshold,
            show_all=True,
        )
        names = {"all": all_filename}
        high_recall_scope = [row for row in scope if row["recall_median"] > args.recall_threshold]
        if high_recall_scope:
            high_recall_filename = f"recall_qps_gt_{args.recall_threshold:.2f}_{stem}.svg"
            render_scope_svg(
                high_recall_scope,
                args.output_dir / high_recall_filename,
                args.source_note,
                recall_threshold=args.recall_threshold,
            )
            names["high_recall"] = high_recall_filename
        plot_names[scope_id] = names
    summary_path = args.output_dir / "benchmark_results.md"
    render_markdown_summary(report, summary_path, plot_names, args.recall_threshold)
    plot_count = sum(len(names) for names in plot_names.values())
    print(
        f"wrote {plot_count} SVG plot(s), {summary_path}, and {validation_path}; "
        f"{len(report['aggregates'])} configuration aggregate(s) validated"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
