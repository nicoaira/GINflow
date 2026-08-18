#!/usr/bin/env python3
"""Render a self-contained live dashboard from the external benchmark cache.

The dashboard deliberately reads partial result JSON files without requiring a
configuration to have all repeats.  Run it with ``--watch`` while the backend
containers are active; the output is replaced atomically so a browser can
refresh the same HTML file safely.
"""
from __future__ import annotations

import argparse
import html
import json
import math
import os
import re
import statistics
import subprocess
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


BACKENDS = ("faiss", "scann", "ngt", "cuvs")
COLORS = {"faiss": "#4f8db3", "scann": "#d5843b", "ngt": "#87954e", "cuvs": "#bd5c91"}
EXPECTED_REPEATS = 3
# Default ``--plan frontier`` counts used for progress estimates.  The runner
# writes terminal JSON only after a configuration finishes, so progress is
# deliberately labelled as an estimate and is never treated as a measurement.
FRONTIER_CONFIGURATION_COUNTS = {"faiss": 23, "scann": 18, "ngt": 3, "cuvs": 40}

INTERACTION_SCRIPT = r"""
<script>
(() => {
  const tooltip = document.getElementById('chart-tooltip');
  const checkbox = document.getElementById('gpu-highlight');
  const paretoCheckbox = document.getElementById('pareto-toggle');
  if (!tooltip || !checkbox || !paretoCheckbox) return;

  const escapeHtml = (value) => String(value).replace(/[&<>"']/g, (character) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
  })[character]);
  const formatValue = (value) => {
    if (value === null || value === undefined || value === '') return '—';
    return typeof value === 'object' ? JSON.stringify(value, null, 2) : String(value);
  };
  const moveTooltip = (event, point) => {
    const bounds = point.getBoundingClientRect();
    const clientX = event.clientX || (bounds.left + bounds.width / 2);
    const clientY = event.clientY || (bounds.top + bounds.height / 2);
    const gap = 14;
    const maxLeft = Math.max(8, window.innerWidth - tooltip.offsetWidth - 8);
    const maxTop = Math.max(8, window.innerHeight - tooltip.offsetHeight - 8);
    tooltip.style.left = `${Math.min(clientX + gap, maxLeft)}px`;
    tooltip.style.top = `${Math.min(clientY + gap, maxTop)}px`;
  };
  const showTooltip = (event) => {
    const point = event.currentTarget;
    let info;
    try { info = JSON.parse(point.dataset.info); }
    catch (_) { info = { error: 'Unable to decode point metadata' }; }
    tooltip.innerHTML = Object.entries(info).map(([key, value]) =>
      `<div class="tooltip-row"><b>${escapeHtml(key)}</b><pre>${escapeHtml(formatValue(value))}</pre></div>`
    ).join('');
    tooltip.hidden = false;
    moveTooltip(event, point);
  };
  const hideTooltip = () => { tooltip.hidden = true; };
  document.querySelectorAll('.chart-point').forEach((point) => {
    point.addEventListener('pointerenter', showTooltip);
    point.addEventListener('pointermove', (event) => moveTooltip(event, point));
    point.addEventListener('pointerleave', hideTooltip);
    point.addEventListener('focus', showTooltip);
    point.addEventListener('blur', hideTooltip);
  });
  const applyGpuHighlight = () => document.body.classList.toggle('gpu-highlight', checkbox.checked);
  const svgNamespace = 'http://www.w3.org/2000/svg';
  const paretoFrontier = (points) => points.filter((point, index) => !points.some((other, otherIndex) => {
    if (index === otherIndex) return false;
    const atLeastAsGood = other.recall >= point.recall && other.qps >= point.qps;
    const strictlyBetter = other.recall > point.recall || other.qps > point.qps;
    return atLeastAsGood && strictlyBetter;
  })).sort((left, right) => left.recall - right.recall);
  const drawParetoCurves = () => {
    document.querySelectorAll('svg.chart').forEach((svg) => {
      svg.querySelectorAll('.pareto-overlay').forEach((overlay) => overlay.remove());
      const byBackend = new Map();
      svg.querySelectorAll('.chart-point').forEach((point) => {
        try {
          const info = JSON.parse(point.dataset.info);
          const recall = Number(info.recall_at_100_median);
          const qps = Number(info.queries_per_second_median);
          if (!info.backend || !Number.isFinite(recall) || !Number.isFinite(qps)) return;
          const circle = point.querySelector('circle');
          if (!circle) return;
          if (!byBackend.has(info.backend)) byBackend.set(info.backend, []);
          byBackend.get(info.backend).push({
            recall, qps,
            x: circle.getAttribute('cx'),
            y: circle.getAttribute('cy'),
            color: circle.getAttribute('fill') || '#ffffff'
          });
        } catch (_) { /* Ignore a malformed point and keep the other curves. */ }
      });
      const overlay = document.createElementNS(svgNamespace, 'g');
      overlay.setAttribute('class', 'pareto-overlay');
      overlay.setAttribute('aria-label', 'Pareto frontiers by backend');
      byBackend.forEach((backendPoints, backend) => {
        const frontier = paretoFrontier(backendPoints);
        if (frontier.length < 2) return;
        const line = document.createElementNS(svgNamespace, 'polyline');
        line.setAttribute('class', `pareto-line pareto-${backend}`);
        line.setAttribute('points', frontier.map((point) => `${point.x},${point.y}`).join(' '));
        line.setAttribute('fill', 'none');
        line.setAttribute('stroke', frontier[0].color);
        line.setAttribute('stroke-width', '3');
        line.setAttribute('stroke-dasharray', '8 5');
        line.setAttribute('stroke-linejoin', 'round');
        line.setAttribute('stroke-linecap', 'round');
        line.setAttribute('opacity', '0.95');
        line.setAttribute('pointer-events', 'none');
        line.setAttribute('data-backend', backend);
        overlay.appendChild(line);
      });
      overlay.style.display = paretoCheckbox.checked ? '' : 'none';
      svg.appendChild(overlay);
    });
  };
  paretoCheckbox.addEventListener('change', drawParetoCurves);
  checkbox.addEventListener('change', applyGpuHighlight);
  applyGpuHighlight();
  drawParetoCurves();
})();
</script>
"""


def _escape(value: Any) -> str:
    return html.escape(str(value), quote=True)


def _bytes(value: Any) -> int | None:
    try:
        if value is None:
            return None
        return max(0, int(value))
    except (TypeError, ValueError):
        return None


def _number(value: Any) -> float | None:
    try:
        if value is None:
            return None
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _fmt_bytes(value: Any) -> str:
    amount = _bytes(value)
    if amount is None:
        return "—"
    return f"{amount / 1024**3:.2f} GiB"


def _fmt_number(value: Any, digits: int = 2) -> str:
    amount = _number(value)
    if amount is None:
        return "—"
    if abs(amount) >= 1000:
        return f"{amount:,.0f}"
    return f"{amount:.{digits}f}"


def _fmt_duration(seconds: Any) -> str:
    amount = _number(seconds)
    if amount is None:
        return "—"
    amount = max(0, int(amount))
    hours, remainder = divmod(amount, 3600)
    minutes, seconds_value = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes:02d}m {seconds_value:02d}s"
    return f"{minutes}m {seconds_value:02d}s"


def _median(values: Iterable[Any]) -> float | None:
    numbers = [number for value in values if (number := _number(value)) is not None]
    return float(statistics.median(numbers)) if numbers else None


def _read_object(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def _cache_state(dataset_dir: Path) -> dict[str, Any]:
    request = _read_object(dataset_dir / "window-cache-request.json") or {}
    flatten = _read_object(dataset_dir / "flat" / "flatten-manifest.json") or {}
    queries = _read_object(dataset_dir / "queries" / "query-selection.json") or {}
    truth = _read_object(dataset_dir / "ground-truth" / "ground-truth.json") or {}
    vectors = flatten.get("vectors") if isinstance(flatten.get("vectors"), dict) else {}
    return {
        "request_state": request.get("state", "not-started"),
        "request_id": request.get("request_id"),
        "window_count": (flatten.get("windows") or {}).get("count") if isinstance(flatten.get("windows"), dict) else None,
        "record_count": (flatten.get("records") or {}).get("count") if isinstance(flatten.get("records"), dict) else None,
        "dimension": (vectors.get("shape") or [None, None])[1] if isinstance(vectors, dict) else None,
        "query_count": queries.get("query_count"),
        "query_selection_id": queries.get("query_selection_id"),
        "ground_truth_cache_id": truth.get("ground_truth_cache_id"),
        "ground_truth_engine": ((truth.get("engine") or {}).get("used") if isinstance(truth.get("engine"), dict) else None),
        "ground_truth_k": truth.get("k"),
    }


def _is_gpu_result(payload: dict[str, Any], parameters: dict[str, Any], provenance: dict[str, Any]) -> bool:
    """Classify whether a result used a GPU, using backend and recorded metadata."""

    backend = str(payload.get("backend", "")).lower()
    if backend == "cuvs":
        return True
    device = str(parameters.get("device", "")).lower()
    if device in {"gpu", "cuda"} or parameters.get("use_gpu") is True:
        return True
    if provenance.get("gpu_device") is not None or provenance.get("gpu") is True:
        return True
    return False


def _result_row(payload: dict[str, Any], source: Path) -> dict[str, Any] | None:
    if payload.get("status") not in {"ok", "skipped", "error"}:
        return None
    measurement = payload.get("measurement") if isinstance(payload.get("measurement"), dict) else {}
    capacity = measurement.get("capacity_preflight") if isinstance(measurement.get("capacity_preflight"), dict) else None
    parameters = payload.get("parameters") if isinstance(payload.get("parameters"), dict) else {}
    provenance = payload.get("provenance") if isinstance(payload.get("provenance"), dict) else {}
    return {
        "status": payload.get("status"),
        "backend": payload.get("backend", "unknown"),
        "parameter_label": payload.get("parameter_label", source.stem),
        "run_id": payload.get("run_id", "unknown"),
        "repeat": payload.get("repeat"),
        "dataset_window_count": payload.get("dataset_window_count"),
        "dimension": payload.get("dimension"),
        "qps": payload.get("qps"),
        "recall": payload.get("recall_at_100"),
        "index_bytes": payload.get("index_bytes"),
        "build_seconds": payload.get("build_seconds"),
        "build_rss": measurement.get("build_peak_rss_bytes"),
        "build_vram": measurement.get("build_peak_vram_bytes"),
        "search_rss": payload.get("peak_rss_bytes"),
        "search_vram": payload.get("peak_vram_bytes"),
        "capacity": capacity,
        "parameters": parameters,
        "provenance": provenance,
        "measurement": measurement,
        "metric": payload.get("metric"),
        "k": payload.get("k"),
        "warmup_queries": payload.get("warmup_queries"),
        "timed_queries": payload.get("timed_queries"),
        "latency": payload.get("latency_ms"),
        "is_gpu": _is_gpu_result(payload, parameters, provenance),
        "error": payload.get("error"),
        "source": str(source),
    }


def _command_arg(command: str, name: str, default: str | None = None) -> str | None:
    match = re.search(rf"(?:^|\s){re.escape(name)}(?:\s+|=)([^\s]+)", command)
    return match.group(1).strip("'\"") if match else default


def _parse_elapsed(value: str) -> int | None:
    value = value.strip()
    try:
        days = 0
        if "-" in value:
            day_text, value = value.split("-", 1)
            days = int(day_text)
        fields = [int(part) for part in value.split(":")]
        if len(fields) == 3:
            hours, minutes, seconds = fields
        elif len(fields) == 2:
            hours, minutes, seconds = 0, *fields
        else:
            return None
        return days * 86400 + hours * 3600 + minutes * 60 + seconds
    except (TypeError, ValueError):
        return None


def _parse_started_at(value: str | None) -> float | None:
    if not value:
        return None
    try:
        return datetime.fromisoformat(value.strip().replace("Z", "+00:00")).timestamp()
    except ValueError:
        return None


def _expected_configuration_count(backend: str, command: str) -> int | None:
    if _command_arg(command, "--plan", "frontier") != "frontier":
        return None
    only = _command_arg(command, "--only")
    if not only:
        return FRONTIER_CONFIGURATION_COUNTS.get(backend)
    requested = {value.strip().lower() for value in only.split(",") if value.strip()}
    if backend == "ngt":
        return len(requested.intersection({"ngt", "qg", "qbg"}))
    if backend == "faiss":
        counts = {"flatip": 1, "ivfflat": 8, "ivfpq": 10, "hnsw": 4}
        return sum(counts.get(value, 0) for value in requested)
    if backend == "cuvs":
        counts = {"cagra": 5, "ivf": 10, "ivf-pq": 25}
        return sum(counts.get(value, 0) for value in requested)
    return None


def _progress_label_allowed(backend: str, label: str, command: str | None) -> bool:
    """Keep progress scoped to the settings selected by the active runner."""

    if not command:
        return True
    only = _command_arg(command, "--only")
    if not only:
        return True
    requested = {value.strip().lower() for value in only.split(",") if value.strip()}
    normalized = label.strip().lower()
    if backend == "ngt":
        return any(normalized == f"pipeline-{value}" for value in requested)
    if backend == "faiss":
        return any(re.search(rf"(?:^|-){re.escape(value)}(?:-|$)", normalized) for value in requested)
    if backend == "cuvs":
        aliases = {"ivfpq": "ivf-pq"}
        return any(normalized.startswith(aliases.get(value, value)) for value in requested)
    return True


def _result_progress(
    dataset_dir: Path,
    backend: str,
    repeats: int,
    expected_configurations: int | None,
    command: str | None = None,
    result_mtime_after: float | None = None,
) -> dict[str, Any]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    record_count = 0
    result_dir = dataset_dir / "results" / backend
    for path in sorted(result_dir.glob("*.json")) if result_dir.is_dir() else ():
        if result_mtime_after is not None:
            try:
                if path.stat().st_mtime < result_mtime_after:
                    continue
            except OSError:
                continue
        row = _result_row(_read_object(path) or {}, path)
        if row is None:
            continue
        if not _progress_label_allowed(backend, str(row["parameter_label"]), command):
            continue
        record_count += 1
        groups[(str(row["parameter_label"]), str(row["run_id"]))].append(row)
    terminal_configurations = 0
    for rows in groups.values():
        ok_count = sum(row["status"] == "ok" for row in rows)
        if ok_count >= repeats or any(row["status"] in {"error", "skipped"} for row in rows):
            terminal_configurations += 1
    progress = None
    if expected_configurations and expected_configurations > 0:
        progress = min(1.0, terminal_configurations / expected_configurations)
    return {
        "terminal_configurations": terminal_configurations,
        "expected_configurations": expected_configurations,
        "record_count": record_count,
        "progress": progress,
    }


def _docker_command(args: list[str], timeout: float = 3.0) -> str | None:
    try:
        completed = subprocess.run(args, capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode:
        return None
    return completed.stdout.strip()


def _running_jobs(cache_dir: Path) -> list[dict[str, Any]]:
    """Discover active backend runner containers and best-effort resource progress."""

    names_output = _docker_command(["docker", "ps", "--format", "{{.Names}}\t{{.Status}}"])
    if not names_output:
        return []
    jobs: list[dict[str, Any]] = []
    runner_pattern = re.compile(r"benchmarks/run_(faiss|scann|ngt|cuvs)\.py")
    for line in names_output.splitlines():
        name, _, status_text = line.partition("\t")
        if not name:
            continue
        command_json = _docker_command(["docker", "inspect", "-f", "{{json .Config.Cmd}}", name])
        try:
            command_parts = json.loads(command_json or "null")
        except json.JSONDecodeError:
            command_parts = None
        if not isinstance(command_parts, list):
            continue
        command = " ".join(str(part) for part in command_parts)
        match = runner_pattern.search(command)
        if not match:
            continue
        backend = match.group(1)
        dataset_id = _command_arg(command, "--dataset", "unknown") or "unknown"
        dataset_dir = cache_dir / dataset_id
        if not dataset_dir.is_dir():
            continue
        repeats = int(_command_arg(command, "--repeats", str(EXPECTED_REPEATS)) or EXPECTED_REPEATS)
        expected = _expected_configuration_count(backend, command)
        started_at = _parse_started_at(
            _docker_command(["docker", "inspect", "-f", "{{.State.StartedAt}}", name])
        )
        result_mtime_after = started_at if "--no-resume" in command.split() else None
        progress = _result_progress(dataset_dir, backend, repeats, expected, command, result_mtime_after)
        top_output = _docker_command(["docker", "top", name, "-eo", "pid,etime,pcpu,pmem,rss,stat,comm,args"])
        processes: list[dict[str, Any]] = []
        if top_output:
            for process_line in top_output.splitlines()[1:]:
                fields = process_line.split(None, 7)
                if len(fields) < 8:
                    continue
                try:
                    processes.append(
                        {
                            "pid": fields[0],
                            "elapsed": fields[1],
                            "elapsed_seconds": _parse_elapsed(fields[1]),
                            "cpu": float(fields[2]),
                            "memory_percent": float(fields[3]),
                            "rss_bytes": int(fields[4]) * 1024,
                            "stat": fields[5],
                            "comm": fields[6],
                            "args": fields[7],
                        }
                    )
                except (TypeError, ValueError):
                    continue
        stage_process = next(
            (process for process in processes if process["comm"] not in {"python", "python3", "bash", "sh"}),
            processes[0] if processes else None,
        )
        image = _docker_command(["docker", "inspect", "-f", "{{.Config.Image}}", name]) or "unknown image"
        jobs.append(
            {
                "container": name,
                "status": status_text or "running",
                "backend": backend,
                "dataset": dataset_id,
                "command": command,
                "image": image,
                "repeats": repeats,
                "progress": progress,
                "processes": processes,
                "stage": stage_process["args"] if stage_process else command,
                "cpu": sum(process["cpu"] for process in processes) if processes else None,
                "rss_bytes": max((process["rss_bytes"] for process in processes), default=None),
                "elapsed_seconds": max((process["elapsed_seconds"] or 0 for process in processes), default=None),
            }
        )
    return jobs


def _collect_dataset(dataset_dir: Path) -> dict[str, Any]:
    grouped: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    raw_counts: Counter[str] = Counter()
    for backend in BACKENDS:
        result_dir = dataset_dir / "results" / backend
        for path in sorted(result_dir.glob("*.json")) if result_dir.is_dir() else ():
            row = _result_row(_read_object(path) or {}, path)
            if row is None:
                continue
            grouped[(str(row["backend"]), str(row["parameter_label"]), str(row["run_id"]))].append(row)
            raw_counts[str(row["status"])] += 1

    configurations: list[dict[str, Any]] = []
    for (backend, label, run_id), rows in sorted(grouped.items()):
        ok = [row for row in rows if row["status"] == "ok"]
        skipped = [row for row in rows if row["status"] == "skipped"]
        errors = [row for row in rows if row["status"] == "error"]
        if ok and len(ok) >= EXPECTED_REPEATS:
            status = "complete"
        elif ok:
            status = "partial"
        elif errors:
            status = "error"
        else:
            status = "skipped"
        first = rows[0]
        representative = next((row for row in rows if row["status"] == "ok"), first)
        unavailable = [row for row in rows if row["status"] in {"error", "skipped"}]
        configurations.append(
            {
                "backend": backend,
                "label": label,
                "run_id": run_id,
                "status": status,
                "ok_repeats": len(ok),
                "skipped_repeats": len(skipped),
                "error_repeats": len(errors),
                "recall": _median(row["recall"] for row in ok),
                "qps": _median(row["qps"] for row in ok),
                "index_bytes": _median(row["index_bytes"] for row in ok),
                "build_seconds": _median(row["build_seconds"] for row in ok),
                "build_rss": _median(row["build_rss"] for row in ok),
                "build_vram": _median(row["build_vram"] for row in ok),
                "search_rss": _median(row["search_rss"] for row in ok),
                "search_vram": _median(row["search_vram"] for row in ok),
                "capacity": next((row["capacity"] for row in unavailable if row.get("capacity")), None),
                "parameters": representative["parameters"],
                "provenance": representative["provenance"],
                "measurement": representative["measurement"],
                "metric": representative["metric"],
                "k": representative["k"],
                "warmup_queries": representative["warmup_queries"],
                "timed_queries": representative["timed_queries"],
                "latency": representative["latency"],
                "is_gpu": any(row["is_gpu"] for row in rows),
                "repeats": [
                    {
                        "status": row["status"],
                        "repeat": row["repeat"],
                        "recall_at_100": row["recall"],
                        "qps": row["qps"],
                        "latency_ms": row["latency"],
                        "index_bytes": row["index_bytes"],
                        "build_seconds": row["build_seconds"],
                        "build_peak_rss_bytes": row["build_rss"],
                        "build_peak_vram_bytes": row["build_vram"],
                        "search_peak_rss_bytes": row["search_rss"],
                        "search_peak_vram_bytes": row["search_vram"],
                        "error": row["error"],
                        "source": row["source"],
                    }
                    for row in rows
                ],
                "raw_database_bytes": (_bytes(first["dataset_window_count"]) or 0)
                * (_bytes(first["dimension"]) or 0)
                * 4,
                "dataset_window_count": first["dataset_window_count"],
                "dimension": first["dimension"],
                "error": next((row["error"] for row in unavailable if row.get("error")), None),
            }
        )

    cache = _cache_state(dataset_dir)
    return {
        "dataset_id": dataset_dir.name,
        "cache": cache,
        "configurations": configurations,
        "record_counts": dict(raw_counts),
    }


def _status_class(status: str) -> str:
    return {"complete": "ok", "partial": "partial", "error": "error", "skipped": "skip"}.get(status, "partial")


def _tooltip_info(dataset: dict[str, Any], row: dict[str, Any]) -> dict[str, Any]:
    """Return the complete configuration/repeat metadata shown on hover."""

    return {
        "dataset": dataset["dataset_id"],
        "backend": row["backend"],
        "GPU": "yes" if row.get("is_gpu") else "no",
        "configuration": row["label"],
        "status": row["status"],
        "run_id": row["run_id"],
        "repeat_summary": f"{row['ok_repeats']}/{EXPECTED_REPEATS} successful; {row['skipped_repeats']} skipped; {row['error_repeats']} errors",
        "indexed_windows": row.get("dataset_window_count"),
        "vector_dimension": row.get("dimension"),
        "cache": dataset.get("cache"),
        "metric": row.get("metric"),
        "k": row.get("k"),
        "recall_at_100_median": row.get("recall"),
        "queries_per_second_median": row.get("qps"),
        "raw_database": _fmt_bytes(row.get("raw_database_bytes")),
        "raw_database_bytes": row.get("raw_database_bytes"),
        "serialized_index": _fmt_bytes(row.get("index_bytes")),
        "serialized_index_bytes": row.get("index_bytes"),
        "build_seconds": row.get("build_seconds"),
        "build_peak_rss": _fmt_bytes(row.get("build_rss")),
        "build_peak_vram": _fmt_bytes(row.get("build_vram")),
        "search_peak_rss": _fmt_bytes(row.get("search_rss")),
        "search_peak_vram": _fmt_bytes(row.get("search_vram")),
        "warmup_queries": row.get("warmup_queries"),
        "timed_queries": row.get("timed_queries"),
        "latency_ms": row.get("latency"),
        "parameters": row.get("parameters"),
        "capacity_preflight": row.get("capacity"),
        "measurement": row.get("measurement"),
        "provenance": row.get("provenance"),
        "repeats": row.get("repeats"),
        "error": row.get("error"),
    }


def _chart(dataset: dict[str, Any]) -> str:
    points = [row for row in dataset["configurations"] if row.get("qps") and row.get("recall") is not None]
    width, height = 860, 350
    left, top, plot_width, plot_height = 68, 24, 585, 250
    right, bottom = left + plot_width, top + plot_height
    if not points:
        return '<div class="empty-chart">No successful search points yet.</div>'
    qps_values = [float(row["qps"]) for row in points if float(row["qps"]) > 0]
    low, high = min(qps_values), max(qps_values)
    log_low = math.floor(math.log10(low))
    log_high = math.ceil(math.log10(high))
    if log_high <= log_low:
        log_high = log_low + 1

    def x(value: float) -> float:
        return left + max(0.0, min(1.0, float(value))) * plot_width

    def y(value: float) -> float:
        return bottom - (math.log10(max(float(value), 10**log_low)) - log_low) / (log_high - log_low) * plot_height

    lines = [
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" aria-label="Queries per second versus recall for {_escape(dataset["dataset_id"])}">',
        f'<line class="axis" x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" />',
        f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" />',
        f'<text x="{left + plot_width / 2}" y="{bottom + 38}" text-anchor="middle">recall@100</text>',
        f'<text x="16" y="{top + plot_height / 2}" transform="rotate(-90 16 {top + plot_height / 2})" text-anchor="middle">queries/s (log scale)</text>',
    ]
    for tick in (0.0, 0.25, 0.5, 0.75, 1.0):
        xx = x(tick)
        lines.append(f'<line class="grid" x1="{xx:.1f}" y1="{top}" x2="{xx:.1f}" y2="{bottom}" />')
        lines.append(f'<text x="{xx:.1f}" y="{bottom + 18}" text-anchor="middle">{tick:.2f}</text>')
    for exponent in range(log_low, log_high + 1):
        tick = 10**exponent
        yy = y(tick)
        lines.append(f'<line class="grid" x1="{left}" y1="{yy:.1f}" x2="{right}" y2="{yy:.1f}" />')
        lines.append(f'<text x="{left - 8}" y="{yy + 4:.1f}" text-anchor="end">{_fmt_number(tick, 0)}</text>')
    for row in points:
        xx, yy = x(float(row["recall"])), y(float(row["qps"]))
        color = COLORS.get(str(row["backend"]), "#777")
        title = _escape(f"{row['backend']} / {row['label']} — recall {row['recall']:.4f}, QPS {row['qps']:.1f}")
        point_class = "gpu-point" if row.get("is_gpu") else "cpu-point"
        info = json.dumps(_tooltip_info(dataset, row), sort_keys=True, separators=(",", ":"), default=str)
        lines.append(
            f'<g class="chart-point {point_class}" tabindex="0" data-info="{_escape(info)}">'
            f'<circle cx="{xx:.1f}" cy="{yy:.1f}" r="5" fill="{color}" stroke="#17202a" stroke-width="1">'
            f'<title>{title}</title></circle></g>'
        )
    legend_x, legend_y = right + 28, top + 8
    lines.append(f'<text x="{legend_x}" y="{legend_y}" class="legend-title">Backend</text>')
    for index, backend in enumerate(BACKENDS):
        if not any(str(row["backend"]) == backend for row in points):
            continue
        yy = legend_y + 22 + index * 22
        color = COLORS[backend]
        lines.append(f'<circle cx="{legend_x + 5}" cy="{yy - 4}" r="5" fill="{color}" stroke="#17202a" stroke-width="1" />')
        lines.append(f'<text x="{legend_x + 16}" y="{yy}" class="legend-label">{_escape(backend.upper())}</text>')
    lines.append("</svg>")
    return "".join(lines)


def _configuration_table(configurations: list[dict[str, Any]]) -> str:
    rows: list[str] = []
    for row in configurations:
        status = _status_class(str(row["status"]))
        resource = " / ".join(
            [
                f"B RSS {_fmt_bytes(row['build_rss'])}",
                f"S RSS {_fmt_bytes(row['search_rss'])}",
                f"B VRAM {_fmt_bytes(row['build_vram'])}",
                f"S VRAM {_fmt_bytes(row['search_vram'])}",
            ]
        )
        capacity = row.get("capacity") or {}
        gate = " / ".join(
            part
            for part in (
                f"gate RAM {_fmt_bytes(capacity.get('estimated_ram_bytes'))}"
                if capacity.get("estimated_ram_bytes") is not None
                else None,
                f"gate VRAM {_fmt_bytes(capacity.get('estimated_vram_bytes'))}"
                if capacity.get("estimated_vram_bytes") is not None
                else None,
                f"free RAM {_fmt_bytes(capacity.get('available_ram_bytes'))}"
                if capacity.get("available_ram_bytes") is not None
                else None,
                f"free VRAM {_fmt_bytes(capacity.get('available_vram_bytes'))}"
                if capacity.get("available_vram_bytes") is not None
                else None,
            )
            if part
        )
        if gate:
            resource = f"{resource} / {gate}"
        message = row.get("error") or ""
        rows.append(
            "<tr>"
            f'<td><span class="pill {status}">{_escape(row["status"])}</span></td>'
            f'<td>{_escape(row["backend"])}</td>'
            f'<td class="mono">{_escape(row["label"])}</td>'
            f'<td>{row["ok_repeats"]}/{EXPECTED_REPEATS} ok</td>'
            f'<td>{_fmt_number(row["recall"], 4)}</td>'
            f'<td>{_fmt_number(row["qps"])}</td>'
            f'<td>{_fmt_bytes(row["raw_database_bytes"])}</td>'
            f'<td>{_fmt_bytes(row["index_bytes"])}</td>'
            f'<td>{_escape(resource)}</td>'
            f'<td>{_escape(message[:180])}</td>'
            "</tr>"
        )
    return (
        '<div class="table-wrap"><table><thead><tr>'
        "<th>Status</th><th>Backend</th><th>Configuration</th><th>Repeats</th>"
        "<th>Recall</th><th>QPS</th><th>Raw DB</th><th>Index</th>"
        "<th>Memory (build / search)</th><th>Reason</th>"
        "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></div>"
    )


def _dataset_html(dataset: dict[str, Any]) -> str:
    cache = dataset["cache"]
    configurations = dataset["configurations"]
    counts = Counter(row["status"] for row in configurations)
    summary = " ".join(
        f'<span class="pill {_status_class(status)}">{_escape(status)} {count}</span>'
        for status, count in sorted(counts.items())
    ) or '<span class="muted">no result rows</span>'
    cache_line = (
        f"request <b>{_escape(cache['request_state'])}</b>"
        + (f" · { _escape(cache['request_id'])}" if cache.get("request_id") else "")
        + f" · windows {cache.get('window_count') or '—'} · queries {cache.get('query_count') or '—'}"
        + (f" · ground truth { _escape(cache['ground_truth_engine'])}" if cache.get("ground_truth_engine") else "")
    )
    return (
        f'<section class="dataset"><div class="dataset-head"><div><h2>{_escape(dataset["dataset_id"])}</h2><p>{cache_line}</p></div><div>{summary}</div></div>'
        f'<div class="chart-wrap"><h3>Recall@100 vs queries/s (all measured points)</h3>{_chart(dataset)}</div>'
        f'<h3>Configuration ledger</h3>{_configuration_table(configurations)}</section>'
    )


def _running_html(jobs: list[dict[str, Any]]) -> str:
    if not jobs:
        return '<section class="running"><h2>Running benchmarks</h2><p class="muted">No active FAISS, ScaNN, NGT, or cuVS runner container detected.</p></section>'
    cards: list[str] = []
    for job in jobs:
        progress = job["progress"]
        estimate = progress.get("progress")
        if estimate is None:
            progress_text = f"{progress['terminal_configurations']} terminal configurations · {progress['record_count']} result records observed"
            progress_width = 0
        else:
            progress_text = (
                f"{progress['terminal_configurations']}/{progress['expected_configurations']} terminal configurations "
                f"({estimate * 100:.0f}% estimated) · {progress['record_count']} result records observed"
            )
            progress_width = estimate * 100
        resources = " · ".join(
            item
            for item in (
                f"elapsed {_fmt_duration(job['elapsed_seconds'])}" if job.get("elapsed_seconds") is not None else None,
                f"CPU {_fmt_number(job['cpu'], 1)}%" if job.get("cpu") is not None else None,
                f"RSS {_fmt_bytes(job['rss_bytes'])}" if job.get("rss_bytes") is not None else None,
            )
            if item
        )
        cards.append(
            '<article class="job-card">'
            f'<div class="job-head"><h3>{_escape(job["backend"].upper())} · {_escape(job["dataset"])}</h3><span class="pill partial">{_escape(job["status"])}</span></div>'
            f'<p><b>Stage:</b> <code>{_escape(job["stage"])}</code></p>'
            f'<p>{_escape(resources)} · container <code>{_escape(job["container"])}</code></p>'
            f'<div class="progress-track" role="progressbar" aria-label="Estimated benchmark progress" aria-valuemin="0" aria-valuemax="100" aria-valuenow="{progress_width:.0f}"><div class="progress-bar" style="width:{progress_width:.1f}%"></div></div>'
            f'<p class="progress-text">{_escape(progress_text)}. Progress is estimated from terminal result JSON files; a configuration is terminal after all repeats, a skip, or an error.</p>'
            f'<p class="muted"><b>Image:</b> {_escape(job["image"])}</p>'
            '</article>'
        )
    return '<section class="running"><h2>Running benchmarks</h2><p class="muted">Active runner containers and best-effort process telemetry.</p>' + "".join(cards) + "</section>"


def render(cache_dir: Path, output: Path, *, refresh_seconds: int) -> None:
    datasets = [_collect_dataset(path) for path in sorted(cache_dir.iterdir()) if path.is_dir() and (path / "results").is_dir()]
    running_jobs = _running_jobs(cache_dir)
    now = datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")
    total_configs = sum(len(dataset["configurations"]) for dataset in datasets)
    total_records = sum(sum(dataset["record_counts"].values()) for dataset in datasets)
    dataset_sections = "".join(_dataset_html(dataset) for dataset in datasets)
    running_section = _running_html(running_jobs)
    html_text = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta http-equiv="refresh" content="{max(5, int(refresh_seconds))}"><meta http-equiv="Cache-Control" content="no-store, no-cache, must-revalidate">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>GINflow benchmark dashboard</title>
<style>
:root {{ color-scheme: dark; --bg:#10151c; --panel:#18212b; --line:#2b3948; --text:#e8eef5; --muted:#9dafc1; --ok:#39c88f; --warn:#e5b34e; --bad:#ef6c6c; --blue:#5ea5d6; }}
* {{ box-sizing:border-box; }} body {{ margin:0; background:var(--bg); color:var(--text); font:14px/1.45 Inter,system-ui,sans-serif; }}
main {{ max-width:1800px; margin:0 auto; padding:24px; }} h1 {{ margin:0 0 6px; font-size:28px; }} h2 {{ margin:0; font-size:22px; }} h3 {{ margin:18px 0 10px; font-size:15px; color:#cbd8e5; }} p {{ color:var(--muted); margin:4px 0; }}
.top {{ display:flex; justify-content:space-between; gap:16px; align-items:flex-end; margin-bottom:18px; }} .meta {{ text-align:right; color:var(--muted); }}
.notice {{ background:#1e2b38; border:1px solid #35506a; border-radius:10px; padding:12px 14px; margin:14px 0 20px; }} code,.mono {{ font-family:"SFMono-Regular",Consolas,monospace; }}
.dataset {{ background:var(--panel); border:1px solid var(--line); border-radius:12px; padding:18px; margin:18px 0 24px; box-shadow:0 8px 28px #0002; }} .dataset-head {{ display:flex; justify-content:space-between; gap:16px; align-items:flex-start; }}
.running {{ background:#172936; border:1px solid #3d6a82; border-radius:12px; padding:16px 18px; margin:0 0 24px; }} .running h2 {{ margin:0 0 4px; }} .job-card {{ background:#11202a; border:1px solid #315064; border-radius:9px; padding:12px 14px; margin:12px 0 0; }} .job-head {{ display:flex; justify-content:space-between; align-items:center; gap:12px; }} .job-card h3 {{ margin:0; color:#e7f2fa; }} .job-card p {{ margin:6px 0; }} .progress-track {{ height:10px; background:#0b1319; border:1px solid #38566a; border-radius:999px; overflow:hidden; margin-top:10px; }} .progress-bar {{ height:100%; background:linear-gradient(90deg,#3fa5d5,#39c88f); border-radius:999px; transition:width .3s ease; }} .progress-text {{ font-size:12px; color:#c2d5e2; }}
.pill {{ display:inline-block; border-radius:999px; padding:2px 8px; margin:2px 3px 2px 0; font-size:12px; border:1px solid var(--line); }} .pill.ok,.pill.complete {{ color:var(--ok); border-color:#287d61; }} .pill.partial {{ color:var(--warn); border-color:#84692a; }} .pill.error {{ color:var(--bad); border-color:#843f48; }} .pill.skip,.pill.skipped {{ color:#b8c5d1; border-color:#536271; }}
.controls {{ position:sticky; top:0; z-index:10; display:flex; align-items:center; gap:14px; background:#22384a; border:2px solid var(--blue); border-radius:9px; padding:11px 14px; margin:0 0 14px; color:#f1f7fc; box-shadow:0 4px 18px #0008; }} .controls label {{ cursor:pointer; font-size:16px; font-weight:700; white-space:nowrap; }} .controls input {{ accent-color:var(--blue); width:18px; height:18px; margin:0 8px 0 0; vertical-align:-3px; }} .controls .muted {{ font-size:12px; }}
.chart-wrap {{ background:#121a23; border:1px solid var(--line); border-radius:9px; padding:8px 12px 2px; margin-top:16px; }} .chart {{ width:100%; max-height:370px; }} .chart text {{ fill:#aebdcb; font-size:11px; }} .chart .legend-title {{ fill:#d8e2ec; font-weight:700; }} .chart .legend-label {{ fill:#c5d0db; font-size:12px; }} .chart .axis {{ stroke:#8294a6; }} .chart .grid {{ stroke:#2b3948; stroke-dasharray:3 4; }} .chart-point {{ cursor:crosshair; outline:none; transition:opacity .15s ease, filter .15s ease; }} .chart-point:focus circle {{ stroke:#ffffff; stroke-width:2; }} .pareto-line {{ transition:opacity .15s ease; }} body.gpu-highlight .cpu-point {{ opacity:.16; filter:saturate(.15); }} body.gpu-highlight .gpu-point {{ opacity:1; filter:none; }} .empty-chart {{ color:var(--muted); padding:30px; }}
.chart-tooltip {{ position:fixed; z-index:20; pointer-events:none; max-width:560px; max-height:78vh; overflow:auto; background:#0b1118f5; border:1px solid #55718a; border-radius:8px; box-shadow:0 10px 28px #000b; padding:9px 11px; color:var(--text); font-size:12px; }} .chart-tooltip[hidden] {{ display:none; }} .tooltip-row {{ display:grid; grid-template-columns:minmax(150px, max-content) minmax(180px, 1fr); gap:9px; border-bottom:1px solid #263746; padding:3px 0; }} .tooltip-row:last-child {{ border-bottom:0; }} .tooltip-row b {{ color:#b9d7ee; }} .tooltip-row pre {{ margin:0; white-space:pre-wrap; overflow-wrap:anywhere; font:12px/1.35 "SFMono-Regular",Consolas,monospace; color:#e3ebf2; }}
.table-wrap {{ overflow:auto; border:1px solid var(--line); border-radius:8px; }} table {{ border-collapse:collapse; width:100%; min-width:1250px; }} th,td {{ border-bottom:1px solid var(--line); padding:7px 8px; text-align:left; vertical-align:top; white-space:nowrap; }} th {{ position:sticky; top:0; background:#202c38; color:#cbd8e5; font-weight:600; }} tr:hover td {{ background:#202a35; }} td:last-child {{ white-space:normal; max-width:420px; color:#d5a6a6; }} .muted {{ color:var(--muted); }} footer {{ color:var(--muted); margin:22px 0 8px; }}
@media (max-width:800px) {{ main {{ padding:12px; }} .top,.dataset-head {{ display:block; }} .meta {{ text-align:left; margin-top:10px; }} }}
</style></head><body><main>
<div class="top"><div><h1>GINflow vector-index benchmark dashboard</h1><p>Live external-cache view · full measured points retained · refreshes every {max(5, int(refresh_seconds))} seconds</p></div><div class="meta">Updated {_escape(now)}<br>{total_configs} configurations · {total_records} result records</div></div>
<div class="notice">Raw results: <code>{_escape(cache_dir)}</code>. Build memory and search memory are shown separately. Capacity-gate estimates are retained for skipped/error rows. Search peaks are currently warm in-process measurements; VRAM is device-wide <code>nvidia-smi</code> usage where applicable. To serve this file: <code>python3 -m http.server --directory {output.parent} 8765</code>.</div>
<div class="controls"><label><input id="gpu-highlight" type="checkbox"> HIGHLIGHT GPU INDEXES</label><label><input id="pareto-toggle" type="checkbox"> SHOW PARETO CURVES</label><span class="muted">Hover or focus a point for complete configuration, repeat, resource, hardware, and provenance details.</span></div>
{running_section}
{dataset_sections or '<div class="notice">No benchmark dataset cache has been discovered yet.</div>'}
<footer>Generated by <code>benchmarks/live_dashboard.py</code>. The HTML is replaced atomically; keep this process running for live updates.</footer>
</main><div id="chart-tooltip" class="chart-tooltip" role="tooltip" hidden></div>{INTERACTION_SCRIPT}</body></html>
'''
    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{output.name}.", suffix=".tmp", dir=output.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(html_text)
        os.replace(temporary, output)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, required=True, help="external benchmark cache root")
    parser.add_argument("--output", type=Path, required=True, help="HTML output path")
    parser.add_argument("--refresh-seconds", type=int, default=15, help="browser refresh and watch interval")
    parser.add_argument("--watch", action="store_true", help="rewrite the dashboard until interrupted")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.refresh_seconds < 5:
        raise SystemExit("--refresh-seconds must be >= 5")
    while True:
        render(args.cache_dir.resolve(), args.output.resolve(), refresh_seconds=args.refresh_seconds)
        if not args.watch:
            return 0
        try:
            time.sleep(args.refresh_seconds)
        except KeyboardInterrupt:
            return 0


if __name__ == "__main__":
    raise SystemExit(main())
