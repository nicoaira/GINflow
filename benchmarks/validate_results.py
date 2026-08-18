#!/usr/bin/env python3
"""Validate GINflow vector-index benchmark result records.

The validator is intentionally strict about successful measurements: a recall
versus throughput point is useful only when its query set, exact ground truth,
hardware scope, warm-up, timed work, and runner provenance are known.  Failed
or skipped configurations remain in the output so a report can state what was
not benchmarked without accidentally plotting them.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from statistics import median
from typing import Any, Iterable

from results_common import (
    COMPARISON_SCOPE_FIELDS,
    LocatedRecord,
    as_int,
    as_number,
    canonical_json,
    comparison_signature,
    configuration_signature,
    hardware_id,
    measurement_value,
    provenance_value,
    read_records,
    stable_id,
)


SUPPORTED_BACKENDS = {"faiss", "scann", "ngt", "cuvs"}
EXPECTED_METRIC = "cosine"
EXPECTED_K = 100
PROVENANCE_FIELDS = (
    "git_commit",
    "runner_version",
    "hardware_id",
    "embedding_cache_id",
    "query_selection_id",
    "ground_truth_cache_id",
)
MEASUREMENT_PROTOCOL = "fixed-query-batches-v1"
LATENCY_UNIT = "milliseconds_per_query_batch"
QPS_SCOPE = "queries_per_second_over_timed_batches"


@dataclass(frozen=True)
class Issue:
    severity: str
    code: str
    message: str
    source: str | None = None

    def as_dict(self) -> dict[str, str]:
        result = {"severity": self.severity, "code": self.code, "message": self.message}
        if self.source:
            result["source"] = self.source
        return result


def _is_nonempty_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _is_sha256(value: Any) -> bool:
    if not _is_nonempty_string(value) or len(value) != 64:
        return False
    return all(character in "0123456789abcdef" for character in value)


def _version_is_v1(value: Any) -> bool:
    if not _is_nonempty_string(value):
        return False
    normalized = value.strip().lower()
    return normalized in {"1", "v1", "ginflow-benchmark/v1", "ginflow-benchmark-v1"}


def _number_in_range(value: Any, lower: float, upper: float | None = None) -> bool:
    number = as_number(value)
    return number is not None and number >= lower and (upper is None or number <= upper)


def _integer_in_range(value: Any, lower: int, upper: int | None = None) -> bool:
    number = as_int(value)
    return number is not None and number >= lower and (upper is None or number <= upper)


def _field_issue(issues: list[Issue], source: str, field: str, expectation: str) -> None:
    issues.append(Issue("error", "invalid_field", f"{field} must be {expectation}", source))


def _validate_identity(payload: dict[str, Any], source: str, issues: list[Issue]) -> None:
    """Validate fields that identify every result, including unavailable runs."""

    if not _version_is_v1(payload.get("schema_version")):
        _field_issue(issues, source, "schema_version", "a supported v1 identifier")
    if payload.get("backend") not in SUPPORTED_BACKENDS:
        _field_issue(issues, source, "backend", f"one of {sorted(SUPPORTED_BACKENDS)}")
    for field in ("dataset_id", "parameter_label", "run_id"):
        if not _is_nonempty_string(payload.get(field)):
            _field_issue(issues, source, field, "a non-empty string")
    if not _integer_in_range(payload.get("dataset_window_count"), 1):
        _field_issue(issues, source, "dataset_window_count", "a positive integer")
    if not _integer_in_range(payload.get("dimension"), 1):
        _field_issue(issues, source, "dimension", "a positive integer")
    if payload.get("metric") != EXPECTED_METRIC:
        _field_issue(issues, source, "metric", repr(EXPECTED_METRIC))
    if as_int(payload.get("k")) != EXPECTED_K:
        _field_issue(issues, source, "k", str(EXPECTED_K))
    if not isinstance(payload.get("parameters"), dict):
        _field_issue(issues, source, "parameters", "an object")
    if not _integer_in_range(payload.get("repeat"), 0):
        _field_issue(issues, source, "repeat", "a non-negative integer")


def _validate_success_measurement(payload: dict[str, Any], source: str, issues: list[Issue]) -> None:
    if not _integer_in_range(payload.get("warmup_queries"), 1):
        _field_issue(issues, source, "warmup_queries", "a positive integer")
    if not _integer_in_range(payload.get("timed_queries"), 1):
        _field_issue(issues, source, "timed_queries", "a positive integer")
    warmup = as_int(payload.get("warmup_queries"))
    timed = as_int(payload.get("timed_queries"))
    if warmup is not None and timed is not None and warmup > timed:
        issues.append(Issue("error", "invalid_measurement", "warmup_queries cannot exceed timed_queries", source))
    if not _number_in_range(payload.get("qps"), 0) or as_number(payload.get("qps")) == 0:
        _field_issue(issues, source, "qps", "a positive finite number")
    if not _number_in_range(payload.get("recall_at_100"), 0, 1):
        _field_issue(issues, source, "recall_at_100", "a finite number in [0, 1]")
    for field in ("query_ids_sha256", "ground_truth_ids_sha256"):
        if not _is_sha256(payload.get(field)):
            _field_issue(issues, source, field, "a 64-character SHA-256 hexadecimal digest")
    if not _integer_in_range(payload.get("index_bytes"), 0):
        _field_issue(issues, source, "index_bytes", "a non-negative integer")
    if not _number_in_range(payload.get("build_seconds"), 0):
        _field_issue(issues, source, "build_seconds", "a non-negative finite number")
    for field in ("peak_rss_bytes", "peak_vram_bytes"):
        value = payload.get(field)
        if value is not None and not _integer_in_range(value, 0):
            _field_issue(issues, source, field, "null or a non-negative integer")
    latency = payload.get("latency_ms")
    if not isinstance(latency, dict):
        _field_issue(issues, source, "latency_ms", "an object with mean, p50, and p95")
    else:
        for field in ("mean", "p50", "p95"):
            if not _number_in_range(latency.get(field), 0) or as_number(latency.get(field)) == 0:
                _field_issue(issues, source, f"latency_ms.{field}", "a positive finite number")
        p50 = as_number(latency.get("p50"))
        p95 = as_number(latency.get("p95"))
        if p50 is not None and p95 is not None and p50 > p95:
            issues.append(Issue("error", "invalid_latency_percentiles", "latency_ms.p50 cannot exceed latency_ms.p95", source))
    measurement = payload.get("measurement")
    if not isinstance(measurement, dict):
        _field_issue(issues, source, "measurement", "an object with query_batch_size")
    else:
        batch_size = measurement_value(payload, "query_batch_size")
        timed_batch_count = measurement_value(payload, "timed_batch_count")
        measured_timed_queries = measurement_value(payload, "timed_queries")
        if measurement.get("protocol") != MEASUREMENT_PROTOCOL:
            _field_issue(issues, source, "measurement.protocol", repr(MEASUREMENT_PROTOCOL))
        if not _integer_in_range(batch_size, 1):
            _field_issue(issues, source, "measurement.query_batch_size", "a positive integer")
        if not _integer_in_range(timed_batch_count, 1):
            _field_issue(issues, source, "measurement.timed_batch_count", "a positive integer")
        if not _integer_in_range(measured_timed_queries, 1):
            _field_issue(issues, source, "measurement.timed_queries", "a positive integer")
        if measurement.get("latency_unit") != LATENCY_UNIT:
            _field_issue(issues, source, "measurement.latency_unit", repr(LATENCY_UNIT))
        if measurement.get("qps_scope") != QPS_SCOPE:
            _field_issue(issues, source, "measurement.qps_scope", repr(QPS_SCOPE))
        if timed is not None and as_int(measured_timed_queries) is not None and as_int(measured_timed_queries) != timed:
            issues.append(
                Issue(
                    "error",
                    "invalid_measurement",
                    "measurement.timed_queries must equal top-level timed_queries",
                    source,
                )
            )
        if (
            as_int(batch_size) is not None
            and as_int(timed_batch_count) is not None
            and as_int(measured_timed_queries) is not None
            and as_int(batch_size) * as_int(timed_batch_count) != as_int(measured_timed_queries)
        ):
            issues.append(
                Issue(
                    "error",
                    "invalid_measurement",
                    "measurement.query_batch_size * measurement.timed_batch_count must equal measurement.timed_queries",
                    source,
                )
            )


def _validate_success_provenance(payload: dict[str, Any], source: str, issues: list[Issue]) -> None:
    provenance = payload.get("provenance")
    if not isinstance(provenance, dict):
        _field_issue(issues, source, "provenance", "an object")
    else:
        for field in PROVENANCE_FIELDS:
            value = hardware_id(payload) if field == "hardware_id" else provenance_value(payload, field)
            if not _is_nonempty_string(value):
                _field_issue(issues, source, f"provenance.{field}", "a non-empty string")


def _validate_success(record: LocatedRecord, issues: list[Issue]) -> None:
    """Validate identity, measurements, and provenance independently for clarity."""

    _validate_identity(record.payload, record.source, issues)
    _validate_success_measurement(record.payload, record.source, issues)
    _validate_success_provenance(record.payload, record.source, issues)


def _validate_non_success(record: LocatedRecord, issues: list[Issue]) -> None:
    """Require enough identity to document unavailable configurations honestly."""

    payload = record.payload
    source = record.source
    _validate_identity(payload, source, issues)
    reason = payload.get("error")
    if not _is_nonempty_string(reason):
        issues.append(Issue("error", "missing_unavailable_reason", "skipped/error record must have an explanatory error/reason", source))


def _scope_dict(signature: tuple[Any, ...]) -> dict[str, Any]:
    result = dict(zip(COMPARISON_SCOPE_FIELDS, signature, strict=True))
    result["scope_id"] = stable_id(signature)
    return result


def aggregate_successful_records(records: Iterable[LocatedRecord]) -> list[dict[str, Any]]:
    """Aggregate repeats, preserving min/max dispersion and provenance references."""

    grouped: dict[tuple[Any, ...], list[LocatedRecord]] = defaultdict(list)
    for record in records:
        grouped[configuration_signature(record.payload)].append(record)

    aggregates: list[dict[str, Any]] = []
    for key in sorted(grouped, key=canonical_json):
        group = sorted(grouped[key], key=lambda item: (as_int(item.payload.get("repeat")) or -1, item.source))
        first = group[0].payload
        qps_values = [as_number(item.payload["qps"]) for item in group]
        recall_values = [as_number(item.payload["recall_at_100"]) for item in group]
        assert all(value is not None for value in qps_values + recall_values)
        qps = [float(value) for value in qps_values if value is not None]
        recalls = [float(value) for value in recall_values if value is not None]
        latency_medians = {
            statistic: float(
                median(
                    float(value)
                    for value in (as_number(item.payload["latency_ms"].get(statistic)) for item in group)
                    if value is not None
                )
            )
            for statistic in ("mean", "p50", "p95")
        }
        optional_medians: dict[str, float | None] = {}
        for field in ("index_bytes", "build_seconds", "peak_rss_bytes", "peak_vram_bytes"):
            values = [as_number(item.payload.get(field)) for item in group]
            usable = [float(value) for value in values if value is not None]
            optional_medians[field] = float(median(usable)) if usable else None
        scope = _scope_dict(comparison_signature(first))
        aggregates.append(
            {
                **scope,
                "backend": first["backend"],
                "parameter_label": first["parameter_label"],
                "parameters": first["parameters"],
                "run_id": first["run_id"],
                "repeat_count": len(group),
                "repeat_indices": [as_int(item.payload["repeat"]) for item in group],
                "qps_median": float(median(qps)),
                "qps_min": min(qps),
                "qps_max": max(qps),
                "recall_median": float(median(recalls)),
                "recall_min": min(recalls),
                "recall_max": max(recalls),
                "latency_ms_median": latency_medians,
                "measurement": first["measurement"],
                **optional_medians,
                "sources": [item.source for item in group],
            }
        )
    return aggregates


def _partition_records(
    records: Iterable[LocatedRecord], issues: list[Issue]
) -> tuple[list[LocatedRecord], list[dict[str, Any]], Counter[str]]:
    """Validate individual rows and split usable measurements from unavailable ones."""

    valid_successes: list[LocatedRecord] = []
    unavailable: list[dict[str, Any]] = []
    status_counts: Counter[str] = Counter()
    for record in records:
        payload = record.payload
        status = payload.get("status")
        if status not in {"ok", "skipped", "error"}:
            issues.append(Issue("error", "invalid_status", "status must be one of ok, skipped, error", record.source))
            continue
        status_counts[str(status)] += 1
        before = len(issues)
        if status == "ok":
            _validate_success(record, issues)
            if not any(issue.severity == "error" for issue in issues[before:]):
                valid_successes.append(record)
            continue
        _validate_non_success(record, issues)
        unavailable.append(
            {
                "status": status,
                "backend": payload.get("backend"),
                "dataset_id": payload.get("dataset_id"),
                "parameter_label": payload.get("parameter_label"),
                "run_id": payload.get("run_id"),
                "error": payload.get("error"),
                "source": record.source,
            }
        )
    return valid_successes, unavailable, status_counts


def _validate_repeat_groups(records: Iterable[LocatedRecord], issues: list[Issue], min_repeats: int) -> None:
    """Reject duplicate repeat IDs and configurations with too little replication."""

    groups: dict[tuple[Any, ...], list[LocatedRecord]] = defaultdict(list)
    for record in records:
        groups[configuration_signature(record.payload)].append(record)
    for group in groups.values():
        repeats = [as_int(record.payload.get("repeat")) for record in group]
        if len(set(repeats)) != len(repeats):
            sources = ", ".join(record.source for record in group)
            issues.append(Issue("error", "duplicate_repeat", f"duplicate repeat index for one configuration: {sources}"))
        if len(group) < min_repeats:
            issues.append(
                Issue(
                    "error",
                    "insufficient_repeats",
                    f"configuration has {len(group)} valid repeat(s); at least {min_repeats} are required",
                    group[0].source,
                )
            )


def _validate_comparison_scopes(records: Iterable[LocatedRecord], issues: list[Issue]) -> None:
    """Ensure a dataset/hardware scope represents one database and exact truth set."""

    dataset_hardware_identities: dict[tuple[Any, Any], set[tuple[Any, ...]]] = defaultdict(set)
    for record in records:
        payload = record.payload
        identity = comparison_signature(payload)
        dataset_hardware_identities[(payload.get("dataset_id"), hardware_id(payload))].add(identity)
    for (dataset_id, machine), identities in sorted(dataset_hardware_identities.items(), key=lambda item: canonical_json(item[0])):
        if len(identities) > 1:
            issues.append(
                Issue(
                    "error",
                    "incomparable_scope",
                    f"dataset_id={dataset_id!r}, hardware_id={machine!r} has {len(identities)} incompatible comparison identities",
                )
            )


def _scope_summaries(aggregates: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = ("scope_id", *COMPARISON_SCOPE_FIELDS)
    scopes: dict[str, dict[str, Any]] = {}
    for aggregate in aggregates:
        scopes.setdefault(aggregate["scope_id"], {key: aggregate[key] for key in keys})
    return [scopes[key] for key in sorted(scopes)]


def validate_records(
    records: list[LocatedRecord],
    *,
    read_issues: Iterable[Any] = (),
    min_repeats: int = 2,
) -> dict[str, Any]:
    """Validate rows and return a machine-readable, auditable validation report."""

    issues: list[Issue] = [
        Issue("error", "input_read_error", issue.message, issue.source) for issue in read_issues
    ]
    if not records:
        issues.append(Issue("error", "no_records", "no benchmark result records were found"))
    valid_successes, unavailable, status_counts = _partition_records(records, issues)
    _validate_repeat_groups(valid_successes, issues, min_repeats)
    _validate_comparison_scopes(valid_successes, issues)
    aggregates = aggregate_successful_records(valid_successes)
    errors = [issue.as_dict() for issue in issues if issue.severity == "error"]
    warnings = [issue.as_dict() for issue in issues if issue.severity == "warning"]
    return {
        "schema_version": "ginflow-benchmark-validation/v1",
        "generated_at": datetime.now(UTC).isoformat(),
        "valid": not errors,
        "status_counts": dict(sorted(status_counts.items())),
        "records_found": len(records),
        "valid_successful_records": len(valid_successes),
        "successful_configurations": len(aggregates),
        "comparison_scopes": _scope_summaries(aggregates),
        "aggregates": aggregates,
        "unavailable": sorted(
            unavailable,
            key=lambda item: tuple(str(item.get(field) or "") for field in ("dataset_id", "backend", "parameter_label", "run_id", "status", "source")),
        ),
        "errors": errors,
        "warnings": warnings,
        "valid_success_sources": sorted(record.source for record in valid_successes),
    }


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", help="result files or directories to scan recursively")
    parser.add_argument("--output", type=Path, help="write the validation report as JSON")
    parser.add_argument(
        "--min-repeats",
        type=int,
        default=2,
        help="minimum valid timed repeat records required per configuration (default: 2)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.min_repeats < 1:
        raise SystemExit("--min-repeats must be at least 1")
    records, read_issues, ignored = read_records(args.inputs)
    report = validate_records(records, read_issues=read_issues, min_repeats=args.min_repeats)
    report["inputs"] = [str(Path(path)) for path in args.inputs]
    report["ignored_files"] = ignored
    if args.output:
        _write_json(args.output, report)
    print(
        "validation: "
        f"{report['valid_successful_records']} successful records, "
        f"{report['successful_configurations']} configurations, "
        f"{len(report['errors'])} error(s), {len(report['warnings'])} warning(s)"
    )
    for issue in report["errors"]:
        location = f" [{issue['source']}]" if "source" in issue else ""
        print(f"ERROR {issue['code']}{location}: {issue['message']}", file=sys.stderr)
    for issue in report["warnings"]:
        location = f" [{issue['source']}]" if "source" in issue else ""
        print(f"WARNING {issue['code']}{location}: {issue['message']}", file=sys.stderr)
    return 0 if report["valid"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
