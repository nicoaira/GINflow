#!/usr/bin/env python3
"""Shared parsing and normalization for GINflow benchmark result records.

The benchmark runners write one JSON object for each timed repeat.  Keeping the
reader in the standard library makes the validation and plotting tools usable
on a fresh checkout, before any index-specific environment is installed.
"""
from __future__ import annotations

import csv
import hashlib
import json
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any


RECORD_SUFFIXES = {".json", ".jsonl", ".ndjson", ".csv", ".tsv"}

# Keep this ordered list next to ``comparison_signature`` so the validator,
# report metadata, and plot grouping cannot silently disagree about what makes
# two recall/QPS points comparable.
COMPARISON_SCOPE_FIELDS = (
    "dataset_id",
    "dataset_window_count",
    "dimension",
    "metric",
    "k",
    "warmup_queries",
    "timed_queries",
    "query_ids_sha256",
    "ground_truth_ids_sha256",
    "embedding_cache_id",
    "query_selection_id",
    "ground_truth_cache_id",
    "hardware_id",
    "git_commit",
    "runner_version",
    "query_batch_size",
    "timed_batch_count",
)


@dataclass(frozen=True)
class LocatedRecord:
    """A decoded record together with a stable source location."""

    payload: dict[str, Any]
    source: str


@dataclass(frozen=True)
class ReadIssue:
    """A malformed input that could not be decoded into benchmark records."""

    source: str
    message: str


def canonical_json(value: Any) -> str:
    """Return deterministic JSON suitable for grouping configuration objects."""

    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def stable_id(parts: Iterable[Any], length: int = 12) -> str:
    """Make a compact, deterministic identifier from comparable fields."""

    digest = hashlib.sha256(canonical_json(list(parts)).encode("utf-8")).hexdigest()
    return digest[:length]


def get_nested(mapping: dict[str, Any], *paths: tuple[str, ...]) -> Any:
    """Return the first non-empty nested value among the supplied paths."""

    for path in paths:
        current: Any = mapping
        for key in path:
            if not isinstance(current, dict) or key not in current:
                current = None
                break
            current = current[key]
        if current not in (None, ""):
            return current
    return None


def hardware_id(record: dict[str, Any]) -> str | None:
    """Extract the stable hardware scope identifier from supported locations."""

    value = get_nested(
        record,
        ("provenance", "hardware_id"),
        ("hardware_id",),
        ("hardware", "id"),
        ("provenance", "hardware", "id"),
    )
    return str(value) if value is not None else None


def provenance_value(record: dict[str, Any], key: str) -> Any:
    """Read a provenance key, accepting the historical top-level location too."""

    return get_nested(record, ("provenance", key), (key,))


def measurement_value(record: dict[str, Any], key: str) -> Any:
    """Read a measurement key, accepting the historical top-level location too."""

    return get_nested(record, ("measurement", key), (key,))


def as_number(value: Any) -> float | None:
    """Parse a finite number without accepting booleans or non-finite values."""

    if isinstance(value, bool):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if result != result or result in (float("inf"), float("-inf")):
        return None
    return result


def as_int(value: Any) -> int | None:
    """Parse an integer without silently rounding fractional values."""

    if isinstance(value, bool):
        return None
    try:
        result = int(value)
    except (TypeError, ValueError):
        return None
    return result if str(result) == str(value).strip() or isinstance(value, int) else None


def is_candidate_record(value: Any) -> bool:
    """Distinguish result rows from adjacent manifests and JSON schemas."""

    if not isinstance(value, dict):
        return False
    # ``validation.json`` also carries a schema_version; only a status or a
    # sufficiently specific benchmark identity makes a JSON object a result.
    return "status" in value or {"backend", "dataset_id", "run_id"}.issubset(value)


def _decode_text_field(value: str) -> Any:
    value = value.strip()
    if not value:
        return None
    if value[0] in "[{":
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            return value
    return value


def _read_tabular(path: Path, delimiter: str) -> Iterator[LocatedRecord]:
    with path.open(newline="", encoding="utf-8") as handle:
        for line_number, row in enumerate(csv.DictReader(handle, delimiter=delimiter), start=2):
            payload = {key: _decode_text_field(value or "") for key, value in row.items() if key}
            yield LocatedRecord(payload, f"{path}:{line_number}")


def _records_from_json(value: Any, source: str) -> Iterator[LocatedRecord]:
    if isinstance(value, list):
        for index, item in enumerate(value):
            if isinstance(item, dict):
                yield LocatedRecord(item, f"{source}[{index}]")
        return
    if isinstance(value, dict) and isinstance(value.get("records"), list):
        for index, item in enumerate(value["records"]):
            if isinstance(item, dict):
                yield LocatedRecord(item, f"{source}.records[{index}]")
        return
    if isinstance(value, dict):
        yield LocatedRecord(value, source)


def _discover_candidate_files(paths: Iterable[str | Path]) -> tuple[list[Path], list[ReadIssue], list[str]]:
    """Expand supplied files/directories without attempting to decode them."""

    files: list[Path] = []
    issues: list[ReadIssue] = []
    ignored: list[str] = []
    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists():
            issues.append(ReadIssue(str(path), "path does not exist"))
        elif path.is_dir():
            files.extend(
                sorted(
                    candidate.resolve()
                    for candidate in path.rglob("*")
                    if candidate.is_file() and candidate.suffix.lower() in RECORD_SUFFIXES
                )
            )
        elif path.suffix.lower() in RECORD_SUFFIXES:
            files.append(path.resolve())
        else:
            ignored.append(str(path))
    return sorted(set(files)), issues, ignored


def _decode_file(path: Path) -> list[LocatedRecord]:
    """Decode one supported result file into raw object rows."""

    suffix = path.suffix.lower()
    if suffix in {".csv", ".tsv"}:
        return list(_read_tabular(path, "\t" if suffix == ".tsv" else ","))
    if suffix in {".jsonl", ".ndjson"}:
        decoded: list[LocatedRecord] = []
        with path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                value = json.loads(line)
                if not isinstance(value, dict):
                    raise ValueError(f"line {line_number} is not a JSON object")
                decoded.append(LocatedRecord(value, f"{path}:{line_number}"))
        return decoded
    with path.open(encoding="utf-8") as handle:
        return list(_records_from_json(json.load(handle), str(path)))


def read_records(paths: Iterable[str | Path]) -> tuple[list[LocatedRecord], list[ReadIssue], list[str]]:
    """Read JSON, JSONL, CSV, and TSV result rows below files or directories.

    Non-result JSON files (for example ``result-schema.json``) are deliberately
    ignored and returned in ``ignored``.  Malformed files are reported as
    issues, never silently skipped.
    """

    records: list[LocatedRecord] = []
    files, issues, ignored = _discover_candidate_files(paths)
    for path in files:
        try:
            decoded = _decode_file(path)
        except (OSError, ValueError, json.JSONDecodeError, csv.Error) as exc:
            issues.append(ReadIssue(str(path), str(exc)))
            continue

        candidates = [record for record in decoded if is_candidate_record(record.payload)]
        if not candidates:
            ignored.append(str(path))
            continue
        records.extend(candidates)
    return records, issues, ignored


def comparison_signature(record: dict[str, Any]) -> tuple[Any, ...]:
    """Fields that must be identical before points can share a benchmark plot."""

    return (
        record.get("dataset_id"),
        record.get("dataset_window_count"),
        record.get("dimension"),
        record.get("metric"),
        record.get("k"),
        record.get("warmup_queries"),
        record.get("timed_queries"),
        record.get("query_ids_sha256"),
        record.get("ground_truth_ids_sha256"),
        provenance_value(record, "embedding_cache_id"),
        provenance_value(record, "query_selection_id"),
        provenance_value(record, "ground_truth_cache_id"),
        hardware_id(record),
        provenance_value(record, "git_commit"),
        provenance_value(record, "runner_version"),
        measurement_value(record, "query_batch_size"),
        measurement_value(record, "timed_batch_count"),
    )


def configuration_signature(record: dict[str, Any]) -> tuple[Any, ...]:
    """Fields identifying one backend/parameter setting inside a comparison scope."""

    return comparison_signature(record) + (
        record.get("backend"),
        record.get("parameter_label"),
        canonical_json(record.get("parameters", {})),
        record.get("run_id"),
    )
