#!/usr/bin/env python3
"""Collapse HSP rows into one BLAST-style result per query-target pair."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


CLUSTER_HEAD = re.compile(
    r"^# cluster\s+(\S+)\s+(\S+)\s+vs\s+(\S+)",
)
LN2 = math.log(2.0)

OUTPUT_COLUMNS = [
    "cluster_id",
    "query_id",
    "target_id",
    "score",
    "total_score",
    "max_score",
    "bit_score",
    "evalue",
    "evalue_pair",
    "log_evalue",
    "alignment_count",
    "cluster_ids",
    "hsp_scores",
    "hsp_spans",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "query_length",
    "target_length",
    "match_count",
    "aligned_columns",
    "ungapped_columns",
    "base_matches",
    "structure_matches",
    "seed_count",
    "max_seed_score",
    "query_sequence",
    "query_structure",
    "target_sequence",
    "target_structure",
    "query_aligned",
    "target_aligned",
]


def parse_float(value: str | None, default: float = math.inf) -> float:
    if value is None or value == "":
        return default
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_int(value: str | None, default: int = 0) -> int:
    if value is None or value == "":
        return default
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def format_float(value: float) -> str:
    if math.isinf(value):
        return "inf"
    if value == 0.0:
        return "0"
    return f"{value:.6g}"


def bit_score(raw_score: float, lam: float, k_value: float) -> float:
    return (lam * raw_score - math.log(k_value)) / LN2


def log_evalue(
    raw_score: float,
    query_length: int,
    search_residues: int,
    lam: float,
    k_value: float,
) -> float:
    return (
        math.log(k_value)
        + math.log(max(query_length, 1))
        + math.log(max(search_residues, 1))
        - lam * raw_score
    )


def evalue_from_log(log_e: float) -> float:
    if log_e > 700:
        return math.inf
    if log_e < -745:
        return 0.0
    return math.exp(log_e)


def load_evd(path: Path | None) -> dict[str, float] | None:
    if path is None:
        return None
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    try:
        lam = float(payload["lambda"])
        k_raw = payload.get("K")
        if k_raw is None:
            k_raw = payload["k"]
        k_value = float(k_raw)
        database_residues = int(payload["database_residues"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{path} is missing valid lambda, K, or database_residues") from exc
    if lam <= 0 or k_value <= 0 or database_residues < 1:
        raise ValueError(f"{path} has non-positive lambda, K, or database_residues")
    return {
        "lambda": lam,
        "K": k_value,
        "database_residues": database_residues,
    }


def load_texts(paths: list[Path]) -> dict[tuple[str, str, str], str]:
    blocks: dict[tuple[str, str, str], list[str]] = {}
    key: tuple[str, str, str] | None = None
    current: list[str] = []
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            continue
        for line in path.read_text().splitlines():
            match = CLUSTER_HEAD.match(line)
            if match:
                if key is not None:
                    blocks[key] = current
                key = (match.group(1), match.group(2), match.group(3))
                current = [line]
                continue
            if key is not None:
                current.append(line)
        if key is not None:
            blocks[key] = current
            key = None
            current = []
    return {item: "\n".join(lines).rstrip() for item, lines in blocks.items()}


def load_rows(paths: list[Path]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for path in paths:
        if not path.is_file():
            continue
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                continue
            rows.extend(row for row in reader if any(row.values()))
    return rows


def row_score(row: dict[str, str]) -> float:
    value = parse_float(row.get("score"), math.nan)
    if not math.isfinite(value):
        value = parse_float(row.get("total_score"), math.nan)
    if not math.isfinite(value):
        raise ValueError(
            f"alignment row for {row.get('query_id', '?')} vs "
            f"{row.get('target_id', '?')} has no finite score"
        )
    return value


def component(row: dict[str, str]) -> dict[str, object]:
    """Normalize either a raw HSP row or a previously collapsed row."""
    count = max(parse_int(row.get("alignment_count"), 1), 1)
    total = parse_float(row.get("total_score"), math.nan)
    if not math.isfinite(total):
        total = row_score(row)
    maximum = parse_float(row.get("max_score"), math.nan)
    if not math.isfinite(maximum):
        maximum = total
    return {"row": row, "count": count, "total": total, "max": maximum}


def optional_sum(rows: list[dict[str, str]], field: str) -> str:
    values = [parse_float(row.get(field), math.nan) for row in rows]
    if not values or not all(math.isfinite(value) for value in values):
        return ""
    return format_float(sum(values))


def aggregate_pair(
    rows: list[dict[str, str]],
    evd: dict[str, float] | None,
) -> dict[str, str]:
    if not rows:
        raise ValueError("cannot aggregate an empty query-target group")
    ordered = sorted(
        rows,
        key=lambda row: (
            parse_int(row.get("query_start")),
            parse_int(row.get("target_start")),
            parse_int(row.get("query_end")),
            parse_int(row.get("target_end")),
            row.get("cluster_id", ""),
        ),
    )
    parts = [component(row) for row in ordered]
    total_score = sum(float(part["total"]) for part in parts)
    max_score = max(float(part["max"]) for part in parts)
    alignment_count = sum(int(part["count"]) for part in parts)

    query_id = ordered[0].get("query_id", "")
    target_id = ordered[0].get("target_id", "")
    query_lengths = {parse_int(row.get("query_length")) for row in ordered}
    target_lengths = {parse_int(row.get("target_length")) for row in ordered}
    if len(query_lengths) != 1 or len(target_lengths) != 1:
        raise ValueError(f"inconsistent sequence lengths for {query_id} vs {target_id}")
    query_length = query_lengths.pop()
    target_length = target_lengths.pop()

    if evd is None:
        if len(parts) != 1 or alignment_count != 1:
            raise ValueError(
                "--evd is required when a query-target pair has multiple alignments"
            )
        first = ordered[0]
        database_evalue = parse_float(first.get("evalue"))
        pair_evalue = parse_float(first.get("evalue_pair"))
        pair_bits = parse_float(first.get("bit_score"), math.nan)
        log_db = parse_float(first.get("log_evalue"), math.nan)
    else:
        lam = evd["lambda"]
        k_value = evd["K"]
        database_residues = int(evd["database_residues"])
        log_db = log_evalue(
            total_score, query_length, database_residues, lam, k_value
        )
        log_pair = log_evalue(
            total_score, query_length, target_length, lam, k_value
        )
        database_evalue = evalue_from_log(log_db)
        pair_evalue = evalue_from_log(log_pair)
        pair_bits = bit_score(total_score, lam, k_value)

    cluster_ids: list[str] = []
    for row in ordered:
        cluster_id = row.get("cluster_id", "")
        if cluster_id not in cluster_ids:
            cluster_ids.append(cluster_id)
    spans = [
        {
            "cluster_id": row.get("cluster_id", ""),
            "query_start": parse_int(row.get("query_start")),
            "query_end": parse_int(row.get("query_end")),
            "target_start": parse_int(row.get("target_start")),
            "target_end": parse_int(row.get("target_end")),
            "score": float(part["total"]),
        }
        for row, part in zip(ordered, parts)
    ]
    result = dict(ordered[0])
    result.update({
        "cluster_id": cluster_ids[0],
        "query_id": query_id,
        "target_id": target_id,
        # score is retained as the public legacy alias for the pair total.
        "score": format_float(total_score),
        "total_score": format_float(total_score),
        "max_score": format_float(max_score),
        "bit_score": format_float(pair_bits),
        "evalue": format_float(database_evalue),
        "evalue_pair": format_float(pair_evalue),
        "log_evalue": format_float(log_db),
        "alignment_count": str(alignment_count),
        "cluster_ids": ",".join(cluster_ids),
        "hsp_scores": json.dumps(
            [round(float(part["total"]), 6) for part in parts], separators=(",", ":")
        ),
        "hsp_spans": json.dumps(spans, separators=(",", ":")),
        "query_start": str(min(parse_int(row.get("query_start")) for row in ordered)),
        "query_end": str(max(parse_int(row.get("query_end")) for row in ordered)),
        "target_start": str(min(parse_int(row.get("target_start")) for row in ordered)),
        "target_end": str(max(parse_int(row.get("target_end")) for row in ordered)),
        "query_length": str(query_length),
        "target_length": str(target_length),
        "match_count": str(sum(parse_int(row.get("match_count")) for row in ordered)),
        "aligned_columns": str(sum(parse_int(row.get("aligned_columns")) for row in ordered)),
        "ungapped_columns": str(sum(parse_int(row.get("ungapped_columns")) for row in ordered)),
        "base_matches": optional_sum(ordered, "base_matches"),
        "structure_matches": optional_sum(ordered, "structure_matches"),
        "seed_count": str(sum(parse_int(row.get("seed_count")) for row in ordered)),
    })
    max_seed_scores = [
        parse_float(row.get("max_seed_score"), math.nan) for row in ordered
    ]
    result["max_seed_score"] = (
        format_float(max(value for value in max_seed_scores if math.isfinite(value)))
        if any(math.isfinite(value) for value in max_seed_scores)
        else ""
    )
    for field in OUTPUT_COLUMNS:
        result.setdefault(field, "")
    return {field: str(result.get(field, "")) for field in OUTPUT_COLUMNS}


def aggregate_rows(
    rows: list[dict[str, str]],
    evd: dict[str, float] | None,
) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        groups.setdefault(
            (row.get("query_id", ""), row.get("target_id", "")), []
        ).append(row)
    merged = [aggregate_pair(group, evd) for group in groups.values()]
    merged.sort(
        key=lambda row: (
            parse_float(row.get("evalue")),
            -parse_float(row.get("total_score"), 0.0),
            -parse_float(row.get("max_score"), 0.0),
            row.get("query_id", ""),
            row.get("target_id", ""),
        )
    )
    return merged


def pair_text(
    row: dict[str, str],
    source_rows: list[dict[str, str]],
    texts: dict[tuple[str, str, str], str],
) -> str:
    query_id = row.get("query_id", "")
    target_id = row.get("target_id", "")
    header = (
        f"# cluster {row.get('cluster_id', '')}  {query_id} vs {target_id}  "
        f"total_score={row.get('total_score', '')}  "
        f"max_score={row.get('max_score', '')}  "
        f"bits={row.get('bit_score', '')}  E={row.get('evalue', '')}  "
        f"E_pair={row.get('evalue_pair', '')}  "
        f"HSPs={row.get('alignment_count', '1')}"
    )
    blocks: list[str] = []
    ordered = sorted(
        source_rows,
        key=lambda item: (
            parse_int(item.get("query_start")),
            parse_int(item.get("target_start")),
            item.get("cluster_id", ""),
        ),
    )
    for index, source in enumerate(ordered, start=1):
        key = (source.get("cluster_id", ""), query_id, target_id)
        raw = texts.get(key, "")
        body_lines = raw.splitlines()
        if body_lines and CLUSTER_HEAD.match(body_lines[0]):
            body_lines = body_lines[1:]
        body = "\n".join(body_lines).strip()
        label = (
            f"# HSP {index}/{len(ordered)} cluster={source.get('cluster_id', '')} "
            f"score={format_float(row_score(source))}"
        )
        blocks.append(label + ("\n" + body if body else ""))
    return header + ("\n\n" + "\n\n".join(blocks) if blocks else "")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignments", type=Path, nargs="*", default=[])
    parser.add_argument("--alignment-text", type=Path, nargs="*", default=[])
    parser.add_argument("--evd", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alignment-text-output", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        rows = load_rows(args.alignments)
        evd = load_evd(args.evd)
        merged = aggregate_rows(rows, evd)
    except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
        print(f"error: {exc}")
        return 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=OUTPUT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(merged)

    if args.alignment_text_output:
        texts = load_texts(args.alignment_text)
        groups: dict[tuple[str, str], list[dict[str, str]]] = {}
        for source in rows:
            groups.setdefault(
                (source.get("query_id", ""), source.get("target_id", "")), []
            ).append(source)
        blocks = []
        for row in merged:
            source_rows = groups[(row.get("query_id", ""), row.get("target_id", ""))]
            blocks.append(pair_text(row, source_rows, texts))
        args.alignment_text_output.write_text(
            "\n\n".join(blocks) + ("\n" if blocks else "")
        )

    print(
        f"merged {len(rows)} HSPs into {len(merged)} query-target pairs "
        f"from {len(args.alignments)} tables"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
