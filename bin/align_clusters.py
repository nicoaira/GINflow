#!/usr/bin/env python3
"""Align seed clusters with GINFINITY-SW on a padded crop of each HSP."""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from ginfinity_sw import Alignment, ScoringParameters, align, format_alignment

LN2 = math.log(2.0)


def bit_score(raw_score: float, lam: float, k_value: float) -> float:
    return (lam * raw_score - math.log(k_value)) / LN2


def log_evalue(raw_score: float, query_length: int, search_residues: int, lam: float, k_value: float) -> float:
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


CLUSTER_REQUIRED = {
    "cluster_id",
    "query_id",
    "target_id",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
}

ALIGNMENT_COLUMNS = [
    "cluster_id",
    "query_id",
    "target_id",
    "score",
    "bit_score",
    "evalue",
    "evalue_pair",
    "query_start",
    "query_end",
    "target_start",
    "target_end",
    "query_length",
    "target_length",
    "match_count",
    "aligned_columns",
    "seed_count",
    "max_seed_score",
]


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path} is not a JSON object")
    return payload


def load_parameters(path: Path) -> ScoringParameters:
    payload = load_json(path)
    return ScoringParameters(**payload.get("scoring_parameters", payload))


def load_clusters(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path} has no header")
        missing = CLUSTER_REQUIRED - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def load_npz_shards(paths: list[Path]) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    for path in paths:
        with np.load(path) as archive:
            for key in archive.files:
                if key in arrays:
                    raise ValueError(f"duplicate embedding id {key!r} in {path}")
                arrays[key] = np.asarray(archive[key])
    return arrays


def load_graph_records(paths: list[Path]) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    for path in paths:
        payload = load_json(path)
        identifiers = payload.get("identifiers")
        sequences = payload.get("sequences")
        structures = payload.get("structures")
        if not (identifiers and sequences and structures):
            raise ValueError(f"{path} is not a GINFINITY graph sidecar")
        if not (len(identifiers) == len(sequences) == len(structures)):
            raise ValueError(f"{path} identifier/sequence/structure length mismatch")
        for identifier, sequence, structure in zip(identifiers, sequences, structures):
            if identifier in records:
                raise ValueError(f"duplicate record id {identifier!r} in {path}")
            records[identifier] = (sequence, structure)
    return records


def load_records_tsv(path: Path) -> dict[str, tuple[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"transcript_id", "sequence", "secondary_structure"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must have columns {sorted(required)}")
        records = {}
        for row in reader:
            records[row["transcript_id"]] = (row["sequence"], row["secondary_structure"])
    return records


def shift_alignment(result: Alignment, query_offset: int, target_offset: int) -> Alignment:
    if not result.columns:
        return result
    columns = tuple(
        (
            query if query < 0 else query + query_offset,
            target if target < 0 else target + target_offset,
        )
        for query, target in result.columns
    )
    return Alignment(
        result.score,
        (result.query_span[0] + query_offset, result.query_span[1] + query_offset),
        (result.target_span[0] + target_offset, result.target_span[1] + target_offset),
        columns,
        result.rows_processed,
    )


def crop_bounds(start: int, end: int, length: int, pad: int) -> tuple[int, int]:
    return max(0, start - pad), min(length, end + pad)


def align_cluster(
    cluster: dict,
    query_emb: np.ndarray,
    target_emb: np.ndarray,
    params: ScoringParameters,
    pad: int,
    max_cells: int,
) -> Alignment:
    q0, q1 = crop_bounds(int(cluster["query_start"]), int(cluster["query_end"]), query_emb.shape[0], pad)
    t0, t1 = crop_bounds(int(cluster["target_start"]), int(cluster["target_end"]), target_emb.shape[0], pad)
    cells = (q1 - q0) * (t1 - t0)
    if cells > max_cells:
        raise ValueError(
            f"cluster {cluster['cluster_id']} crop {q1 - q0}x{t1 - t0} "
            f"exceeds --max-cells {max_cells}"
        )
    result = align(query_emb[q0:q1], target_emb[t0:t1], params=params, max_cells=max_cells)
    return shift_alignment(result, q0, t0)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alignment-text", type=Path)
    parser.add_argument("--stats-json", type=Path)
    parser.add_argument("--query-embeddings", type=Path, nargs="+", required=True)
    parser.add_argument("--query-metadata", type=Path, nargs="+", required=True)
    parser.add_argument("--target-embeddings", type=Path, nargs="*")
    parser.add_argument("--target-metadata", type=Path, nargs="*")
    parser.add_argument("--database", type=Path, help="Packed faiss/ directory with embeddings.npz and records.tsv")
    parser.add_argument("--evd", type=Path, required=True, help="evd.json from estimate_evd.py")
    parser.add_argument("--pad", type=int, default=32)
    parser.add_argument("--max-cells", type=int, default=16_777_216)
    return parser.parse_args(argv)


def resolve_targets(args: argparse.Namespace) -> tuple[dict[str, np.ndarray], dict[str, tuple[str, str]]]:
    if args.database:
        embeddings = args.database / "embeddings.npz"
        records = args.database / "records.tsv"
        if not embeddings.is_file() or not records.is_file():
            raise ValueError(
                f"{args.database} is missing embeddings.npz/records.tsv; "
                "rebuild the database with this pipeline version"
            )
        return load_npz_shards([embeddings]), load_records_tsv(records)
    if not args.target_embeddings or not args.target_metadata:
        raise ValueError("provide --database or both --target-embeddings and --target-metadata")
    return load_npz_shards(args.target_embeddings), load_graph_records(args.target_metadata)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.pad < 0:
        print("error: --pad must be >= 0", file=sys.stderr)
        return 2
    try:
        clusters = load_clusters(args.clusters)
        params = load_parameters(args.parameters)
        query_emb = load_npz_shards(args.query_embeddings)
        query_meta = load_graph_records(args.query_metadata)
        target_emb, target_meta = resolve_targets(args)
        evd = load_json(args.evd)
        lam = float(evd["lambda"])
        k_value = float(evd["K"])
        db_residues = int(evd["database_residues"])
        if lam <= 0 or k_value <= 0 or db_residues < 1:
            raise ValueError(f"{args.evd} has non-positive λ, K, or database_residues")
    except (OSError, ValueError, TypeError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    rows = []
    texts = []
    skipped = 0
    for cluster in clusters:
        query_id = cluster["query_id"]
        target_id = cluster["target_id"]
        if query_id not in query_emb or query_id not in query_meta:
            print(f"error: query {query_id!r} missing from query embeddings/metadata", file=sys.stderr)
            return 1
        if target_id not in target_emb or target_id not in target_meta:
            print(f"error: target {target_id!r} missing from target embeddings/metadata", file=sys.stderr)
            return 1
        try:
            result = align_cluster(
                cluster,
                query_emb[query_id],
                target_emb[target_id],
                params,
                args.pad,
                args.max_cells,
            )
        except ValueError as exc:
            print(f"warning: {exc}", file=sys.stderr)
            skipped += 1
            continue
        if result.score <= 0 or not result.columns:
            skipped += 1
            continue
        query_length = int(query_emb[query_id].shape[0])
        target_length = int(target_emb[target_id].shape[0])
        bits = bit_score(result.score, lam, k_value)
        log_db = log_evalue(result.score, query_length, db_residues, lam, k_value)
        log_pair = log_evalue(result.score, query_length, target_length, lam, k_value)
        evalue_db = evalue_from_log(log_db)
        rows.append({
            "cluster_id": cluster["cluster_id"],
            "query_id": query_id,
            "target_id": target_id,
            "score": result.score,
            "bit_score": bits,
            "evalue": evalue_db,
            "evalue_pair": evalue_from_log(log_pair),
            "log_evalue": log_db,
            "query_start": result.query_span[0],
            "query_end": result.query_span[1],
            "target_start": result.target_span[0],
            "target_end": result.target_span[1],
            "query_length": query_length,
            "target_length": target_length,
            "match_count": result.match_count,
            "aligned_columns": len(result.columns),
            "seed_count": cluster.get("seed_count", ""),
            "max_seed_score": cluster.get("max_score", ""),
            "result": result,
        })

    rows.sort(key=lambda item: (
        item["log_evalue"],
        -item["score"],
        item["query_id"],
        item["target_id"],
        str(item["cluster_id"]),
    ))
    formatted_rows = []
    texts = []
    for item in rows:
        formatted_rows.append({
            "cluster_id": item["cluster_id"],
            "query_id": item["query_id"],
            "target_id": item["target_id"],
            "score": f"{item['score']:.6f}",
            "bit_score": f"{item['bit_score']:.3f}",
            "evalue": f"{item['evalue']:.6g}",
            "evalue_pair": f"{item['evalue_pair']:.6g}",
            "query_start": item["query_start"],
            "query_end": item["query_end"],
            "target_start": item["target_start"],
            "target_end": item["target_end"],
            "query_length": item["query_length"],
            "target_length": item["target_length"],
            "match_count": item["match_count"],
            "aligned_columns": item["aligned_columns"],
            "seed_count": item["seed_count"],
            "max_seed_score": item["max_seed_score"],
        })
        if args.alignment_text:
            q_seq, q_struct = query_meta[item["query_id"]]
            t_seq, t_struct = target_meta[item["target_id"]]
            header = (
                f"# cluster {item['cluster_id']}  {item['query_id']} vs {item['target_id']}  "
                f"score={item['score']:.4f}  bits={item['bit_score']:.2f}  E={item['evalue']:.3g}"
            )
            texts.append(header + "\n" + format_alignment(item["result"], q_seq, q_struct, t_seq, t_struct))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ALIGNMENT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(formatted_rows)
    if args.alignment_text:
        args.alignment_text.write_text("\n\n".join(texts) + ("\n" if texts else ""))
    stats = {
        "clusters": len(clusters),
        "alignments": len(formatted_rows),
        "skipped": skipped,
        "lambda": lam,
        "K": k_value,
        "database_residues": db_residues,
    }
    if args.stats_json:
        args.stats_json.write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
