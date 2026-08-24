#!/usr/bin/env python3
"""Align seed clusters with GINFINITY-SW on a padded crop of each HSP.

This process deliberately emits one row per seed cluster. MERGE_ALIGNMENTS is
the pair-level boundary: it combines all rows with the same query and target
into one BLAST-style result.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np
from ginfinity_sw import Alignment, ScoringParameters, align, format_alignment, transform_scores

sys.path.insert(0, str(Path(__file__).resolve().parent))
from record_pack import load_embedding_files, load_residue_embeddings  # noqa: E402
from sw_batch import align_score_matrices, pin_blas_threads  # noqa: E402

pin_blas_threads()

LN2 = math.log(2.0)
SLICE_ID_RE = re.compile(r"^(?P<base>.+):(?P<start>\d+)-(?P<end>\d+)$")


def parse_slice_id(identifier: str) -> tuple[int, int] | None:
    match = SLICE_ID_RE.fullmatch(identifier)
    if not match:
        return None
    return int(match["start"]), int(match["end"])


def pair_table(structure: str) -> list[int]:
    partners = [-1] * len(structure)
    stack: list[int] = []
    for index, char in enumerate(structure):
        if char == "(":
            stack.append(index)
        elif char == ")":
            if not stack:
                continue
            opening = stack.pop()
            partners[opening] = index
            partners[index] = opening
    return partners


def close_window_structure(structure: str, start: int, end: int) -> str:
    """Keep only pairs whose both ends lie inside ``[start, end)``."""
    partners = pair_table(structure)
    out: list[str] = []
    for index in range(start, end):
        partner = partners[index]
        if partner < start or partner >= end:
            out.append(".")
        else:
            out.append("(" if partner > index else ")")
    return "".join(out)


def subject_sequence(identifier: str, sequence: str, structure: str) -> tuple[str, str]:
    """Return the independent subject (core window, or the full molecule)."""
    if len(sequence) != len(structure):
        raise ValueError(f"{identifier} sequence/structure length mismatch")
    span = parse_slice_id(identifier)
    if span is None:
        return sequence, structure
    start, end = span
    if end - start == len(sequence):
        return sequence, close_window_structure(structure, 0, len(structure))
    if not (0 <= start < end <= len(sequence)):
        raise ValueError(
            f"{identifier} slice [{start}, {end}) is outside a {len(sequence)} nt sequence"
        )
    return sequence[start:end], close_window_structure(structure, start, end)


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
    "total_score",
    "max_score",
    "bit_score",
    "evalue",
    "evalue_pair",
    "log_evalue",
    "alignment_count",
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


def aligned_strings(result: Alignment, query_seq: str, target_seq: str) -> tuple[str, str]:
    query_parts = []
    target_parts = []
    for query, target in result.columns:
        query_parts.append("-" if query < 0 else query_seq[query])
        target_parts.append("-" if target < 0 else target_seq[target])
    return "".join(query_parts), "".join(target_parts)


def structure_state(character: str) -> int:
    return 0 if character == "(" else (2 if character == ")" else 1)


def identity_counts(
    result: Alignment,
    query_seq: str,
    query_structure: str,
    target_seq: str,
    target_structure: str,
) -> tuple[int, int, int]:
    """Return base matches, structure-state matches, and ungapped columns."""
    base_matches = 0
    structure_matches = 0
    ungapped = 0
    for query, target in result.columns:
        if query < 0 or target < 0:
            continue
        ungapped += 1
        base_matches += int(
            query_seq[query].upper().replace("T", "U")
            == target_seq[target].upper().replace("T", "U")
        )
        structure_matches += int(
            structure_state(query_structure[query])
            == structure_state(target_structure[target])
        )
    return base_matches, structure_matches, ungapped


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


def normalize_residue_dict(arrays: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    normalized: dict[str, np.ndarray] = {}
    for identifier, matrix in arrays.items():
        values = np.ascontiguousarray(matrix, dtype=np.float64)
        norms = np.linalg.norm(values, axis=1, keepdims=True)
        normalized[identifier] = values / np.maximum(norms, 1e-12)
    return normalized


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
            records[identifier] = subject_sequence(identifier, sequence, structure)
    return records


def load_records_tsv(path: Path) -> dict[str, tuple[str, str]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"transcript_id", "sequence", "secondary_structure"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must have columns {sorted(required)}")
        records = {}
        for row in reader:
            identifier = row["transcript_id"]
            records[identifier] = subject_sequence(
                identifier, row["sequence"], row["secondary_structure"]
            )
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
    parser.add_argument(
        "--cpus",
        type=int,
        default=1,
        help="independent-pair SW workers (Numba threads)",
    )
    return parser.parse_args(argv)


def resolve_targets(
    args: argparse.Namespace,
    ids: set[str] | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, tuple[str, str]]]:
    if args.database:
        records = args.database / "records.tsv"
        packed_npy = args.database / "embeddings.vectors.npy"
        packed_npz = args.database / "embeddings.npz"
        if not records.is_file() or not (packed_npy.is_file() or packed_npz.is_file()):
            raise ValueError(
                f"{args.database} is missing records.tsv or residue embeddings; "
                "rebuild the database with this pipeline version"
            )
        return load_residue_embeddings(args.database, ids=ids), load_records_tsv(records)
    if not args.target_embeddings or not args.target_metadata:
        raise ValueError("provide --database or both --target-embeddings and --target-metadata")
    return load_embedding_files(args.target_embeddings, ids=ids), load_graph_records(args.target_metadata)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.pad < 0:
        print("error: --pad must be >= 0", file=sys.stderr)
        return 2
    if args.cpus < 1:
        print("error: --cpus must be >= 1", file=sys.stderr)
        return 2
    load_started = time.perf_counter()
    try:
        clusters = load_clusters(args.clusters)
        params = load_parameters(args.parameters)
        query_ids = {cluster["query_id"] for cluster in clusters}
        target_ids = {cluster["target_id"] for cluster in clusters}
        query_emb = normalize_residue_dict(
            load_embedding_files(args.query_embeddings, ids=query_ids or None)
        )
        query_meta = load_graph_records(args.query_metadata)
        target_emb, target_meta = resolve_targets(args, ids=target_ids or None)
        target_emb = normalize_residue_dict(target_emb)
        evd = load_json(args.evd)
        lam = float(evd["lambda"])
        k_value = float(evd["K"])
        db_residues = int(evd["database_residues"])
        if lam <= 0 or k_value <= 0 or db_residues < 1:
            raise ValueError(f"{args.evd} has non-positive λ, K, or database_residues")
    except (OSError, ValueError, TypeError, KeyError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    load_seconds = time.perf_counter() - load_started

    jobs: list[tuple[dict, int, int, np.ndarray]] = []
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
        query_matrix = query_emb[query_id]
        target_matrix = target_emb[target_id]
        q0, q1 = crop_bounds(int(cluster["query_start"]), int(cluster["query_end"]), query_matrix.shape[0], args.pad)
        t0, t1 = crop_bounds(int(cluster["target_start"]), int(cluster["target_end"]), target_matrix.shape[0], args.pad)
        cells = (q1 - q0) * (t1 - t0)
        if cells > args.max_cells:
            print(
                f"warning: cluster {cluster['cluster_id']} crop {q1 - q0}x{t1 - t0} "
                f"exceeds --max-cells {args.max_cells}",
                file=sys.stderr,
            )
            skipped += 1
            continue
        if cells == 0:
            skipped += 1
            continue
        scores = transform_scores(query_matrix[q0:q1] @ target_matrix[t0:t1].T, params)
        jobs.append((cluster, q0, t0, scores))

    align_started = time.perf_counter()
    alignments = align_score_matrices(
        [scores for _, _, _, scores in jobs],
        params,
        traceback=True,
        max_cells=args.max_cells,
        workers=args.cpus,
    ) if jobs else []
    align_seconds = time.perf_counter() - align_started

    rows = []
    for (cluster, q0, t0, _scores), raw in zip(jobs, alignments):
        result = shift_alignment(raw, q0, t0)
        if result.score <= 0 or not result.columns:
            skipped += 1
            continue
        query_id = cluster["query_id"]
        target_id = cluster["target_id"]
        query_length = int(query_emb[query_id].shape[0])
        target_length = int(target_emb[target_id].shape[0])
        bits = bit_score(result.score, lam, k_value)
        log_db = log_evalue(result.score, query_length, db_residues, lam, k_value)
        log_pair = log_evalue(result.score, query_length, target_length, lam, k_value)
        evalue_db = evalue_from_log(log_db)
        query_aligned, target_aligned = aligned_strings(
            result, query_meta[query_id][0], target_meta[target_id][0]
        )
        base_matches, structure_matches, ungapped_columns = identity_counts(
            result,
            query_meta[query_id][0],
            query_meta[query_id][1],
            target_meta[target_id][0],
            target_meta[target_id][1],
        )
        rows.append({
            "cluster_id": cluster["cluster_id"],
            "query_id": query_id,
            "target_id": target_id,
            "score": result.score,
            "total_score": result.score,
            "max_score": result.score,
            "bit_score": bits,
            "evalue": evalue_db,
            "evalue_pair": evalue_from_log(log_pair),
            "log_evalue": log_db,
            "alignment_count": 1,
            "query_start": result.query_span[0],
            "query_end": result.query_span[1],
            "target_start": result.target_span[0],
            "target_end": result.target_span[1],
            "query_length": query_length,
            "target_length": target_length,
            "match_count": result.match_count,
            "aligned_columns": len(result.columns),
            "ungapped_columns": ungapped_columns,
            "base_matches": base_matches,
            "structure_matches": structure_matches,
            "seed_count": cluster.get("seed_count", ""),
            "max_seed_score": cluster.get("max_score", ""),
            "query_sequence": query_meta[query_id][0],
            "query_structure": query_meta[query_id][1],
            "target_sequence": target_meta[target_id][0],
            "target_structure": target_meta[target_id][1],
            "query_aligned": query_aligned,
            "target_aligned": target_aligned,
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
            "total_score": f"{item['total_score']:.6f}",
            "max_score": f"{item['max_score']:.6f}",
            "bit_score": f"{item['bit_score']:.3f}",
            "evalue": f"{item['evalue']:.6g}",
            "evalue_pair": f"{item['evalue_pair']:.6g}",
            "log_evalue": f"{item['log_evalue']:.6g}",
            "alignment_count": item["alignment_count"],
            "query_start": item["query_start"],
            "query_end": item["query_end"],
            "target_start": item["target_start"],
            "target_end": item["target_end"],
            "query_length": item["query_length"],
            "target_length": item["target_length"],
            "match_count": item["match_count"],
            "aligned_columns": item["aligned_columns"],
            "ungapped_columns": item["ungapped_columns"],
            "base_matches": item["base_matches"],
            "structure_matches": item["structure_matches"],
            "seed_count": item["seed_count"],
            "max_seed_score": item["max_seed_score"],
            "query_sequence": item["query_sequence"],
            "query_structure": item["query_structure"],
            "target_sequence": item["target_sequence"],
            "target_structure": item["target_structure"],
            "query_aligned": item["query_aligned"],
            "target_aligned": item["target_aligned"],
        })
        if args.alignment_text:
            q_seq, q_struct = query_meta[item["query_id"]]
            t_seq, t_struct = target_meta[item["target_id"]]
            header = (
                f"# cluster {item['cluster_id']}  {item['query_id']} vs {item['target_id']}  "
                f"hsp_score={item['score']:.4f}  bits={item['bit_score']:.2f}  "
                f"E_hsp={item['evalue']:.3g}"
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
        "hsp_alignments": len(formatted_rows),
        "skipped": skipped,
        "lambda": lam,
        "K": k_value,
        "database_residues": db_residues,
        "workers": args.cpus,
        "load_seconds": round(load_seconds, 3),
        "align_seconds": round(align_seconds, 3),
    }
    if args.stats_json:
        args.stats_json.write_text(json.dumps(stats, indent=2) + "\n")
    print(json.dumps(stats))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
