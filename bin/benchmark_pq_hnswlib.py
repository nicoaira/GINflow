#!/usr/bin/env python3
"""Benchmark node-level product-quantized windows with custom C++ HNSW.

PQ is trained on the original 128-dimensional node embeddings. Each window
stores the packed PQ code for every node, and the C++ distance sums
subquantizer lookup scores at matching window positions. Final ranking is
optionally recomputed with the original normalized windows.
"""
from __future__ import annotations

import argparse
import gc
import json
import resource
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np

from benchmark_hnswlib import load_queries, load_structures, normalize_rows, source_embeddings
from benchmark_faiss_hnsw_int8 import recall_at


def parse_int_list(value: str) -> list[int]:
    values = [int(item.strip()) for item in str(value).split(",") if item.strip()]
    if not values or any(item < 1 for item in values):
        raise ValueError(f"expected positive comma-separated integers, got {value!r}")
    return values


def unpack_pq_codes(packed: np.ndarray, pq_m: int, nbits: int) -> np.ndarray:
    """Expand FAISS's little-endian packed PQ codes to one byte per subcode."""
    values = np.asarray(packed, dtype=np.uint8)
    if values.ndim != 2:
        raise ValueError(f"packed PQ codes must be 2-D, got {values.shape}")
    bits = np.unpackbits(values, axis=1, bitorder="little")
    needed = pq_m * nbits
    if bits.shape[1] < needed:
        raise ValueError("packed PQ code buffer is too short")
    bit_rows = bits[:, :needed].reshape(values.shape[0], pq_m, nbits)
    weights = (np.uint16(1) << np.arange(nbits, dtype=np.uint16))[None, None, :]
    return np.asarray((bit_rows * weights).sum(axis=2), dtype=np.uint8)


def pack_pq_codes(codes: np.ndarray, pq_m: int, nbits: int) -> np.ndarray:
    """Pack one-byte subcodes into FAISS-compatible little-endian bytes."""
    values = np.asarray(codes, dtype=np.uint8)
    if values.ndim != 2 or values.shape[1] != pq_m:
        raise ValueError(f"PQ codes must have shape (n, {pq_m}), got {values.shape}")
    if values.size and int(values.max()) >= (1 << nbits):
        raise ValueError("PQ code exceeds the requested subquantizer range")
    bits = ((values[:, :, None] >> np.arange(nbits, dtype=np.uint8)) & 1).reshape(
        values.shape[0], pq_m * nbits
    )
    return np.packbits(bits, axis=1, bitorder="little")


def normalize_nodes(values: np.ndarray, original_dtype: str) -> np.ndarray:
    source = np.asarray(values, dtype=np.float16 if original_dtype == "float16" else np.float32)
    return normalize_rows(source)


def rerank_candidates_allow_missing(
    database_windows: np.ndarray,
    query_windows: np.ndarray,
    labels: np.ndarray,
    output_k: int,
    batch_size: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Rerank C++ results while preserving rows with fewer than k neighbors."""
    if labels.ndim != 2 or labels.shape[0] != query_windows.shape[0]:
        raise ValueError("candidate labels and queries have incompatible shapes")
    if output_k < 1 or output_k > labels.shape[1]:
        raise ValueError("output_k must be within the candidate-label width")
    final_labels = np.empty((labels.shape[0], output_k), dtype=np.int64)
    final_scores = np.empty((labels.shape[0], output_k), dtype=np.float32)
    started = time.perf_counter()
    for start in range(0, labels.shape[0], batch_size):
        stop = min(start + batch_size, labels.shape[0])
        batch_labels = np.asarray(labels[start:stop], dtype=np.int64)
        valid = (batch_labels >= 0) & (batch_labels < database_windows.shape[0])
        safe_labels = np.where(valid, batch_labels, 0)
        candidates = np.asarray(database_windows[safe_labels], dtype=np.float32)
        scores = np.einsum(
            "bd,bkd->bk", query_windows[start:stop], candidates, optimize=True
        )
        scores[~valid] = -np.inf
        order = np.argsort(-scores, axis=1, kind="stable")[:, :output_k]
        rows = np.arange(stop - start)[:, None]
        final_labels[start:stop] = batch_labels[rows, order]
        final_scores[start:stop] = scores[rows, order]
    return final_labels, final_scores, time.perf_counter() - started


def pq_codebook_array(pq: Any, pq_m: int, ksub: int, dimension: int) -> np.ndarray:
    import faiss

    dsub = dimension // pq_m
    values = faiss.vector_to_array(pq.centroids)
    return np.asarray(values, dtype=np.float32).reshape(pq_m, ksub, dsub)


def train_and_encode_nodes(
    embeddings: np.ndarray,
    pq_m: int,
    nbits: int,
    sample_size: int,
    niter: int,
    seed: int,
    batch_size: int,
    original_dtype: str,
    output_dir: Path,
) -> dict[str, Any]:
    import faiss

    dimension = int(embeddings.shape[1])
    if dimension % pq_m != 0:
        raise ValueError(f"embedding dimension {dimension} is not divisible by PQ M={pq_m}")
    ksub = 1 << nbits
    node_code_bytes = (pq_m * nbits + 7) // 8
    node_codes_path = output_dir / "node_codes.npy"
    codebook_path = output_dir / "codebook.float32.npy"
    similarity_path = output_dir / "similarity.float32.bin"
    metadata_path = output_dir / "pq.json"
    expected_meta = {
        "pq_m": pq_m,
        "nbits": nbits,
        "ksub": ksub,
        "dimension": dimension,
        "n_nodes": int(embeddings.shape[0]),
        "node_code_bytes": node_code_bytes,
        "sample_size": min(sample_size, int(embeddings.shape[0])),
        "niter": niter,
        "seed": seed,
        "original_dtype": original_dtype,
    }
    if metadata_path.exists() and node_codes_path.exists() and codebook_path.exists() and similarity_path.exists():
        try:
            metadata = json.loads(metadata_path.read_text())
            cached_codes = np.load(node_codes_path, mmap_mode="r")
            if (
                all(metadata.get(key) == value for key, value in expected_meta.items())
                and cached_codes.shape == (embeddings.shape[0], node_code_bytes)
                and cached_codes.dtype == np.uint8
            ):
                return metadata
        except (OSError, ValueError, TypeError):
            pass

    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    sample_count = min(sample_size, int(embeddings.shape[0]))
    sample_indices = rng.choice(embeddings.shape[0], size=sample_count, replace=False)
    training = normalize_nodes(embeddings[sample_indices], original_dtype)
    pq = faiss.ProductQuantizer(dimension, pq_m, nbits)
    pq.cp.niter = int(niter)
    pq.cp.nredo = 1
    pq.cp.verbose = False
    train_started = time.perf_counter()
    pq.train(np.ascontiguousarray(training, dtype=np.float32))
    train_seconds = time.perf_counter() - train_started
    codebook = pq_codebook_array(pq, pq_m, ksub, dimension)
    similarity = np.asarray(np.matmul(codebook, np.swapaxes(codebook, 1, 2)), dtype=np.float32)
    np.save(codebook_path, codebook)
    np.ascontiguousarray(similarity, dtype=np.float32).tofile(similarity_path)

    node_codes = np.lib.format.open_memmap(
        node_codes_path,
        mode="w+",
        dtype=np.uint8,
        shape=(embeddings.shape[0], node_code_bytes),
    )
    encode_started = time.perf_counter()
    for start in range(0, embeddings.shape[0], batch_size):
        stop = min(start + batch_size, embeddings.shape[0])
        values = normalize_nodes(embeddings[start:stop], original_dtype)
        packed = np.asarray(pq.compute_codes(np.ascontiguousarray(values, dtype=np.float32)), dtype=np.uint8)
        if packed.shape != (stop - start, node_code_bytes):
            raise ValueError(f"FAISS returned packed codes with unexpected shape {packed.shape}")
        node_codes[start:stop] = packed
    node_codes.flush()
    encode_seconds = time.perf_counter() - encode_started
    del node_codes, pq, training
    gc.collect()

    metadata = {
        **expected_meta,
        "train_seconds": train_seconds,
        "encode_seconds": encode_seconds,
        "codebook_bytes": int(codebook.nbytes),
        "similarity_bytes": int(similarity.nbytes),
        "node_code_storage_bytes": int(embeddings.shape[0] * node_code_bytes),
        "node_code_path": str(node_codes_path),
        "codebook_path": str(codebook_path),
        "similarity_path": str(similarity_path),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def pack_window_codes(
    node_codes: np.ndarray,
    records: list[Any],
    window_size: int,
    stride: int,
    pq_m: int,
    nbits: int,
    output_path: Path,
) -> int:
    n_windows = int(records[-1].window_start + records[-1].n_windows)
    node_code_bytes = (pq_m * nbits + 7) // 8
    expected_shape = (n_windows, window_size * node_code_bytes)
    cache_metadata_path = output_path.with_name(output_path.name + ".meta.json")
    cache_metadata = {
        "layout": "node-major-packed-v2",
        "n_windows": n_windows,
        "window_size": window_size,
        "stride": stride,
        "pq_m": pq_m,
        "nbits": nbits,
        "node_code_bytes": node_code_bytes,
    }
    if output_path.exists() and output_path.stat().st_size == int(np.prod(expected_shape)):
        try:
            if json.loads(cache_metadata_path.read_text()) == cache_metadata:
                return n_windows
        except (OSError, ValueError):
            pass
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output = np.memmap(output_path, mode="w+", dtype=np.uint8, shape=expected_shape)
    for record in records:
        if record.n_windows == 0:
            continue
        values = np.asarray(node_codes[record.node_start : record.node_start + record.length])
        views = np.lib.stride_tricks.sliding_window_view(values, window_size, axis=0)[::stride]
        rows = views.transpose(0, 2, 1).reshape(record.n_windows, window_size * node_code_bytes)
        output[record.window_start : record.window_start + record.n_windows] = rows
    output.flush()
    del output
    cache_metadata_path.write_text(json.dumps(cache_metadata, indent=2) + "\n")
    return n_windows


def query_window_codes(
    node_codes: np.ndarray,
    queries: list[Any],
    records: list[Any],
    window_size: int,
    stride: int,
    pq_m: int,
    nbits: int,
) -> np.ndarray:
    by_id = {record.identifier: record for record in records}
    node_code_bytes = (pq_m * nbits + 7) // 8
    output = np.empty((len(queries), window_size * node_code_bytes), dtype=np.uint8)
    for row, query in enumerate(queries):
        record = by_id[query.transcript_id]
        start = record.node_start + query.window_offset * stride
        output[row] = np.asarray(node_codes[start : start + window_size]).reshape(-1)
    return output


def run_timed(command: list[str], output_dir: Path, stem: str) -> tuple[float, float]:
    started = time.perf_counter()
    time_binary = shutil.which("/usr/bin/time") or shutil.which("time")
    if time_binary:
        rss_path = output_dir / f"{stem}.time"
        subprocess.run([time_binary, "-v", "-o", str(rss_path), *command], check=True)
    else:
        process = subprocess.Popen(command)
        peak_kib = 0.0
        while process.poll() is None:
            try:
                for line in Path(f"/proc/{process.pid}/status").read_text().splitlines():
                    if line.startswith("VmHWM:"):
                        peak_kib = max(peak_kib, float(line.split()[1]))
                        break
            except (OSError, ValueError):
                pass
            time.sleep(0.01)
        if process.returncode:
            raise subprocess.CalledProcessError(process.returncode, command)
        try:
            for line in Path(f"/proc/{process.pid}/status").read_text().splitlines():
                if line.startswith("VmHWM:"):
                    peak_kib = max(peak_kib, float(line.split()[1]))
                    break
        except (OSError, ValueError):
            pass
        elapsed = time.perf_counter() - started
        if peak_kib:
            return elapsed, peak_kib / 1024.0
        peak_kib = float(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        return elapsed, peak_kib / 1024.0
    elapsed = time.perf_counter() - started
    if time_binary:
        peak_kib = 0.0
        for line in rss_path.read_text().splitlines():
            if "Maximum resident set size" in line:
                peak_kib = float(line.rsplit(":", 1)[1].strip())
                break
        return elapsed, peak_kib * 1024.0 / (1024.0 * 1024.0)
    raise RuntimeError("timing helper did not collect memory")


def run_cpp_search(
    executable: Path,
    index_path: Path,
    similarity_path: Path,
    query_codes_path: Path,
    query_count: int,
    window_size: int,
    pq_m: int,
    nbits: int,
    candidate_k: int,
    ef_search: int,
    threads: int,
    output_dir: Path,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    labels_path = output_dir / f"labels_k{candidate_k}_ef{ef_search}.u64"
    distances_path = output_dir / f"distances_k{candidate_k}_ef{ef_search}.f32"
    metadata_path = output_dir / f"search_k{candidate_k}_ef{ef_search}.json"
    metrics_path = output_dir / f"search_k{candidate_k}_ef{ef_search}.metrics.json"
    expected_labels = query_count * candidate_k * np.dtype("<u8").itemsize
    expected_distances = query_count * candidate_k * np.dtype("<f4").itemsize
    signature = {
        "layout": "node-major-packed-v2",
        "query_count": query_count,
        "window_size": window_size,
        "pq_m": pq_m,
        "nbits": nbits,
        "candidate_k": candidate_k,
        "ef_search": ef_search,
        "index_bytes": index_path.stat().st_size,
    }
    cache_valid = False
    if labels_path.exists() and distances_path.exists() and metadata_path.exists() and metrics_path.exists():
        try:
            cache_valid = json.loads(metadata_path.read_text()) == signature
        except (OSError, ValueError):
            cache_valid = False
    if cache_valid and labels_path.stat().st_size == expected_labels and distances_path.stat().st_size == expected_distances:
        try:
            metrics = json.loads(metrics_path.read_text())
            elapsed = float(metrics["elapsed_seconds"])
            peak_rss_mib = float(metrics["peak_rss_mib"])
        except (KeyError, OSError, TypeError, ValueError):
            cache_valid = False
        else:
            labels = np.fromfile(labels_path, dtype="<u8").reshape(query_count, candidate_k)
            distances = np.fromfile(distances_path, dtype="<f4").reshape(query_count, candidate_k)
            return labels, distances, elapsed, peak_rss_mib
    command = [
        str(executable), "search",
        "--codes", str(query_codes_path),
        "--similarity", str(similarity_path),
        "--query-count", str(query_count),
        "--window-size", str(window_size),
        "--pq-m", str(pq_m),
        "--nbits", str(nbits),
        "--index", str(index_path),
        "--k", str(candidate_k),
        "--ef-search", str(ef_search),
        "--threads", str(threads),
        "--labels-out", str(labels_path),
        "--distances-out", str(distances_path),
    ]
    elapsed, peak_rss_mib = run_timed(command, output_dir, f"search_k{candidate_k}_ef{ef_search}")
    metadata_path.write_text(json.dumps(signature, indent=2) + "\n")
    metrics_path.write_text(json.dumps({
        "elapsed_seconds": elapsed,
        "peak_rss_mib": peak_rss_mib,
    }, indent=2) + "\n")
    labels = np.fromfile(labels_path, dtype="<u8").reshape(query_count, candidate_k)
    distances = np.fromfile(distances_path, dtype="<f4").reshape(query_count, candidate_k)
    return labels, distances, elapsed, peak_rss_mib


def evaluate(args: argparse.Namespace) -> dict[str, Any]:
    pq_ms = parse_int_list(args.pq_m_values)
    nbits_values = parse_int_list(args.nbits_values)
    candidates = parse_int_list(args.candidate_k_values)
    efs = parse_int_list(args.ef_search_values)
    if any(bits > 8 for bits in nbits_values):
        raise ValueError("nbits must be <= 8")
    if any(ef < candidate for ef in efs for candidate in candidates):
        raise ValueError("every ef-search value must be >= every candidate-k value")
    if args.output_k > min(candidates):
        raise ValueError("output-k must not exceed the smallest candidate-k")

    records, node_count, window_count = load_structures(args.structures, args.window_size, args.stride)
    queries, query_stats = load_queries(args.query_selections, records, args.max_queries)
    embeddings = source_embeddings(args.embeddings, node_count, args.original_dtype)
    database_windows = np.load(args.database_windows, mmap_mode="r")
    query_windows = np.asarray(np.load(args.queries), dtype=np.float32)
    query_windows = query_windows[: len(queries)]
    reference = np.asarray(np.load(args.reference_labels), dtype=np.int64)[: len(queries)] if args.reference_labels else None
    if database_windows.shape != (window_count, args.window_size * embeddings.shape[1]):
        raise ValueError(f"database windows shape {database_windows.shape} is incompatible with structures")
    if query_windows.shape[0] != len(queries) or query_windows.shape[1] != database_windows.shape[1]:
        raise ValueError("query windows and query selections are incompatible")
    if reference is not None and reference.shape[0] != len(queries):
        raise ValueError("reference labels and query selections are incompatible")

    args.outdir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, Any]] = []
    executable = args.executable.resolve()
    for pq_m in pq_ms:
        if embeddings.shape[1] % pq_m != 0:
            continue
        for nbits in nbits_values:
            config_dir = args.outdir / f"pq_m{pq_m}_b{nbits}"
            config_dir.mkdir(parents=True, exist_ok=True)
            model = train_and_encode_nodes(
                embeddings, pq_m, nbits, args.sample_size, args.pq_niter, args.seed,
                args.batch_size, args.original_dtype, config_dir,
            )
            node_codes = np.load(model["node_code_path"], mmap_mode="r")
            window_codes_path = config_dir / "database_windows.codes.bin"
            query_codes_path = config_dir / "query_windows.codes.bin"
            pack_started = time.perf_counter()
            actual_windows = pack_window_codes(
                node_codes, records, args.window_size, args.stride, pq_m, nbits, window_codes_path
            )
            query_codes = query_window_codes(
                node_codes, queries, records, args.window_size, args.stride, pq_m, nbits
            )
            query_codes.tofile(query_codes_path)
            pack_seconds = time.perf_counter() - pack_started
            if actual_windows != window_count:
                raise ValueError("packed window count does not match structures")

            for hnsw_m in args.hnsw_m_values:
                for ef_construction in args.ef_construction_values:
                    index_dir = config_dir / f"hnsw_m{hnsw_m}_efc{ef_construction}"
                    index_dir.mkdir(parents=True, exist_ok=True)
                    index_path = index_dir / "index.bin"
                    index_metadata_path = index_dir / "index.json"
                    build_metrics_path = index_dir / "build.metrics.json"
                    index_signature = {
                        "layout": "node-major-packed-v2",
                        "n_windows": window_count,
                        "window_size": args.window_size,
                        "pq_m": pq_m,
                        "nbits": nbits,
                        "hnsw_m": hnsw_m,
                        "ef_construction": ef_construction,
                    }
                    build_peak_rss_mib = 0.0
                    index_valid = False
                    if index_path.exists() and index_metadata_path.exists():
                        try:
                            index_valid = json.loads(index_metadata_path.read_text()) == index_signature
                        except (OSError, ValueError):
                            index_valid = False
                    if not index_valid:
                        build_command = [
                            str(executable), "build",
                            "--codes", str(window_codes_path),
                            "--similarity", str(config_dir / "similarity.float32.bin"),
                            "--count", str(window_count),
                            "--window-size", str(args.window_size),
                            "--pq-m", str(pq_m),
                            "--nbits", str(nbits),
                            "--index", str(index_path),
                            "--m", str(hnsw_m),
                            "--ef-construction", str(ef_construction),
                            "--ef-search", str(max(efs)),
                            "--threads", str(args.threads),
                            "--random-seed", str(args.seed),
                        ]
                        build_seconds, build_peak_rss_mib = run_timed(build_command, index_dir, "build")
                        index_metadata_path.write_text(json.dumps(index_signature, indent=2) + "\n")
                        build_metrics_path.write_text(json.dumps({
                            "elapsed_seconds": build_seconds,
                            "peak_rss_mib": build_peak_rss_mib,
                        }, indent=2) + "\n")
                    else:
                        build_seconds = None
                        build_peak_rss_mib = None
                        if build_metrics_path.exists():
                            try:
                                build_metrics = json.loads(build_metrics_path.read_text())
                                build_seconds = float(build_metrics["elapsed_seconds"])
                                build_peak_rss_mib = float(build_metrics["peak_rss_mib"])
                            except (KeyError, OSError, TypeError, ValueError):
                                pass
                    for candidate_k in candidates:
                        for ef_search in efs:
                            labels, distances, search_seconds, search_peak_rss_mib = run_cpp_search(
                                executable, index_path, config_dir / "similarity.float32.bin", query_codes_path,
                                len(queries), args.window_size, pq_m, nbits, candidate_k, ef_search,
                                args.threads, index_dir,
                            )
                            final_labels, final_scores, rerank_seconds = rerank_candidates_allow_missing(
                                database_windows, query_windows, labels, args.output_k, args.rerank_batch_size
                            )
                            row: dict[str, Any] = {
                                "backend": "hnswlib_cpp_node_pq",
                                "pq_m": pq_m,
                                "nbits": nbits,
                                "ksub": 1 << nbits,
                                "node_code_bytes": int(model["node_code_bytes"]),
                                "window_code_bytes": int(args.window_size * model["node_code_bytes"]),
                                "codebook_bytes": int(model["codebook_bytes"]),
                                "similarity_bytes": int(model["similarity_bytes"]),
                                "pq_train_seconds": model["train_seconds"],
                                "pq_encode_seconds": model["encode_seconds"],
                                "window_pack_seconds": pack_seconds,
                                "hnsw_m": hnsw_m,
                                "ef_construction": ef_construction,
                                "ef_search": ef_search,
                                "candidate_k": candidate_k,
                                "output_k": args.output_k,
                                "n_windows": window_count,
                                "n_queries": len(queries),
                                "search_seconds": search_seconds,
                                "exact_rerank_seconds": rerank_seconds,
                                "search_plus_rerank_seconds": search_seconds + rerank_seconds,
                                "build_seconds": build_seconds,
                                "build_peak_rss_mib": build_peak_rss_mib,
                                "search_peak_rss_mib": search_peak_rss_mib,
                                "index_bytes": index_path.stat().st_size,
                                "candidate_mean_score": float(np.mean(-distances[:, 0])),
                                "final_mean_score": float(np.mean(final_scores[:, 0])),
                                "index_path": str(index_path),
                            }
                            if reference is not None:
                                for k in (1, 5, 10, 50):
                                    row[f"candidate_recall_at_{k}"] = recall_at(reference, labels, k)
                                    row[f"final_recall_at_{k}"] = recall_at(reference, final_labels, k)
                            results.append(row)
                            print(json.dumps(row), flush=True)
    benchmark = {
        "backend": "hnswlib_cpp_node_pq",
        "embeddings": str(args.embeddings),
        "structures": str(args.structures),
        "query_selections": str(args.query_selections),
        "database_windows": str(args.database_windows),
        "queries": str(args.queries),
        "reference_labels": str(args.reference_labels) if args.reference_labels else None,
        "original_dtype": args.original_dtype,
        "window_size": args.window_size,
        "stride": args.stride,
        "node_count": node_count,
        "n_windows": window_count,
        "n_queries": len(queries),
        "query_stats": query_stats,
        "results": results,
    }
    (args.outdir / args.results_name).write_text(json.dumps(benchmark, indent=2) + "\n")
    return benchmark


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--query-selections", type=Path, required=True)
    parser.add_argument("--database-windows", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--reference-labels", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--results-name", default="results.json")
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--pq-m-values", default="4,8,16")
    parser.add_argument("--nbits-values", default="4,8")
    parser.add_argument("--hnsw-m-values", type=parse_int_list, default=[32])
    parser.add_argument("--ef-construction-values", type=parse_int_list, default=[200])
    parser.add_argument("--ef-search-values", default="200,800")
    parser.add_argument("--candidate-k-values", default="50,100")
    parser.add_argument("--output-k", type=int, default=50)
    parser.add_argument("--sample-size", type=int, default=250000)
    parser.add_argument("--pq-niter", type=int, default=20)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=32768)
    parser.add_argument("--rerank-batch-size", type=int, default=32)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--window-size", type=int, default=11)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--original-dtype", choices=("float16", "float32"), default="float16")
    parser.add_argument("--max-queries", type=int)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(evaluate(parse_args()), indent=2))
