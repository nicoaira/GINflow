#!/usr/bin/env python3
"""Build configurable FAISS indices for window vectors."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np


def _ensure_numpy_faiss_compat() -> None:
    version_parts = tuple(int(part) for part in np.__version__.split('.')[:2])
    if version_parts >= (2, 0):
        raise SystemExit(
            "Detected NumPy "
            f"{np.__version__}"
            " – FAISS Python wheels bundled with the pipeline target NumPy < 2. "
            "Re-run with the 'conda' or 'docker' profile, or install numpy<2 in the active environment."
        )


_ensure_numpy_faiss_compat()

import faiss
import pandas as pd

GPU_MEMORY_ERRORS: tuple[type[BaseException], ...] = ()
if hasattr(faiss, "GpuMemoryException"):
    GPU_MEMORY_ERRORS = GPU_MEMORY_ERRORS + (faiss.GpuMemoryException,)

OOM_EXCEPTIONS: tuple[type[BaseException], ...] = (MemoryError,) + GPU_MEMORY_ERRORS


def log(message: str) -> None:
    """Emit progress logs to stderr without affecting JSON output."""
    print(message, file=sys.stderr, flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a FAISS index over window vectors")
    parser.add_argument("--vectors", required=True, help="Path to .npy array with database window vectors")
    parser.add_argument("--metadata", required=True, help="TSV mapping rows to transcript/position metadata")
    parser.add_argument("--index-type", default="flat_ip",
                        choices=["flat_ip", "flat_l2", "ivf", "ivfpq", "opq_ivfpq", "hnsw", "hnswsq8"],
                        help="FAISS index variant to construct")
    parser.add_argument("--metric", default="ip", choices=["ip", "l2"],
                        help="Distance metric for the index (inner-product or L2)")
    parser.add_argument("--nlist", type=int, default=1024, help="Number of coarse clusters for IVF-based indices")
    parser.add_argument("--pq-m", type=int, default=32, help="Number of subquantizers for (I)VFPQ indices")
    parser.add_argument("--pq-bits", type=int, default=8, help="Bits per subquantizer for (I)VFPQ indices")
    parser.add_argument("--opq-m", type=int, default=64, help="OPQ matrix size for opq_ivfpq")
    parser.add_argument("--hnsw-m", type=int, default=32, help="Connectivity parameter for HNSW")
    parser.add_argument("--hnsw-efc", type=int, default=200, help="efConstruction for HNSW")
    parser.add_argument("--output-index", required=True, help="Destination path for the FAISS index")
    parser.add_argument("--output-mapping", required=True, help="Destination path for the metadata TSV copy")
    parser.add_argument("--stats-json", default=None, help="Optional JSON file with index metadata")
    parser.add_argument("--use-gpu", action="store_true", help="Move the index build to GPU when available")
    parser.add_argument(
        "--fallback-index-type",
        default="ivfpq",
        choices=["flat_ip", "flat_l2", "ivf", "ivfpq", "opq_ivfpq", "hnsw", "hnswsq8", "none"],
        help="Index type to retry with if a flat index runs out of memory (use 'none' to disable)",
    )
    parser.add_argument(
        "--add-batch-size",
        type=int,
        default=0,
        help="Number of vectors to add per FAISS add() call; 0 loads everything at once",
    )
    parser.add_argument(
        "--train-size",
        type=int,
        default=250_000,
        help="Maximum number of vectors to sample for IVF/IVFPQ training (0 uses all vectors)",
    )
    parser.add_argument(
        "--train-seed",
        type=int,
        default=0,
        help="RNG seed for selecting the training subset",
    )
    return parser.parse_args()


def load_vectors(path: str) -> np.ndarray:
    vectors = np.load(path, mmap_mode='r')
    if vectors.ndim != 2:
        raise SystemExit(f"Expected 2D array in {path}, found shape {vectors.shape}")
    if vectors.dtype != np.float32:
        vectors = vectors.astype(np.float32)
    return vectors


def metric_id(metric: str) -> int:
    if metric == "ip":
        return faiss.METRIC_INNER_PRODUCT
    return faiss.METRIC_L2


def ensure_unit_norm(vectors: np.ndarray) -> None:
    """Warn if vectors are not approximately unit-length (for cosine/IP searches)."""
    norms = np.linalg.norm(vectors[: min(len(vectors), 1000)], axis=1)
    avg = float(norms.mean()) if len(norms) else 0.0
    if not np.isclose(avg, 1.0, atol=1e-2):
        print(f"WARNING: average vector norm ≈ {avg:.3f}; cosine/IP assumptions may be violated", file=sys.stderr)


def _base_index_for_gpu_options(index: faiss.Index) -> faiss.Index:
    """Return the underlying FAISS index when wrapped in pre/post transforms."""
    if hasattr(faiss, "IndexPreTransform") and isinstance(index, faiss.IndexPreTransform):
        return faiss.downcast_index(index.index)
    return index


def maybe_to_gpu(index: faiss.Index, use_gpu: bool) -> faiss.Index:
    if not use_gpu:
        return index
    if not hasattr(faiss, "StandardGpuResources"):
        print("WARNING: FAISS build lacks GPU support; continuing on CPU", file=sys.stderr)
        return index
    res = faiss.StandardGpuResources()
    options = None
    if hasattr(faiss, "GpuClonerOptions"):
        options = faiss.GpuClonerOptions()
        base_index = _base_index_for_gpu_options(index)
        if hasattr(faiss, "IndexIVFPQ") and isinstance(base_index, faiss.IndexIVFPQ):
            if hasattr(options, "useFloat16LookupTables"):
                options.useFloat16LookupTables = True
    try:
        if options is not None:
            return faiss.index_cpu_to_gpu(res, 0, index, options)
        return faiss.index_cpu_to_gpu(res, 0, index)
    except RuntimeError as err:
        message = str(err)
        if "shared memory" in message and "requires" in message:
            print(
                "WARNING: GPU index initialisation failed due to shared memory limits; continuing on CPU",
                file=sys.stderr,
            )
            return index
        raise


def select_training_vectors(vectors: np.ndarray, limit: int, seed: int) -> np.ndarray:
    if limit <= 0 or len(vectors) <= limit:
        return np.asarray(vectors, dtype=np.float32)
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(len(vectors), size=limit, replace=False))
    return np.asarray(vectors[idx], dtype=np.float32)


def build_index(args: argparse.Namespace, vectors: np.ndarray) -> tuple[faiss.Index, str]:
    dim = vectors.shape[1]
    metric = metric_id(args.metric)

    if args.metric == "ip":
        ensure_unit_norm(vectors)

    train_vectors = None
    requires_training = args.index_type in {"ivf", "ivfpq", "opq_ivfpq", "hnswsq8"}
    if requires_training:
        train_vectors = select_training_vectors(vectors, args.train_size, args.train_seed)
        if len(train_vectors) == 0:
            raise SystemExit("Training set empty – cannot build trained index")
        log(
            f"Prepared {len(train_vectors)} training vectors "
            f"(limit {args.train_size or 'all'})"
        )

    if args.index_type == "flat_ip":
        index = faiss.IndexFlatIP(dim)
    elif args.index_type == "flat_l2":
        index = faiss.IndexFlatL2(dim)
    elif args.index_type == "ivf":
        quantizer = faiss.IndexFlatIP(dim) if args.metric == "ip" else faiss.IndexFlatL2(dim)
        index = faiss.IndexIVFFlat(quantizer, dim, args.nlist, metric)
    elif args.index_type == "ivfpq":
        quantizer = faiss.IndexFlatIP(dim) if args.metric == "ip" else faiss.IndexFlatL2(dim)
        index = faiss.IndexIVFPQ(quantizer, dim, args.nlist, args.pq_m, args.pq_bits, metric)
    elif args.index_type == "opq_ivfpq":
        opq = faiss.OPQMatrix(dim, args.opq_m)
        quantizer = faiss.IndexFlatIP(dim) if args.metric == "ip" else faiss.IndexFlatL2(dim)
        base = faiss.IndexIVFPQ(quantizer, dim, args.nlist, args.pq_m, args.pq_bits, metric)
        index = faiss.IndexPreTransform(opq, base)
    elif args.index_type == "hnsw":
        index = faiss.IndexHNSWFlat(dim, args.hnsw_m)
        index.metric_type = metric
        index.hnsw.efConstruction = args.hnsw_efc
    elif args.index_type == "hnswsq8":
        if not hasattr(faiss, "IndexHNSWSQ"):
            raise SystemExit("FAISS build lacks IndexHNSWSQ required for hnswsq8")
        if not hasattr(faiss, "ScalarQuantizer"):
            raise SystemExit("FAISS build lacks ScalarQuantizer support required for hnswsq8")
        index = faiss.IndexHNSWSQ(dim, faiss.ScalarQuantizer.QT_8bit, args.hnsw_m)
        index.metric_type = metric
        index.hnsw.efConstruction = args.hnsw_efc
    else:  # pragma: no cover - guarded by argparse choices
        raise SystemExit(f"Unsupported index type: {args.index_type}")

    cpu_index = index
    index = maybe_to_gpu(index, args.use_gpu)
    if requires_training:
        device = "GPU" if index is not cpu_index and args.use_gpu else "CPU"
        log(
            f"Training {args.index_type} index on {device} "
            f"with {len(train_vectors)} vectors"
        )
        train_start = time.perf_counter()
        index.train(train_vectors)
        train_elapsed = time.perf_counter() - train_start
        log(
            f"Finished training {args.index_type} index in {train_elapsed:.2f}s"
        )
    try:
        add_vectors(index, vectors, args.add_batch_size)
    except OOM_EXCEPTIONS as err:
        if args.index_type in {"flat_ip", "flat_l2"} and args.fallback_index_type != "none":
            is_gpu_index = index is not cpu_index
            is_gpu_oom = isinstance(err, GPU_MEMORY_ERRORS) or (is_gpu_index and args.use_gpu)
            resource = "GPU VRAM" if is_gpu_oom else "host RAM"
            log(
                f"WARNING: flat index ran out of {resource}; retrying with {args.fallback_index_type}"
            )
            fallback_args = argparse.Namespace(**vars(args))
            fallback_args.index_type = args.fallback_index_type
            fallback_args.fallback_index_type = "none"
            return build_index(fallback_args, vectors)
        raise err
    if index is not cpu_index:
        index = faiss.index_gpu_to_cpu(index)
    return index, args.index_type


def add_vectors(index: faiss.Index, vectors: np.ndarray, batch_size: int) -> None:
    """Add vectors to the index in RAM-friendly chunks."""
    total = vectors.shape[0]
    if batch_size <= 0:
        log(f"Adding all {total} vectors in a single batch")
        index.add(np.asarray(vectors, dtype=np.float32))
        log("Finished adding vectors to index")
        return

    batches = (total + batch_size - 1) // batch_size
    log(f"Adding {total} vectors in {batches} batches (batch size {batch_size})")
    processed = 0
    for start in range(0, total, batch_size):
        end = min(start + batch_size, total)
        chunk = np.asarray(vectors[start:end], dtype=np.float32)
        index.add(chunk)
        processed = end
        log(f"  Added batch {(start // batch_size) + 1}/{batches} ({processed}/{total} vectors)")
    log("Finished adding vectors to index")


def main() -> None:
    args = parse_args()
    log("Starting FAISS index build")
    log(f"Loading vectors from {args.vectors}")
    vectors = load_vectors(args.vectors)
    metadata = Path(args.metadata)
    if not metadata.exists():
        raise SystemExit(f"Metadata file {metadata} not found")
    log(f"Loaded {vectors.shape[0]} vectors (dimension {vectors.shape[1]})")
    log(f"Reading metadata from {metadata}")

    if vectors.shape[0] == 0:
        raise SystemExit("No database windows available – aborting index build")

    log(
        "Building FAISS index with parameters: "
        f"type={args.index_type}, metric={args.metric}, use_gpu={args.use_gpu}"
    )
    index, effective_type = build_index(args, vectors)
    log(f"Constructed {effective_type} index; writing outputs")
    log(f"Writing FAISS index to {args.output_index}")
    faiss.write_index(index, args.output_index)

    # Copy metadata alongside the index for downstream lookups.
    log(f"Copying metadata to {args.output_mapping}")
    metadata_df = pd.read_csv(metadata, sep='\t')
    metadata_df.to_csv(args.output_mapping, sep='\t', index=False)

    stats = {
        "index_type": effective_type,
        "metric": args.metric,
        "vectors": int(vectors.shape[0]),
        "dimension": int(vectors.shape[1]),
    }
    if args.stats_json:
        log(f"Writing index statistics to {args.stats_json}")
        Path(args.stats_json).write_text(json.dumps(stats, indent=2) + "\n")

    log("FAISS index build complete")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
