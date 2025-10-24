#!/usr/bin/env python3
"""Calculate EVD parameters (lambda and K) for E-value estimation.

This script calculates database-wide Extreme Value Distribution parameters
that can be reused across multiple query alignments. The parameters are
estimated by performing random alignments on shuffled sequences and fitting
a Gumbel distribution to the resulting score distribution.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from multiprocessing import shared_memory
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd

NEG_INF = -1e9

TRACE_STOP = 0
TRACE_FROM_M = 1
TRACE_FROM_X = 2
TRACE_FROM_Y = 3
TRACE_GAP_FROM_M = 0
TRACE_GAP_FROM_SELF = 1

EMBEDDING_CHUNK_SIZE = 200_000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate EVD parameters for E-value estimation"
    )
    parser.add_argument(
        "--embeddings", required=True, help="Node-level embeddings TSV"
    )
    parser.add_argument(
        "--id-column",
        default="transcript_id",
        help="Identifier column in embeddings table",
    )
    parser.add_argument(
        "--node-index-column",
        default="node_index",
        help="Column giving node order within each transcript",
    )
    parser.add_argument(
        "--embedding-column",
        default="embedding_vector",
        help="Column containing comma-separated embedding vectors",
    )
    parser.add_argument(
        "--background-samples",
        type=int,
        default=10000,
        help="Number of random pairs to sample for background statistics",
    )
    parser.add_argument(
        "--evd-samples",
        type=int,
        default=1000,
        help="Number of random alignment samples for EVD parameter estimation",
    )
    parser.add_argument(
        "--gamma",
        type=float,
        default=1.5,
        help="Scaling factor for similarity scores",
    )
    parser.add_argument(
        "--band-width", type=int, default=96, help="Band width for banded alignment"
    )
    parser.add_argument(
        "--gap-open", type=float, default=12, help="Gap opening penalty"
    )
    parser.add_argument(
        "--gap-extend", type=float, default=2, help="Gap extension penalty"
    )
    parser.add_argument(
        "--xdrop", type=float, default=50, help="X-drop threshold for alignment"
    )
    parser.add_argument(
        "--score-min", type=float, default=-4, help="Minimum score clamp"
    )
    parser.add_argument(
        "--score-max", type=float, default=8, help="Maximum score clamp"
    )
    parser.add_argument(
        "--random-seed", type=int, default=42, help="Random seed for reproducibility"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Number of parallel workers (default: all available CPUs)",
    )
    parser.add_argument(
        "--sampled-sequences",
        type=int,
        default=50000,
        help="Number of sequences to randomly sample for EVD calculation (default: 50000, use 0 for no limit)",
    )
    parser.add_argument(
        "--output", required=True, help="Output JSON file with EVD parameters"
    )
    return parser.parse_args()


def collect_sequence_ids(path: str, id_col: str) -> List[str]:
    """Collect all unique sequence IDs from embeddings file.

    This is a lightweight first pass to get IDs without loading embeddings.
    """
    unique_ids = set()
    chunk_iter = pd.read_csv(
        path, sep="\t", chunksize=EMBEDDING_CHUNK_SIZE, dtype={id_col: str}, usecols=[id_col]
    )

    for chunk_df in chunk_iter:
        unique_ids.update(chunk_df[id_col].astype(str).unique())

    return list(unique_ids)


def sample_sequence_ids(all_ids: List[str], max_sequences: int, seed: int) -> Set[str]:
    """Randomly sample sequence IDs for EVD calculation.

    Args:
        all_ids: List of all available sequence IDs
        max_sequences: Maximum number to sample (0 = no limit)
        seed: Random seed for reproducibility

    Returns:
        Set of sampled sequence IDs (or all IDs if max_sequences=0 or >= len(all_ids))
    """
    if max_sequences <= 0 or max_sequences >= len(all_ids):
        return set(all_ids)

    rng = np.random.default_rng(seed)
    sampled = rng.choice(all_ids, size=max_sequences, replace=False)
    return set(sampled)


class SharedEmbeddings:
    """Wrapper for embeddings stored in shared memory.

    This allows multiple processes to access the same embeddings without
    copying the data, reducing memory usage from (1+workers)× to ~1×.
    """

    def __init__(self, shm_name: str, metadata: Dict):
        """Initialize accessor to shared memory embeddings.

        Args:
            shm_name: Name of the shared memory block
            metadata: Dict containing 'ids', 'offsets', 'lengths', 'embedding_dim', 'total_size'
        """
        self.shm_name = shm_name
        self.ids = metadata['ids']
        self.offsets = metadata['offsets']
        self.lengths = metadata['lengths']
        self.embedding_dim = metadata['embedding_dim']
        self.total_size = metadata['total_size']

        # Attach to existing shared memory
        self.shm = shared_memory.SharedMemory(name=shm_name)
        # Create numpy array view (no copy)
        self.flat_data = np.ndarray(
            (self.total_size,), dtype=np.float32, buffer=self.shm.buf
        )

    def __getitem__(self, seq_id: str) -> np.ndarray:
        """Get embeddings for a sequence ID.

        Returns a copy to avoid shared memory corruption.
        """
        if seq_id not in self.offsets:
            raise KeyError(f"Sequence ID not found: {seq_id}")

        offset = self.offsets[seq_id]
        length = self.lengths[seq_id]

        # Extract flat data
        flat = self.flat_data[offset:offset + length]

        # Reshape based on embedding dimension
        if self.embedding_dim == 0 or length % self.embedding_dim != 0:
            # Something went wrong, return as is
            return flat.copy()

        num_rows = length // self.embedding_dim
        data = flat.reshape(num_rows, self.embedding_dim)
        return data.copy()  # Return copy to avoid shared memory issues

    def __contains__(self, seq_id: str) -> bool:
        return seq_id in self.offsets

    def keys(self):
        return self.ids

    def items(self):
        for seq_id in self.ids:
            yield seq_id, self[seq_id]

    def __len__(self):
        return len(self.ids)

    def cleanup(self):
        """Close and unlink shared memory."""
        self.shm.close()

    def unlink(self):
        """Unlink shared memory (call from main process only)."""
        self.shm.unlink()


def embeddings_to_shared_memory(embeddings: Dict[str, np.ndarray]) -> Tuple[SharedEmbeddings, shared_memory.SharedMemory]:
    """Convert embeddings dictionary to shared memory format.

    Args:
        embeddings: Dictionary mapping sequence IDs to embedding arrays

    Returns:
        Tuple of (SharedEmbeddings accessor, SharedMemory object for cleanup)
    """
    if not embeddings:
        raise ValueError("Cannot create shared memory from empty embeddings")

    print(f"  Creating shared memory for {len(embeddings)} sequences...", file=sys.stderr, flush=True)

    # Determine embedding dimension - find first 2D array
    embedding_dim = None
    first_shape = None
    for emb in embeddings.values():
        first_shape = emb.shape
        if len(emb.shape) == 2:
            embedding_dim = emb.shape[1]
            break
        elif len(emb.shape) == 1:
            # This is a single-node sequence, dimension is the length
            embedding_dim = len(emb)
            break

    if embedding_dim is None:
        raise ValueError("Could not determine embedding dimension")

    print(f"  Embedding dimension: {embedding_dim}, first shape: {first_shape}", file=sys.stderr, flush=True)

    # Calculate total size needed and build metadata
    ids = []
    offsets = {}
    lengths = {}
    current_offset = 0

    for seq_id, emb in embeddings.items():
        ids.append(seq_id)
        offsets[seq_id] = current_offset

        # Ensure consistent shapes
        if len(emb.shape) == 1:
            # Single node: reshape to (1, embedding_dim)
            if len(emb) != embedding_dim:
                raise ValueError(f"Inconsistent embedding dimension for {seq_id}: expected {embedding_dim}, got {len(emb)}")
            emb_length = embedding_dim
        elif len(emb.shape) == 2:
            # Multiple nodes: (num_nodes, embedding_dim)
            if emb.shape[1] != embedding_dim:
                raise ValueError(f"Inconsistent embedding dimension for {seq_id}: expected {embedding_dim}, got {emb.shape[1]}")
            emb_length = emb.shape[0] * emb.shape[1]
        else:
            raise ValueError(f"Invalid embedding shape for {seq_id}: {emb.shape}")

        lengths[seq_id] = emb_length
        current_offset += emb_length

    total_size = current_offset
    memory_needed = total_size * 4  # 4 bytes per float32

    print(f"  Total memory needed: {memory_needed / (1024**2):.2f} MB", file=sys.stderr, flush=True)

    # Check if /dev/shm has enough space (important for Docker containers)
    try:
        import shutil
        shm_stats = shutil.disk_usage('/dev/shm')
        print(f"  /dev/shm available: {shm_stats.free / (1024**2):.2f} MB", file=sys.stderr, flush=True)
        if shm_stats.free < memory_needed * 1.1:  # Add 10% safety margin
            raise ValueError(
                f"Insufficient /dev/shm space: need {memory_needed / (1024**2):.2f} MB, "
                f"have {shm_stats.free / (1024**2):.2f} MB. "
                f"If running in Docker, increase --shm-size to at least {int(memory_needed / (1024**3) + 1)}g"
            )
    except ValueError:
        # Re-raise ValueError to be caught by caller
        raise
    except Exception as e:
        print(f"  Warning: Could not check /dev/shm space: {e}", file=sys.stderr, flush=True)

    # Create shared memory block
    try:
        shm = shared_memory.SharedMemory(create=True, size=memory_needed)
        print(f"  Created shared memory block '{shm.name}' ({memory_needed / (1024**2):.2f} MB)", file=sys.stderr, flush=True)
    except Exception as e:
        raise RuntimeError(f"Failed to create shared memory: {e}. If running in Docker, try increasing --shm-size")

    # Create numpy array view
    flat_data = np.ndarray((total_size,), dtype=np.float32, buffer=shm.buf)

    # Copy embeddings into shared memory
    print(f"  Copying {total_size} floats ({total_size * 4} bytes) to shared memory...", file=sys.stderr, flush=True)
    try:
        for i, seq_id in enumerate(ids):
            offset = offsets[seq_id]
            length = lengths[seq_id]
            emb = embeddings[seq_id]

            # Validate before copy
            if offset + length > total_size:
                raise ValueError(f"Offset overflow for {seq_id}: offset={offset}, length={length}, total={total_size}")

            # Flatten and copy
            flat_emb = emb.flatten()
            if len(flat_emb) != length:
                raise ValueError(f"Length mismatch for {seq_id}: expected {length}, got {len(flat_emb)}")

            flat_data[offset:offset + length] = flat_emb

            if (i + 1) % 20 == 0:
                print(f"    Copied {i + 1}/{len(ids)} sequences...", file=sys.stderr, flush=True)

        print(f"  Successfully copied all embeddings to shared memory", file=sys.stderr, flush=True)

    except Exception as e:
        shm.close()
        shm.unlink()
        raise RuntimeError(f"Failed to copy embeddings to shared memory: {e}")

    # Create metadata dict
    metadata = {
        'ids': ids,
        'offsets': offsets,
        'lengths': lengths,
        'embedding_dim': embedding_dim,
        'total_size': total_size,
    }

    # Create accessor
    print(f"  Creating SharedEmbeddings accessor...", file=sys.stderr, flush=True)
    try:
        accessor = SharedEmbeddings(shm.name, metadata)
        print(f"  SharedEmbeddings created successfully with name '{shm.name}'", file=sys.stderr, flush=True)
    except Exception as e:
        shm.close()
        shm.unlink()
        raise RuntimeError(f"Failed to create SharedEmbeddings accessor: {e}")

    return accessor, shm


def load_embeddings(
    path: str,
    id_col: str,
    node_idx_col: str,
    emb_col: str,
    required_ids: Optional[Set[str]] = None,
) -> Dict[str, np.ndarray]:
    """Load embeddings from TSV into a dictionary keyed by transcript ID."""
    embeddings: Dict[str, List[np.ndarray]] = {}
    chunk_iter = pd.read_csv(
        path, sep="\t", chunksize=EMBEDDING_CHUNK_SIZE, dtype={id_col: str}
    )

    for chunk_df in chunk_iter:
        chunk_df = chunk_df.sort_values([id_col, node_idx_col])
        for tid, group in chunk_df.groupby(id_col, sort=False):
            tid_str = str(tid)
            if required_ids is not None and tid_str not in required_ids:
                continue
            vectors = [
                np.fromstring(row[emb_col], sep=",", dtype=np.float32)
                for _, row in group.iterrows()
            ]
            if tid_str not in embeddings:
                embeddings[tid_str] = []
            embeddings[tid_str].extend(vectors)

    # Convert lists to numpy arrays
    final_embeddings = {tid: np.array(vec_list) for tid, vec_list in embeddings.items()}
    return final_embeddings


def _sample_background_batch(args_tuple):
    """Helper function for parallel background sampling.

    Args:
        args_tuple: Tuple containing (batch_idx, batch_size, embeddings, ids, base_seed)

    Returns:
        List of similarity scores for this batch
    """
    batch_idx, batch_size, embeddings, ids, base_seed = args_tuple

    # Create a unique RNG for this batch to ensure reproducibility
    rng = np.random.default_rng(base_seed + batch_idx)

    sims: List[float] = []
    for _ in range(batch_size):
        qid = rng.choice(ids)
        tid = rng.choice(ids)
        qvec = embeddings[qid][rng.integers(len(embeddings[qid]))]
        tvec = embeddings[tid][rng.integers(len(embeddings[tid]))]
        sims.append(float(np.dot(qvec, tvec)))

    return sims


def _sample_background_batch_shm(args_tuple):
    """Helper function for parallel background sampling using shared memory.

    Args:
        args_tuple: Tuple containing (batch_idx, batch_size, shm_name, metadata, ids, base_seed)

    Returns:
        List of similarity scores for this batch
    """
    batch_idx, batch_size, shm_name, metadata, ids, base_seed = args_tuple

    # Create shared embeddings accessor (no data copying)
    embeddings = SharedEmbeddings(shm_name, metadata)

    # Create a unique RNG for this batch to ensure reproducibility
    rng = np.random.default_rng(base_seed + batch_idx)

    sims: List[float] = []
    for _ in range(batch_size):
        qid = rng.choice(ids)
        tid = rng.choice(ids)
        qvec = embeddings[qid][rng.integers(len(embeddings[qid]))]
        tvec = embeddings[tid][rng.integers(len(embeddings[tid]))]
        sims.append(float(np.dot(qvec, tvec)))

    # Cleanup
    embeddings.cleanup()

    return sims


def sample_background(
    embeddings: Dict[str, np.ndarray], sample_count: int, seed: int, workers: Optional[int] = None, use_shared_memory: bool = True
) -> Tuple[float, float]:
    """Sample random node pairs to estimate background similarity distribution.

    Args:
        embeddings: Dictionary or SharedEmbeddings object containing sequence embeddings
        sample_count: Number of random pairs to sample
        seed: Random seed for reproducibility
        workers: Number of parallel workers (None = use all CPUs)
        use_shared_memory: If True, use shared memory to reduce memory overhead (default: True)

    Returns:
        Tuple of (mu, sigma) for background distribution
    """
    ids = [tid for tid, vec in embeddings.items() if len(vec) > 0]
    if not ids:
        raise SystemExit("No embeddings available for background sampling")

    # Determine number of workers
    if workers is None:
        workers = mp.cpu_count()
    workers = max(1, workers)

    # If single worker, no benefit from shared memory
    if workers == 1:
        use_shared_memory = False

    # Split work into batches for parallel processing
    batch_size = max(1, sample_count // workers)

    all_sims: List[float] = []

    if use_shared_memory:
        # Convert to shared memory if not already
        try:
            if isinstance(embeddings, SharedEmbeddings):
                shm_embeddings = embeddings
                shm = None  # Don't own the shared memory
            else:
                shm_embeddings, shm = embeddings_to_shared_memory(embeddings)
        except Exception as e:
            print(f"  WARNING: Shared memory initialization failed: {e}", file=sys.stderr, flush=True)
            print(f"  Falling back to standard multiprocessing (higher memory usage)", file=sys.stderr, flush=True)
            use_shared_memory = False

    if use_shared_memory:
        try:
            # Build metadata dict for workers
            metadata = {
                'ids': list(shm_embeddings.ids),
                'offsets': shm_embeddings.offsets,
                'lengths': shm_embeddings.lengths,
                'embedding_dim': shm_embeddings.embedding_dim,
                'total_size': shm_embeddings.total_size,
            }

            args_list = [
                (i, batch_size if i < workers - 1 else sample_count - i * batch_size, shm_embeddings.shm_name, metadata, ids, seed)
                for i in range(workers)
            ]

            # Run background sampling in parallel with shared memory
            with mp.Pool(processes=workers) as pool:
                for batch_sims in pool.imap_unordered(_sample_background_batch_shm, args_list):
                    all_sims.extend(batch_sims)

        finally:
            # Cleanup shared memory if we created it
            if shm is not None:
                shm_embeddings.cleanup()
                shm.unlink()
    else:
        # Original implementation without shared memory
        args_list = [
            (i, batch_size if i < workers - 1 else sample_count - i * batch_size, embeddings, ids, seed)
            for i in range(workers)
        ]

        with mp.Pool(processes=workers) as pool:
            for batch_sims in pool.imap_unordered(_sample_background_batch, args_list):
                all_sims.extend(batch_sims)

    mu = float(np.mean(all_sims))
    sigma = float(np.std(all_sims))
    if sigma < 1e-6:
        sigma = 1.0
    return mu, sigma


def compute_score(
    dot: float, mu: float, sigma: float, gamma: float, smin: float, smax: float
) -> float:
    """Convert dot product to scaled score."""
    z = gamma * ((dot - mu) / sigma)
    return float(np.clip(z, smin, smax))


def shuffle_embeddings(embeddings: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle the order of embedding vectors to create a random sequence.

    This preserves the composition (same vectors) but destroys sequential information.

    Args:
        embeddings: Array of shape (n_nodes, embedding_dim)
        rng: Random number generator

    Returns:
        Shuffled copy of embeddings
    """
    shuffled = embeddings.copy()
    rng.shuffle(shuffled)
    return shuffled


def smith_waterman(
    query: np.ndarray,
    target: np.ndarray,
    mu: float,
    sigma: float,
    gamma: float,
    band: int,
    diag_offset: int,
    gap_open: float,
    gap_extend: float,
    xdrop: float,
    clamp_min: float,
    clamp_max: float,
    store_matrix: bool = False,
) -> Tuple[float, List[Tuple[str, int, int]], Optional[np.ndarray]]:
    """Banded Smith-Waterman alignment with affine gaps.

    This is a simplified version for EVD parameter estimation that only
    needs to return the alignment score.

    Args:
        query: Query embedding array (m x d)
        target: Target embedding array (n x d)
        mu: Background mean for scoring
        sigma: Background std for scoring
        gamma: Scaling factor for scores
        band: Band width
        diag_offset: Diagonal offset for band center
        gap_open: Gap opening penalty
        gap_extend: Gap extension penalty
        xdrop: X-drop threshold
        clamp_min: Minimum score clamp
        clamp_max: Maximum score clamp
        store_matrix: Whether to store full DP matrix (unused for EVD)

    Returns:
        Tuple of (score, path, matrix)
    """
    m, n = len(query), len(target)
    if m == 0 or n == 0:
        return 0.0, [], None

    band = max(1, int(band))
    band_half = band // 2

    M_prev = np.zeros(n + 1, dtype=np.float32)
    X_prev = np.full(n + 1, NEG_INF, dtype=np.float32)
    Y_prev = np.full(n + 1, NEG_INF, dtype=np.float32)
    M_curr = np.zeros(n + 1, dtype=np.float32)
    X_curr = np.full(n + 1, NEG_INF, dtype=np.float32)
    Y_curr = np.full(n + 1, NEG_INF, dtype=np.float32)

    band_capacity = min(n, band + 1) + 2
    trace_M = np.zeros((m + 1, band_capacity), dtype=np.uint8)
    trace_X = np.zeros((m + 1, band_capacity), dtype=np.uint8)
    trace_Y = np.zeros((m + 1, band_capacity), dtype=np.uint8)
    row_offsets = np.zeros(m + 1, dtype=np.int32)
    row_lengths = np.zeros(m + 1, dtype=np.int32)

    matrix_scores: Optional[np.ndarray] = None
    if store_matrix:
        matrix_scores = np.full((m + 1, n + 1), np.nan, dtype=np.float32)
        matrix_scores[0, :] = 0.0
        matrix_scores[:, 0] = 0.0

    best_score = 0.0
    best_state = ("M", 0, 0)

    for i in range(1, m + 1):
        M_curr.fill(NEG_INF)
        X_curr.fill(NEG_INF)
        Y_curr.fill(NEG_INF)
        M_curr[0] = 0.0

        row_best = NEG_INF
        center = i + diag_offset
        j_start = max(1, center - band_half)
        j_end = min(n, center + band_half)
        if j_start > j_end:
            row_offsets[i] = 1
            row_lengths[i] = 0
            continue

        row_offsets[i] = j_start
        span = j_end - j_start + 1
        row_lengths[i] = span

        for j in range(j_start, j_end + 1):
            i_global = i - 1
            j_global = j - 1
            local_idx = j - j_start
            if local_idx < 0 or local_idx >= band_capacity:
                continue
            dot = float(np.dot(query[i_global], target[j_global]))
            score = compute_score(dot, mu, sigma, gamma, clamp_min, clamp_max)

            # Match/mismatch state
            candidates = [
                (0.0, None),
                (float(M_prev[j - 1]), ("M", i - 1, j - 1)),
                (float(X_prev[j - 1]), ("X", i - 1, j - 1)),
                (float(Y_prev[j - 1]), ("Y", i - 1, j - 1)),
            ]
            base_val, base_ptr = max(candidates, key=lambda item: item[0])
            match_val = base_val + score
            M_curr[j] = match_val
            if base_ptr is None:
                trace_M[i, local_idx] = TRACE_STOP
            else:
                prev_state = base_ptr[0]
                if prev_state == "M":
                    trace_M[i, local_idx] = TRACE_FROM_M
                elif prev_state == "X":
                    trace_M[i, local_idx] = TRACE_FROM_X
                else:
                    trace_M[i, local_idx] = TRACE_FROM_Y

            # Gap in query (insertion in target)
            left_m = M_curr[j - 1] - gap_open
            left_x = X_curr[j - 1] - gap_extend
            if left_m >= left_x:
                X_curr[j] = left_m
                trace_X[i, local_idx] = TRACE_GAP_FROM_M
            else:
                X_curr[j] = left_x
                trace_X[i, local_idx] = TRACE_GAP_FROM_SELF

            # Gap in target (insertion in query)
            up_m = M_prev[j] - gap_open
            up_y = Y_prev[j] - gap_extend
            if up_m >= up_y:
                Y_curr[j] = up_m
                trace_Y[i, local_idx] = TRACE_GAP_FROM_M
            else:
                Y_curr[j] = up_y
                trace_Y[i, local_idx] = TRACE_GAP_FROM_SELF

            # Track best score
            curr_best = max(M_curr[j], X_curr[j], Y_curr[j])
            if curr_best > best_score:
                best_score = curr_best
                if M_curr[j] >= X_curr[j] and M_curr[j] >= Y_curr[j]:
                    best_state = ("M", i, j)
                elif X_curr[j] >= Y_curr[j]:
                    best_state = ("X", i, j)
                else:
                    best_state = ("Y", i, j)

            if curr_best > row_best:
                row_best = curr_best

            if store_matrix:
                matrix_scores[i, j] = curr_best

        # X-drop pruning
        if best_score - row_best > xdrop and best_score > 0:
            # Continue for a few more rows to ensure we don't prematurely stop
            pass

        M_prev, M_curr = M_curr, M_prev
        X_prev, X_curr = X_curr, X_prev
        Y_prev, Y_curr = Y_curr, Y_prev

    # For EVD estimation, we only need the score, not the traceback
    return float(best_score), [], matrix_scores


def _align_random_pair(args_tuple):
    """Helper function for parallel alignment of random sequence pairs.

    Args:
        args_tuple: Tuple containing (sample_idx, embeddings_dict, seq_ids, mu, sigma,
                    gamma, band, gap_open, gap_extend, xdrop, clamp_min, clamp_max, seed)

    Returns:
        Alignment score (float) or None if score <= 0
    """
    (
        sample_idx,
        embeddings,
        seq_ids,
        mu,
        sigma,
        gamma,
        band,
        gap_open,
        gap_extend,
        xdrop,
        clamp_min,
        clamp_max,
        base_seed,
    ) = args_tuple

    # Create a unique RNG for this sample to ensure reproducibility
    rng = np.random.default_rng(base_seed + sample_idx)

    # Sample two random sequences
    id1, id2 = rng.choice(seq_ids, size=2, replace=True)

    # Get embeddings and shuffle them
    emb1 = shuffle_embeddings(embeddings[id1], rng)
    emb2 = shuffle_embeddings(embeddings[id2], rng)

    # Align with Smith-Waterman
    score, _, _ = smith_waterman(
        emb1,
        emb2,
        mu,
        sigma,
        gamma,
        band,
        diag_offset=0,  # Random offset for null model
        gap_open=gap_open,
        gap_extend=gap_extend,
        xdrop=xdrop,
        clamp_min=clamp_min,
        clamp_max=clamp_max,
        store_matrix=False,
    )

    return float(score) if score > 0 else None


def _align_random_pair_shm(args_tuple):
    """Helper function for parallel alignment of random sequence pairs using shared memory.

    Args:
        args_tuple: Tuple containing (sample_idx, shm_name, metadata, seq_ids, mu, sigma,
                    gamma, band, gap_open, gap_extend, xdrop, clamp_min, clamp_max, seed)

    Returns:
        Alignment score (float) or None if score <= 0
    """
    (
        sample_idx,
        shm_name,
        metadata,
        seq_ids,
        mu,
        sigma,
        gamma,
        band,
        gap_open,
        gap_extend,
        xdrop,
        clamp_min,
        clamp_max,
        base_seed,
    ) = args_tuple

    # Create shared embeddings accessor (no data copying)
    embeddings = SharedEmbeddings(shm_name, metadata)

    # Create a unique RNG for this sample to ensure reproducibility
    rng = np.random.default_rng(base_seed + sample_idx)

    # Sample two random sequences
    id1, id2 = rng.choice(seq_ids, size=2, replace=True)

    # Get embeddings and shuffle them
    emb1 = shuffle_embeddings(embeddings[id1], rng)
    emb2 = shuffle_embeddings(embeddings[id2], rng)

    # Align with Smith-Waterman
    score, _, _ = smith_waterman(
        emb1,
        emb2,
        mu,
        sigma,
        gamma,
        band,
        diag_offset=0,  # Random offset for null model
        gap_open=gap_open,
        gap_extend=gap_extend,
        xdrop=xdrop,
        clamp_min=clamp_min,
        clamp_max=clamp_max,
        store_matrix=False,
    )

    # Cleanup
    embeddings.cleanup()

    return float(score) if score > 0 else None


def estimate_evd_parameters(
    embeddings: Dict[str, np.ndarray],
    mu: float,
    sigma: float,
    gamma: float,
    band: int,
    gap_open: float,
    gap_extend: float,
    xdrop: float,
    clamp_min: float,
    clamp_max: float,
    num_samples: int = 1000,
    seed: int = 42,
    workers: Optional[int] = None,
    use_shared_memory: bool = False,
) -> Tuple[float, float]:
    """Estimate lambda and K parameters for BLAST-like E-value calculation.

    Generates random sequence pairs by shuffling, aligns them to get a null
    distribution of scores, and fits an Extreme Value Distribution (Gumbel).

    Args:
        embeddings: Dictionary mapping sequence IDs to embedding arrays
        mu: Background mean for scoring
        sigma: Background std for scoring
        gamma: Scaling factor for scores
        band: Band width for alignment
        gap_open: Gap opening penalty
        gap_extend: Gap_extension penalty
        xdrop: X-drop threshold
        clamp_min: Minimum score clamp
        clamp_max: Maximum score clamp
        num_samples: Number of random pairs to generate
        seed: Random seed
        workers: Number of parallel workers (None = use all CPUs)
        use_shared_memory: If True, use shared memory to reduce memory overhead (default: True)

    Returns:
        Tuple of (lambda, K) parameters
    """
    from scipy.stats import gumbel_r

    # Get list of sequences that have sufficient length
    seq_ids = [tid for tid, vec in embeddings.items() if len(vec) >= 20]
    if len(seq_ids) < 2:
        # Not enough data, return default values
        print(
            "Warning: Not enough sequences for EVD estimation. Using default parameters.",
            file=sys.stderr,
            flush=True,
        )
        return 0.1, 0.01

    # Determine number of workers
    if workers is None:
        workers = mp.cpu_count()
    workers = max(1, workers)

    # If single worker, no benefit from shared memory
    if workers == 1:
        use_shared_memory = False

    print(
        f"Estimating EVD parameters from {num_samples} random alignments using {workers} workers...",
        file=sys.stderr,
        flush=True,
    )

    null_scores: List[float] = []

    if use_shared_memory:
        # Convert to shared memory if not already
        if isinstance(embeddings, SharedEmbeddings):
            shm_embeddings = embeddings
            shm = None  # Don't own the shared memory
        else:
            print("  Converting embeddings to shared memory...", file=sys.stderr, flush=True)
            shm_embeddings, shm = embeddings_to_shared_memory(embeddings)

        try:
            # Build metadata dict for workers
            metadata = {
                'ids': list(shm_embeddings.ids),
                'offsets': shm_embeddings.offsets,
                'lengths': shm_embeddings.lengths,
                'embedding_dim': shm_embeddings.embedding_dim,
                'total_size': shm_embeddings.total_size,
            }

            # Prepare arguments for parallel processing
            args_list = [
                (
                    i,
                    shm_embeddings.shm_name,
                    metadata,
                    seq_ids,
                    mu,
                    sigma,
                    gamma,
                    band,
                    gap_open,
                    gap_extend,
                    xdrop,
                    clamp_min,
                    clamp_max,
                    seed,
                )
                for i in range(num_samples)
            ]

            # Run alignments in parallel with shared memory
            with mp.Pool(processes=workers) as pool:
                # Use imap_unordered for better memory efficiency and progress tracking
                for i, score in enumerate(pool.imap_unordered(_align_random_pair_shm, args_list, chunksize=10)):
                    if score is not None:
                        null_scores.append(score)

                    if (i + 1) % 100 == 0:
                        print(
                            f"  Completed {i + 1}/{num_samples} random alignments...",
                            file=sys.stderr,
                            flush=True,
                        )

        finally:
            # Cleanup shared memory if we created it
            if shm is not None:
                shm_embeddings.cleanup()
                shm.unlink()
    else:
        # Original implementation without shared memory
        # Prepare arguments for parallel processing
        args_list = [
            (
                i,
                embeddings,
                seq_ids,
                mu,
                sigma,
                gamma,
                band,
                gap_open,
                gap_extend,
                xdrop,
                clamp_min,
                clamp_max,
                seed,
            )
            for i in range(num_samples)
        ]

        # Run alignments in parallel
        with mp.Pool(processes=workers) as pool:
            # Use imap_unordered for better memory efficiency and progress tracking
            for i, score in enumerate(pool.imap_unordered(_align_random_pair, args_list, chunksize=10)):
                if score is not None:
                    null_scores.append(score)

                if (i + 1) % 100 == 0:
                    print(
                        f"  Completed {i + 1}/{num_samples} random alignments...",
                        file=sys.stderr,
                        flush=True,
                    )

    if len(null_scores) < 10:
        print(
            "Warning: Very few random alignments produced positive scores. Using default parameters.",
            file=sys.stderr,
            flush=True,
        )
        return 0.1, 0.01

    # Fit Gumbel distribution (right-skewed EVD)
    # The Gumbel distribution parameters are (loc=mu, scale=beta)
    fit_mu, fit_beta = gumbel_r.fit(null_scores)

    # Convert to Karlin-Altschul parameters
    lambda_param = 1.0 / fit_beta
    K = np.exp(-lambda_param * fit_mu)

    print(f"EVD fitting complete:", file=sys.stderr, flush=True)
    print(f"  Gumbel location (μ) = {fit_mu:.4f}", file=sys.stderr, flush=True)
    print(f"  Gumbel scale (β) = {fit_beta:.4f}", file=sys.stderr, flush=True)
    print(f"  Lambda (λ) = {lambda_param:.4f}", file=sys.stderr, flush=True)
    print(f"  K = {K:.6f}", file=sys.stderr, flush=True)
    print(f"  Mean null score = {np.mean(null_scores):.4f}", file=sys.stderr, flush=True)
    print(f"  Std null score = {np.std(null_scores):.4f}", file=sys.stderr, flush=True)

    return float(lambda_param), float(K)


def main() -> None:
    args = parse_args()

    # Collect all sequence IDs and optionally sample them
    sampled_ids = None
    total_sequences = None
    if args.sampled_sequences > 0:
        print("Collecting sequence IDs...", file=sys.stderr, flush=True)
        all_ids = collect_sequence_ids(args.embeddings, args.id_column)
        total_sequences = len(all_ids)
        print(f"Found {total_sequences} total sequences in database", file=sys.stderr, flush=True)

        sampled_ids = sample_sequence_ids(all_ids, args.sampled_sequences, args.random_seed)
        if len(sampled_ids) < total_sequences:
            print(
                f"Randomly sampling {len(sampled_ids)} sequences for EVD calculation",
                file=sys.stderr,
                flush=True,
            )
        else:
            print("Using all sequences (no sampling needed)", file=sys.stderr, flush=True)

    print("Loading embeddings...", file=sys.stderr, flush=True)
    embeddings = load_embeddings(
        args.embeddings,
        args.id_column,
        args.node_index_column,
        args.embedding_column,
        required_ids=sampled_ids,
    )

    print(
        f"Loaded {len(embeddings)} sequences with embeddings", file=sys.stderr, flush=True
    )

    print("Sampling background statistics...", file=sys.stderr, flush=True)
    mu0, sigma0 = sample_background(embeddings, args.background_samples, args.random_seed, workers=args.workers)
    print(f"  Background μ = {mu0:.4f}", file=sys.stderr, flush=True)
    print(f"  Background σ = {sigma0:.4f}", file=sys.stderr, flush=True)

    # Estimate EVD parameters
    lambda_param, K_param = estimate_evd_parameters(
        embeddings,
        mu0,
        sigma0,
        args.gamma,
        args.band_width,
        args.gap_open,
        args.gap_extend,
        args.xdrop,
        args.score_min,
        args.score_max,
        num_samples=args.evd_samples,
        seed=args.random_seed,
        workers=args.workers,
    )

    # Write output JSON
    output_data = {
        "evd_lambda": float(lambda_param),
        "evd_K": float(K_param),
        "background_mu": float(mu0),
        "background_sigma": float(sigma0),
        "gamma": args.gamma,
        "band_width": args.band_width,
        "gap_open": args.gap_open,
        "gap_extend": args.gap_extend,
        "xdrop": args.xdrop,
        "score_min": args.score_min,
        "score_max": args.score_max,
        "evd_samples": args.evd_samples,
        "background_samples": args.background_samples,
        "random_seed": args.random_seed,
        "num_sequences": len(embeddings),
        "sampled_sequences": args.sampled_sequences,
        "total_sequences_in_database": total_sequences if total_sequences is not None else len(embeddings),
        "sequences_sampled": total_sequences is not None and len(embeddings) < total_sequences,
    }

    Path(args.output).write_text(json.dumps(output_data, indent=2) + "\n")

    print(f"\nEVD parameters written to {args.output}", file=sys.stderr, flush=True)
    print(json.dumps(output_data, indent=2))


if __name__ == "__main__":
    main()
