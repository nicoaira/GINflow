#!/usr/bin/env python3
"""Run the resumable 30k Rouskin PQ-HNSW/CAGRA benchmark matrix.

This driver is deliberately a thin subprocess orchestrator.  It writes one
log and one result directory per configuration, so an interrupted overnight
run can be resumed without rebuilding completed PQ/HNSW artifacts.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run(command: list[str], log_path: Path, env: dict[str, str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as log:
        log.write("\n$ " + " ".join(command) + "\n")
        log.flush()
        completed = subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}); see {log_path}")


def run_pq(args: argparse.Namespace, env: dict[str, str]) -> None:
    output = args.outdir / "pq-hnsw"
    result = output / args.results_name
    if result.exists() and not args.force:
        return
    command = [
        sys.executable,
        str(ROOT / "bin/benchmark_pq_hnswlib.py"),
        "--embeddings", str(args.embeddings),
        "--structures", str(args.structures),
        "--query-selections", str(args.query_selections),
        "--database-windows", str(args.database_windows),
        "--queries", str(args.queries),
        "--reference-labels", str(args.reference_labels),
        "--outdir", str(output),
        "--results-name", args.results_name,
        "--executable", str(args.pq_executable),
        "--pq-m-values", args.pq_m_values,
        "--nbits-values", args.nbits_values,
        "--hnsw-m-values", args.hnsw_m_values,
        "--ef-construction-values", args.ef_construction_values,
        "--ef-search-values", args.ef_search_values,
        "--candidate-k-values", args.candidate_k_values,
        "--output-k", str(args.output_k),
        "--recall-k-values", "100,500",
        "--sample-size", str(args.sample_size),
        "--pq-niter", str(args.pq_niter),
        "--batch-size", str(args.encode_batch_size),
        "--rerank-batch-size", str(args.rerank_batch_size),
        "--candidate-batch-size", str(args.candidate_batch_size),
        "--rerank-workers", str(args.rerank_workers),
        "--rerank-device", "cpu",
        "--threads", str(args.hnsw_threads),
        "--window-size", str(args.window_size),
        "--stride", str(args.stride),
        "--original-dtype", "float16",
    ]
    run(command, args.outdir / "pq-hnsw.log", env)


def run_cagra(args: argparse.Namespace, env: dict[str, str]) -> None:
    if not args.cagra:
        return
    for candidate in [int(value) for value in args.candidate_k_values.split(",") if value.strip()]:
        output = args.outdir / f"cagra-int8-k{candidate}"
        result = output / "result.json"
        if result.exists() and not args.force:
            continue
        command = [
            sys.executable,
            str(ROOT / "bin/benchmark_cagra_gpu.py"),
            "--data", str(args.database_windows),
            "--queries", str(args.queries),
            "--outdir", str(output),
            "--index", str(args.cagra_index),
            "--reference", str(args.reference_labels),
            "--dtype", "int8",
            "--int8-scale", str(args.cagra_int8_scale),
            "--k", str(candidate),
            "--output-k", str(args.output_k),
            "--itopk-size", str(max(args.cagra_itopk_size, candidate)),
            "--metric", "inner_product",
            "--search-batch-size", str(args.cagra_search_batch_size),
            "--rerank-batch-size", str(args.rerank_batch_size),
            "--candidate-batch-size", str(args.candidate_batch_size),
            "--rerank-device", "cuda",
            "--recall-k-values", "100,500",
        ]
        run(command, args.outdir / "cagra-int8.log", env)


def run_window_ivf(args: argparse.Namespace, env: dict[str, str]) -> None:
    if not args.window_ivf:
        return
    executable = args.window_ivf_executable or args.outdir.parent / "tooling/window_ivf"
    command = [
        sys.executable,
        str(ROOT / "bin/benchmark_window_ivf.py"),
        "--embeddings", str(args.embeddings),
        "--structures", str(args.structures),
        "--query-selections", str(args.query_selections),
        "--database-windows", str(args.database_windows),
        "--queries", str(args.queries),
        "--reference-labels", str(args.reference_labels),
        "--executable", str(executable),
        "--outdir", str(args.outdir / "window-ivf"),
        "--pq-m-values", args.window_ivf_pq_m_values,
        "--nbits-values", args.window_ivf_nbits_values,
        "--nlist-values", args.window_ivf_nlist_values,
        "--nprobe-values", args.window_ivf_nprobe_values,
        "--candidate-k-values", args.candidate_k_values,
        "--output-k", str(args.output_k),
        "--recall-k-values", "100,500",
        "--sample-size", str(args.sample_size),
        "--pq-niter", str(args.pq_niter),
        "--batch-size", str(args.encode_batch_size),
        "--rerank-batch-size", str(args.rerank_batch_size),
        "--candidate-batch-size", str(args.candidate_batch_size),
        "--rerank-workers", str(args.rerank_workers),
        "--rerank-device", "cpu",
        "--threads", str(args.hnsw_threads),
        "--window-size", str(args.window_size),
        "--stride", str(args.stride),
        "--original-dtype", "float16",
    ]
    result = args.outdir / "window-ivf" / "results.json"
    if not result.exists() or args.force:
        run(command, args.outdir / "window-ivf.log", env)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, required=True)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--query-selections", type=Path, required=True)
    parser.add_argument("--database-windows", type=Path, required=True)
    parser.add_argument("--queries", type=Path, required=True)
    parser.add_argument("--reference-labels", type=Path, required=True)
    parser.add_argument("--pq-executable", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--results-name", default="results.json")
    parser.add_argument("--pq-m-values", default="8,16")
    parser.add_argument("--nbits-values", default="4,8")
    parser.add_argument("--hnsw-m-values", default="32,64")
    parser.add_argument("--ef-construction-values", default="200,400")
    parser.add_argument("--ef-search-values", default="1000,5000,10000")
    parser.add_argument("--candidate-k-values", default="1000,5000,10000")
    parser.add_argument("--output-k", type=int, default=500)
    parser.add_argument("--sample-size", type=int, default=500_000)
    parser.add_argument("--pq-niter", type=int, default=25)
    parser.add_argument("--encode-batch-size", type=int, default=32768)
    parser.add_argument("--rerank-batch-size", type=int, default=8)
    parser.add_argument("--candidate-batch-size", type=int, default=1024)
    parser.add_argument("--rerank-workers", type=int, default=8)
    parser.add_argument("--hnsw-threads", type=int, default=16)
    parser.add_argument("--window-size", type=int, default=11)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--cagra", action="store_true")
    parser.add_argument("--cagra-index", type=Path, required=False, default=Path("cagra.index"))
    parser.add_argument("--cagra-int8-scale", type=float, default=850.0)
    parser.add_argument("--cagra-itopk-size", type=int, default=256)
    parser.add_argument("--cagra-search-batch-size", type=int, default=8)
    parser.add_argument("--window-ivf", action="store_true")
    parser.add_argument("--window-ivf-executable", type=Path)
    parser.add_argument("--window-ivf-pq-m-values", default="16")
    parser.add_argument("--window-ivf-nbits-values", default="4")
    parser.add_argument("--window-ivf-nlist-values", default="16,32")
    parser.add_argument("--window-ivf-nprobe-values", default="1,4,16,32")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.cagra_index == Path("cagra.index"):
        args.cagra_index = args.outdir / "cagra-int8.index"
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    env.setdefault("MKL_NUM_THREADS", "1")
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    try:
        run_pq(args, env)
        run_cagra(args, env)
        run_window_ivf(args, env)
    except (OSError, RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    (args.outdir / "run-complete.json").write_text(
        json.dumps({"completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}, indent=2) + "\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
