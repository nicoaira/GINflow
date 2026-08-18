from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from live_dashboard import render


def _result(*, backend: str, label: str, device: str, qps: float, recall: float) -> dict:
    return {
        "status": "ok",
        "backend": backend,
        "parameter_label": label,
        "run_id": f"{backend}-run",
        "repeat": 0,
        "dataset_window_count": 10,
        "dimension": 4,
        "metric": "cosine",
        "k": 100,
        "parameters": {"device": device, "index_type": label},
        "recall_at_100": recall,
        "qps": qps,
        "latency_ms": {"mean": 1.0, "p50": 1.0, "p95": 1.2},
        "index_bytes": 160,
        "build_seconds": 2.0,
        "peak_rss_bytes": 1024,
        "peak_vram_bytes": 2048 if device == "gpu" else None,
        "warmup_queries": 4,
        "timed_queries": 8,
        "measurement": {
            "build_peak_rss_bytes": 2048,
            "build_peak_vram_bytes": 4096 if device == "gpu" else None,
            "protocol": "fixed-query-batches-v1",
        },
        "provenance": {
            "hardware_id": "fixture",
            "hardware": {"gpus": [{"name": "fixture-gpu"}]} if device == "gpu" else {"gpus": []},
        },
        "query_ids_sha256": "q",
        "ground_truth_ids_sha256": "g",
    }


class TestLiveDashboard(unittest.TestCase):
    def test_interactive_points_include_metadata_and_gpu_classes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result_dir = root / "dataset" / "results"
            for backend, label, device, qps, recall in (
                ("faiss", "cpu-flat", "cpu", 100.0, 0.9),
                ("cuvs", "gpu-cagra", "gpu", 200.0, 0.95),
            ):
                path = result_dir / backend / f"{label}.json"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(json.dumps(_result(backend=backend, label=label, device=device, qps=qps, recall=recall)), encoding="utf-8")
            output = root / "dashboard" / "index.html"
            render(root, output, refresh_seconds=15)
            html = output.read_text(encoding="utf-8")
            self.assertIn("Running benchmarks", html)
            self.assertIn('id="gpu-highlight"', html)
            self.assertIn('id="pareto-toggle"', html)
            self.assertIn('id="chart-tooltip"', html)
            self.assertIn('class="chart-point gpu-point"', html)
            self.assertIn('class="chart-point cpu-point"', html)
            self.assertIn("pointerenter", html)
            self.assertIn("paretoFrontier", html)
            self.assertIn("parameters", html)
            self.assertIn("provenance", html)


if __name__ == "__main__":
    unittest.main()
