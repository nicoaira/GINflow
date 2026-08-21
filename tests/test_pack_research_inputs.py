import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parents[1] / "bin"))
from pack_research_inputs import build_reference, pack_nodes  # noqa: E402


class TestPackResearchInputs(unittest.TestCase):
    def test_blockwise_reference_matches_flat_inner_product(self):
        rng = np.random.default_rng(17)
        database = rng.normal(size=(73, 12)).astype(np.float32)
        queries = rng.normal(size=(5, 12)).astype(np.float32)
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "labels.npy"
            metadata = build_reference(database, queries, output, 9, 11)
            expected = np.argsort(-(queries @ database.T), axis=1)[:, :9]
            np.testing.assert_array_equal(np.load(output), expected)
            self.assertEqual(metadata["index_type"], "blockwise-IndexFlatIP")

    def test_nodes_follow_requested_identifier_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            database = root / "database"
            database.mkdir()
            np.savez(
                database / "embeddings.npz",
                first=np.asarray([[1, 2], [3, 4]], dtype=np.float32),
                second=np.asarray([[5, 6]], dtype=np.float32),
            )
            output = root / "nodes.npy"
            values = pack_nodes(database, ["second", "first"], output, np.dtype("float16"))
            np.testing.assert_array_equal(
                np.asarray(values),
                np.asarray([[5, 6], [1, 2], [3, 4]], dtype=np.float16),
            )


if __name__ == "__main__":
    unittest.main()
