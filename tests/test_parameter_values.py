#!/usr/bin/env python3
"""Ensure every string-valued schema choice is lowercase."""
from __future__ import annotations

import json
import unittest
from pathlib import Path
from typing import Any, Iterator


SCHEMA = Path(__file__).resolve().parents[1] / "nextflow_schema.json"


def schema_nodes(value: Any) -> Iterator[dict[str, Any]]:
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from schema_nodes(child)
    elif isinstance(value, list):
        for child in value:
            yield from schema_nodes(child)


class TestParameterValues(unittest.TestCase):
    def test_string_enum_values_are_lowercase(self) -> None:
        schema = json.loads(SCHEMA.read_text())
        for node in schema_nodes(schema):
            for value in node.get("enum", []):
                if isinstance(value, str):
                    with self.subTest(value=value):
                        self.assertEqual(value, value.lower())


if __name__ == "__main__":
    raise SystemExit(unittest.main())
