#!/usr/bin/env python3
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
"""Require an explicit assertion-test mapping for every public LINALG header.

This checks the inventory and CMake registration, not numerical correctness or
instrumented line coverage. Run LINALG-test under both backends for those checks.
"""
import json
from pathlib import Path
import re

repository = Path(__file__).resolve().parents[3]
source = repository / "projects/LINALG/src"
tests = repository / "projects/LINALG/test"
inventory = json.loads((tests / "coverage.json").read_text())
public_headers = {path.name for path in source.glob("*.hpp")}
assert set(inventory) == public_headers, {
    "unmapped": sorted(public_headers - set(inventory)),
    "stale": sorted(set(inventory) - public_headers),
}
block = re.search(r"set\(TEST_SOURCES\s+(.*?)\)",
                  (tests / "CMakeLists.txt").read_text(), re.S).group(1)
registered = set(re.findall(r"\w+\.cpp", block))
assertions = re.compile(r"\b(?:REQUIRE|CHECK|STATIC_REQUIRE)(?:_\w+)?\s*\(")
for header, sources in inventory.items():
    assert sources, header
    for filename in sources:
        assert filename in registered, (header, filename, "not registered")
        contents = (tests / filename).read_text()
        assert "TEST_CASE(" in contents or "TEST_CASE (" in contents, (header, filename)
        assert assertions.search(contents), (header, filename, "no assertions")
        assert filename != "LINALG_Tutorial.cpp", (header, "tutorial output is not a test")
print("All {} public LINALG headers map to registered assertion tests.".format(len(public_headers)))
