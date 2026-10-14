#!/usr/bin/env python3
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
"""Check repository-owned CMake files with the pinned cmake-format version."""
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[3]
paths = [root / "CMakeLists.txt"]
for directory in ("projects", "share"):
    paths.extend(p for p in (root / directory).rglob("*")
                 if p.name == "CMakeLists.txt" or p.name.endswith((".cmake", ".cmake.in")))
subprocess.run(["cmake-format", "--check", *map(str, sorted(paths))], cwd=str(root), check=True)
