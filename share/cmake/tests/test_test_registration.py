#!/usr/bin/env python3
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
"""Check serial fallback, MPI counts, commands, and CTest resource accounting."""
import json
from pathlib import Path
import subprocess
import tempfile

cmake_dir = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="moris-registration-") as temporary:
    root = Path(temporary)
    (root / "main.cpp").write_text("int main() { return 0; }\n")
    (root / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.17)
project(Registration LANGUAGES CXX)
enable_testing()
include("${MORIS_CMAKE_DIR}/utilities/moris_testing.cmake")
set(MORIS_EXECUTE_COMMAND "${CMAKE_COMMAND}")
add_executable(COM-test main.cpp)
moris_register_unit_tests(COM-test 1 2 4)
add_executable(HMR-test main.cpp)
moris_register_unit_tests(HMR-test 1 2 4 12)
''')
    for parallel in ("ON", "OFF", "ON"):
        build = root / "build"
        subprocess.run(["cmake", "-S", str(root), "-B", str(build),
                        "-DMORIS_CMAKE_DIR=" + str(cmake_dir),
                        "-DMORIS_HAVE_PARALLEL_TESTS=" + parallel], check=True,
                       stdout=subprocess.PIPE)
        # CTest resolves target paths only after executables have been built.
        subprocess.run(["cmake", "--build", str(build)], check=True, stdout=subprocess.PIPE)
        result = subprocess.check_output(["ctest", "--test-dir", str(build), "--show-only=json-v1"])
        tests = json.loads(result)["tests"]
        expected = {module + "-test-" + str(n) + "-procs": n
                    for module, counts in (("COM", (1, 2, 4)), ("HMR", (1, 2, 4, 12)))
                    for n in (counts if parallel == "ON" else (1,))}
        assert {t["name"] for t in tests} == set(expected), tests
        for test in tests:
            count = expected[test["name"]]
            assert Path(test["command"][0]).name == "cmake", test
            assert test["command"][1:3] == ["-n", str(count)], test
            properties = {p["name"]: p["value"] for p in test["properties"]}
            assert properties.get("PROCESSORS", 1) == count, properties
print("Serial fallback, MPI registration, cache transitions, and processor counts passed.")
