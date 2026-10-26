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

# Configure the actual vessel registration with lightweight TPL/target stubs.
# These cases share VIS temporary filenames and must keep all output isolated,
# including the legacy Amesos tests that can run beside the new BDDC cases.
with tempfile.TemporaryDirectory(prefix="moris-vessel-registration-") as temporary:
    root = Path(temporary)
    (root / "tpl").mkdir()
    (root / "tpl/trilinos_new.cmake").write_text("# No TPL discovery needed for registration.\n")
    (root / "petsc").mkdir()
    (root / "petsc/petscconf.h").write_text("#define PETSC_HAVE_MUMPS 1\n")
    (root / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.17)
project(VesselRegistration LANGUAGES CXX)
enable_testing()
set(MORIS_HAVE_PETSC ON)
set(MORIS_HAVE_PARALLEL_TESTS ON)
set(MORIS_EXECUTE_COMMAND "${CMAKE_COMMAND}")
set(MORIS_PACKAGE_DIR "${MORIS_SOURCE_DIR}/projects")
set(MORIS_TPL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tpl")
set(PETSC_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/petsc")
set(BIN bin)
set(EXE_EXT .exe)
function(dynamic_link_input target output source)
    add_library(${target} SHARED "${source}")
endfunction()
add_subdirectory("${MORIS_SOURCE_DIR}/projects/EXA/structure/linear/Thick_Pressure_Vessel" vessel)
''')
    build = root / "build"
    repository = cmake_dir.parents[1]
    subprocess.run(["cmake", "-S", str(root), "-B", str(build),
                    "-DMORIS_SOURCE_DIR=" + str(repository)], check=True,
                   stdout=subprocess.PIPE)
    # CTest resolves executable paths while listing tests. Supply an inert
    # executable; no MORIS source or third-party library is compiled here.
    executable = build / "vessel/bin/Main_Thick_Pressure_Vessel.exe"
    executable.parent.mkdir(parents=True, exist_ok=True)
    executable.write_text("#!/bin/sh\nexit 0\n")
    executable.chmod(0o755)
    result = subprocess.check_output(["ctest", "--test-dir", str(build), "--show-only=json-v1"])
    tests = json.loads(result)["tests"]
    expected = {"Thick_Pressure_Vessel-" + str(n) + "-procs": n for n in (1, 4)}
    expected.update({example + "-BDDC-" + str(n) + "-procs": n
                     for example in ("Pressure_Vessel_3D", "Pressure_Vessel_3D_Immersed")
                     for n in (1, 4)})
    assert {test["name"] for test in tests} == set(expected), tests
    directories = set()
    for test in tests:
        properties = {p["name"]: p["value"] for p in test["properties"]}
        directory = Path(properties["WORKING_DIRECTORY"])
        assert directory == build / "vessel/cases" / test["name"], properties
        assert directory.is_dir(), directory
        assert directory not in directories, directory
        directories.add(directory)
        assert properties.get("PROCESSORS", 1) == expected[test["name"]], properties
print("Pressure-vessel output directories are isolated for all Amesos/BDDC process counts.")
