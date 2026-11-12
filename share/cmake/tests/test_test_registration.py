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

# Grouped launches must preserve quoting, isolation, and failure propagation.
with tempfile.TemporaryDirectory(prefix="moris-test-groups-") as temporary:
    root = Path(temporary)
    (root / "worker.cmake").write_text('''
if(NOT PAYLOAD STREQUAL "argument with spaces and ]=] [=[ brackets")
    message(FATAL_ERROR "Argument quoting was lost: ${PAYLOAD}")
endif()
file(WRITE marker "ran")
if(FAIL)
    message(FATAL_ERROR "Intentional child failure")
endif()
''')
    (root / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.17)
project(TestGroups LANGUAGES NONE)
enable_testing()
include("${MORIS_CMAKE_DIR}/utilities/moris_testing.cmake")
foreach(group passing failing)
    set(cases "")
    foreach(case first second)
        set(directory "${CMAKE_CURRENT_BINARY_DIR}/${group} ${case}")
        file(MAKE_DIRECTORY "${directory}")
        set(fail OFF)
        if(group STREQUAL "failing" AND case STREQUAL "first")
            set(fail ON)
        endif()
        moris_append_test_group_case(cases "${group}-${case}" "${directory}" 3
            "${CMAKE_COMMAND}" "-DFAIL=${fail}" "-DPAYLOAD=argument with spaces and ]=] [=[ brackets"
            -P "${CMAKE_CURRENT_SOURCE_DIR}/worker.cmake"
        )
    endforeach()
    moris_register_test_group(${group} "${cases}" PROCESSORS 2 TIMEOUT 20 LABELS grouped)
endforeach()
''')
    build = root / "build"
    subprocess.run(["cmake", "-S", str(root), "-B", str(build),
                    "-DMORIS_CMAKE_DIR=" + str(cmake_dir)], check=True, stdout=subprocess.PIPE)
    passing = subprocess.run(["ctest", "--test-dir", str(build), "-R", "^passing$",
                              "--output-on-failure"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert passing.returncode == 0, passing.stdout.decode()
    failure = subprocess.run(["ctest", "--test-dir", str(build), "-R", "^failing$",
                              "--output-on-failure"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    assert failure.returncode != 0, failure.stdout.decode()
    assert "failing-first" in failure.stdout.decode(), failure.stdout.decode()
    for group in ("passing", "failing"):
        for case in ("first", "second"):
            assert (build / (group + " " + case) / "marker").read_text() == "ran"
print("Grouped cases preserve arguments and working directories, continue after failure, and fail CTest.")

# Configure the actual EXA registrations with lightweight TPL/target stubs.
with tempfile.TemporaryDirectory(prefix="moris-exa-registration-") as temporary:
    root = Path(temporary)
    (root / "tpl").mkdir()
    (root / "tpl/trilinos_new.cmake").write_text("# No TPL discovery needed for registration.\n")
    (root / "petsc").mkdir()
    (root / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.17)
project(ExampleRegistration LANGUAGES CXX)
enable_testing()
include("${MORIS_SOURCE_DIR}/share/cmake/utilities/moris_testing.cmake")
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
add_subdirectory("${MORIS_SOURCE_DIR}/projects/EXA/structure/linear/bddc" bddc)
''')
    repository = cmake_dir.parents[1]
    for petsc, parallel, mumps in ((True, True, True), (True, False, True),
                                   (True, True, False), (False, True, False),
                                   (False, False, False)):
        (root / "petsc/petscconf.h").write_text("#define PETSC_HAVE_MUMPS 1\n" if mumps else "// No MUMPS\n")
        build = root / ("build-" + str(int(petsc)) + str(int(parallel)) + str(int(mumps)))
        subprocess.run(["cmake", "-S", str(root), "-B", str(build),
                        "-DMORIS_SOURCE_DIR=" + str(repository),
                        "-DMORIS_HAVE_PETSC=" + ("ON" if petsc else "OFF"),
                        "-DMORIS_HAVE_PARALLEL_TESTS=" + ("ON" if parallel else "OFF")],
                       check=True, stdout=subprocess.PIPE)
        result = subprocess.check_output(["ctest", "--test-dir", str(build), "--show-only=json-v1"])
        tests = json.loads(result)["tests"]
        expected = {"Thick_Pressure_Vessel-" + str(n) + "-procs": n
                    for n in ((1, 4) if parallel else (1,))}
        if petsc:
            expected["EXA-structure-linear-bddc"] = 8 if parallel and mumps else 1
        assert {test["name"] for test in tests} == set(expected), tests
        directories = set()
        for test in tests:
            properties = {p["name"]: p["value"] for p in test["properties"]}
            assert properties.get("PROCESSORS", 1) == expected[test["name"]], properties
            script = Path(test["command"][2]).read_text()
            cases = [line for line in script.splitlines() if line.startswith("moris_run_test_case(")]
            if test["name"] == "EXA-structure-linear-bddc":
                expected_cases = 12 if parallel and mumps else 4
            else:
                count = expected[test["name"]]
                bddc = petsc and (count == 1 or mumps)
                expected_cases = 3 if bddc else 1
                assert ("Pressure_Vessel_3D_BDDC" in script) == bddc, script
                assert ("Pressure_Vessel_3D_Immersed_BDDC" in script) == bddc, script
            assert len(cases) == expected_cases, script
            assert properties["TIMEOUT"] == 600 * expected_cases + 30, properties
            for case in cases:
                directory = Path(case.split("[=[")[2].split("]=]")[0])
                assert directory.is_dir(), directory
                assert directory not in directories, directory
                directories.add(directory)
print("EXA groups retain all cases, isolate outputs, reserve MPI slots, and honor PETSc/MUMPS/serial gating.")
