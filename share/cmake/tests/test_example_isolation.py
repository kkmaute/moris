#!/usr/bin/env python3
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
"""Exercise EXA startup, owned-output cleanup, and configuration-aware runs."""
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

repository = Path(__file__).resolve().parents[3]


def run(*args, **kwargs):
    result = subprocess.run(list(map(str, args)), capture_output=True,
                            text=True, **kwargs)
    if result.returncode:
        print(result.stdout)
        print(result.stderr)
    result.check_returncode()
    return result


with tempfile.TemporaryDirectory(prefix="moris exa isolation ") as temporary:
    root = Path(temporary)
    stubs = root / "stubs"
    stubs.mkdir()
    (stubs / "catch.hpp").write_text(
        "namespace Catch { struct Session { int run(int, char**) { return 0; } }; }\n")
    (stubs / "cl_Communication_Manager.hpp").write_text("""
namespace moris {
struct Comm_Manager {
    Comm_Manager() = default;
    Comm_Manager(int*, char***) {}
    void finalize() {}
};
inline int par_rank() { return 0; }
}
""")
    (stubs / "cl_Logger.hpp").write_text(
        "namespace moris { struct Logger { void initialize(int) {} }; }\n")
    (stubs / "banner.hpp").write_text(
        "namespace moris { inline void print_banner(int, char**) {} }\n")
    (stubs / "Kokkos_Core.hpp").write_text("// ArborX disabled in this fixture.\n")
    compiler = shutil.which("g++-13") or shutil.which("c++")
    startup = root / "startup"
    run(compiler, "-std=c++17", "-I", stubs,
        repository / "projects/EXA/src/example_main.cpp", "-o", startup)

    # Listing tests must preserve inputs and unrelated output-like filenames.
    protected = ("restart_input.hdf5", "notes_exodus.txt", "result.exo")
    for name in protected:
        (root / name).write_text("preserve\n")
    run(startup, "--list-tests", cwd=root)
    assert all((root / name).read_text() == "preserve\n" for name in protected)

    (stubs / "cl_Communication_Tools.hpp").write_text("""
#pragma once
namespace moris {
using sint = int;
inline int rank = 0, size = 1, reductions = 0;
inline int par_rank() { return rank; }
inline int par_size() { return size; }
inline int max_all(int value) { ++reductions; return value; }
}
""")
    (stubs / "fn_assert.hpp").write_text("""
#include <stdexcept>
#define MORIS_ERROR(condition, ...) do { if (!(condition)) throw std::runtime_error("cleanup failed"); } while (false)
""")
    cleanup_source = root / "cleanup.cpp"
    cleanup_source.write_text("""
#include "fn_EXA_Test_Utilities.hpp"
#include <cassert>
int main(int argc, char** argv) {
    moris::size = std::stoi(argv[1]);
    moris::rank = std::stoi(argv[2]);
    moris::exa::remove_test_outputs({"result.exo", "sensitivity.hdf5"});
    assert(moris::reductions == 1);
}
""")
    cleanup = root / "cleanup"
    run(compiler, "-std=c++17", "-I", stubs, "-I",
        repository / "projects/EXA/src", cleanup_source, "-o", cleanup)
    for count, rank in ((1, 0), (2, 0), (2, 1), (12, 3)):
        directory = root / ("cleanup-{}-{}".format(count, rank))
        directory.mkdir()
        exodus = "result.exo" if count == 1 else "result.exo.{}.{}".format(
            count, str(rank).zfill(len(str(count))))
        hdf5 = "sensitivity.hdf5" if count == 1 else "sensitivity_{}.{}.hdf5".format(count, rank)
        for name in (exodus, hdf5, "restart_input.hdf5", "notes_exodus.txt"):
            (directory / name).write_text("fixture\n")
        run(cleanup, count, rank, cwd=directory)
        assert not (directory / exodus).exists()
        assert not (directory / hdf5).exists()
        assert (directory / "restart_input.hdf5").exists()
        assert (directory / "notes_exodus.txt").exists()
        # Missing outputs on a repeat cleanup are harmless.
        run(cleanup, count, rank, cwd=directory)

    fixture = root / "fixture"
    fixture.mkdir()
    (fixture / "input.cpp").write_text("""
#ifndef MORIS_EXAMPLE_INPUT_LIBRARY
#error Missing callback library path
#endif
extern "C" const char* input_path() { return MORIS_EXAMPLE_INPUT_LIBRARY; }
""")
    (fixture / "main.cpp").write_text("""
#include <dlfcn.h>
#include <fstream>
#include <string>
int main(int argc, char** argv) {
    const char* path = std::string(argv[1]) == "first"
        ? MORIS_EXAMPLE_INPUT_first : MORIS_EXAMPLE_INPUT_second;
    void* library = dlopen(path, RTLD_NOW);
    if (!library) return 1;
    auto callback = reinterpret_cast<const char*(*)()>(dlsym(library, "input_path"));
    if (!callback || std::string(callback()) != path) return 2;
    dlclose(library);
    std::ofstream("result.exo") << argv[1];
}
""")
    (fixture / "CMakeLists.txt").write_text("""
cmake_minimum_required(VERSION 3.17)
project(ExampleIsolation LANGUAGES CXX)
enable_testing()
include("${MORIS_SOURCE_DIR}/projects/EXA/example_test_helpers.cmake")
set(MORIS_PACKAGE_DIR "${MORIS_SOURCE_DIR}/projects")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_library(first SHARED input.cpp)
add_library(second SHARED input.cpp)
add_executable(driver main.cpp)
target_link_libraries(driver PRIVATE ${CMAKE_DL_LIBS})
add_test(NAME first-1-procs COMMAND driver first)
add_test(NAME second-2-procs COMMAND driver second)
moris_isolate_example_tests(driver first second)
""")
    for generator in ("Unix Makefiles", "Ninja Multi-Config"):
        build = root / generator
        run("cmake", "-S", fixture, "-B", build, "-G", generator,
            "-DMORIS_SOURCE_DIR=" + str(repository), "-DCMAKE_BUILD_TYPE=Debug")
        configurations = ("Debug", "Release") if "Multi-Config" in generator else ("Debug",)
        for configuration in configurations:
            # Building just the driver must also build both loaded libraries.
            run("cmake", "--build", build, "--config", configuration, "--target", "driver", "-j2")
            tests = json.loads(run("ctest", "--test-dir", build, "-C", configuration,
                                   "--show-only=json-v1").stdout)["tests"]
            directories = set()
            for test in tests:
                properties = {p["name"]: p["value"] for p in test["properties"]}
                directory = Path(properties["WORKING_DIRECTORY"])
                assert directory == build / "runs" / configuration / test["name"]
                assert directory.is_dir()
                assert directory not in directories
                assert properties.get("PROCESSORS", 1) == int(test["name"].split("-")[-2])
                directories.add(directory)
            run("ctest", "--test-dir", build, "-C", configuration, "-j2", "--output-on-failure")
            for test in tests:
                directory = build / "runs" / configuration / test["name"]
                assert (directory / "result.exo").read_text() == test["name"].split("-")[0]

print("EXA startup preservation, exact output cleanup, target paths, and isolated single/multi-config tests passed.")
