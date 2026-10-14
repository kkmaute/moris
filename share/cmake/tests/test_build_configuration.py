#!/usr/bin/env python3
"""Check compiler policy and plugin helpers without MORIS third-party libraries.

Run: python3 share/cmake/tests/test_build_configuration.py [path/to/cmake]
"""
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import textwrap


def run(*args):
    subprocess.run(args, check=True)


def main():
    cmake = sys.argv[1] if len(sys.argv) > 1 else "cmake"
    cmake_dir = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="moris-build-policy-") as temporary:
        root = Path(temporary)
        source = root / "source"
        (source / "headers").mkdir(parents=True)
        (source / "extra_headers").mkdir()
        (source / "mpi_headers").mkdir()
        (source / "CMakeLists.txt").write_text(textwrap.dedent('''
            cmake_minimum_required(VERSION 3.17...4.4)
            project(BuildPolicy LANGUAGES C CXX)
            include("${MORIS_CMAKE_DIR}/moris_options.cmake")
            include("${MORIS_CMAKE_DIR}/moris_compiler.cmake")
            include("${MORIS_CMAKE_DIR}/utilities/dynamic_link_input.cmake")
            if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE STREQUAL EXPECTED_CONFIG)
                message(FATAL_ERROR "Unexpected build configuration: ${CMAKE_BUILD_TYPE}")
            endif()
            if(NOT CMAKE_CXX_FLAGS STREQUAL "-DUSER_CMAKE_FLAG=7")
                message(FATAL_ERROR "The toolchain/user flags were overwritten")
            endif()
            add_library(payload STATIC payload.cpp)
            set(MORIS moris)
            add_library(moris::mpi INTERFACE IMPORTED)
            set_target_properties(moris::mpi PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/mpi_headers")
            # A static archive must remain usable from a shared plugin.
            set(plugin_includes "${CMAKE_CURRENT_SOURCE_DIR}/headers" "${CMAKE_CURRENT_SOURCE_DIR}/extra_headers")
            # Production callers expand lists; verify directories beyond the first.
            dynamic_link_input(input_plugin policy_plugin plugin.cpp ${plugin_includes})
            target_link_libraries(input_plugin PRIVATE payload)
            add_executable(check main.cpp)
            target_link_libraries(check PRIVATE input_plugin)
        '''))
        (source / "headers/value.hpp").write_text("int value();\nint plugin_value();\n")
        (source / "extra_headers/extra.hpp").write_text('#include "value.hpp"\n')
        (source / "mpi_headers/mpi.h").write_text('// Stand-in for MPI usage requirements.\n')
        (source / "payload.cpp").write_text(textwrap.dedent('''
            #include <optional>
            #if USER_CMAKE_FLAG != 7 || USER_MORIS_FLAG != 11
            #error Missing user-supplied compiler flags
            #endif
            int payload = 40;
            int value() { return std::optional<int>(payload).value() + 2; }
        '''))
        (source / "plugin.cpp").write_text('#include "extra.hpp"\n#include <mpi.h>\nint plugin_value() { return value(); }\n')
        (source / "main.cpp").write_text('int plugin_value();\nint main() { return plugin_value() == 42 ? 0 : 1; }\n')
        cases = [
            ("debug", "Debug", ["-DCMAKE_BUILD_TYPE=Debug", "-DMORIS_USE_LTO=OFF"]),
            ("symbolic", "Debug", ["-DMORIS_HAVE_SYMBOLIC=ON", "-DMORIS_HAVE_SYMBOLIC_STRONG=ON"]),
            ("release-ipo", "Release", ["-DCMAKE_BUILD_TYPE=Release", "-DMORIS_USE_LTO=ON"]),
        ]
        if shutil.which("ninja"):
            cases.append(("multi-config", "Release", ["-G", "Ninja Multi-Config", "-DMORIS_USE_LTO=OFF"]))
        for name, config, options in cases:
            build = root / name
            run(cmake, "-S", str(source), "-B", str(build),
                "-DMORIS_CMAKE_DIR=" + str(cmake_dir),
                "-DCMAKE_CXX_FLAGS=-DUSER_CMAKE_FLAG=7",
                "-DMORIS_CXX_FLAGS=-DUSER_MORIS_FLAG=11",
                "-DEXPECTED_CONFIG=" + config, *options)
            run(cmake, "--build", str(build), "--config", config, "--parallel", "2")
            executable = build / config / "check" if name == "multi-config" else build / "check"
            run(str(executable))
    print("Compiler flags, C++17, PIC, IPO, configurations, and plugin include checks passed.")


if __name__ == "__main__":
    main()
