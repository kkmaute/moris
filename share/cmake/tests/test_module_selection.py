#!/usr/bin/env python3
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
"""Check ALL/NONE and legacy selectors without compilers or third-party libraries."""
from pathlib import Path
import subprocess
import tempfile


def main():
    cmake_dir = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory(prefix="moris-module-selection-") as temporary:
        root = Path(temporary)
        (root / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.17)
project(Selection NONE)
include("${MORIS_CMAKE_DIR}/moris_options.cmake")
include("${MORIS_CMAKE_DIR}/utilities/moris_select_modules.cmake")
moris_select_modules()
if(DEFINED EXPECT_PROFILING AND NOT "${EXPECT_PROFILING}" STREQUAL "${MORIS_USE_GPERFTOOLS}")
    message(FATAL_ERROR "Unexpected profiling option precedence")
endif()
if(DEFINED EXPECT_SDF_EXE AND NOT "${EXPECT_SDF_EXE}" STREQUAL "${BUILD_SDF_EXE}")
    message(FATAL_ERROR "Unexpected SDF executable selection")
endif()
if(EXPECT_ALL)
    foreach(component IN LISTS MORIS_BUILD_MODULES)
        if(NOT component IN_LIST MORIS_ALL_EXCLUSIONS AND NOT BUILD_${component})
            message(FATAL_ERROR "ALL omitted ${component}")
        endif()
    endforeach()
    if(BUILD_STK OR BUILD_MOD OR NOT BUILD_MAIN OR NOT MORIS_USE_EXAMPLES)
        message(FATAL_ERROR "Incorrect ALL selection")
    endif()
endif()
if(EXPECT_NONE)
    foreach(component IN LISTS MORIS_BUILD_MODULES)
        if(BUILD_${component})
            message(FATAL_ERROR "NONE left ${component} enabled")
        endif()
    endforeach()
    if(BUILD_MAIN OR BUILD_SDF_EXE OR BUILD_TUTORIALS OR MORIS_USE_EXAMPLES)
        message(FATAL_ERROR "NONE left executables enabled")
    endif()
endif()
if(BUILD_ALL OR BUILD_NONE)
    message(FATAL_ERROR "Selection controls did not reset")
endif()
if(NOT "${BUILD_EXA}" STREQUAL "${MORIS_USE_EXAMPLES}")
    message(FATAL_ERROR "Example selectors disagree")
endif()
if(DEFINED EXPECT_EXAMPLES AND NOT "${EXPECT_EXAMPLES}" STREQUAL "${MORIS_USE_EXAMPLES}")
    message(FATAL_ERROR "Unexpected examples selection")
endif()
if(DEFINED EXPECT_GUI AND NOT "${EXPECT_GUI}" STREQUAL "${BUILD_GUI}")
    message(FATAL_ERROR "Unexpected GUI selection")
endif()
''')

        def configure(name, *options, failure=None):
            result = subprocess.run(
                ["cmake", "-S", str(root), "-B", str(root / name),
                 "-DMORIS_CMAKE_DIR=" + str(cmake_dir), *options],
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
            if failure:
                assert result.returncode and failure in result.stdout, result.stdout
            elif result.returncode:
                raise RuntimeError(result.stdout)

        configure("sdf-default", "-DEXPECT_SDF_EXE=OFF")
        configure("sdf-enabled", "-DBUILD_SDF_EXE=ON", "-DEXPECT_SDF_EXE=ON")
        configure("sdf-enabled", "-DBUILD_ALL=ON", "-DEXPECT_SDF_EXE=ON")
        configure("sdf-enabled", "-DBUILD_NONE=ON", "-DEXPECT_SDF_EXE=OFF")
        configure("all", "-DEXPECT_SDF_EXE=OFF", "-DBUILD_ALL=ON", "-DEXPECT_ALL=ON", "-DEXPECT_GUI=OFF")
        configure("gui", "-DBUILD_ALL=ON", "-DBUILD_GUI=ON", "-DEXPECT_ALL=ON", "-DEXPECT_GUI=ON")
        configure("all", "-DBUILD_NONE=ON", "-DEXPECT_ALL=OFF", "-DEXPECT_NONE=ON")
        configure("all", "-DEXPECT_NONE=ON")  # Reconfiguration remains stable.
        configure("conflict", "-DBUILD_ALL=ON", "-DBUILD_NONE=ON",
                  failure="BUILD_ALL and BUILD_NONE cannot both be ON")
        configure("alias", "-DBUILD_EXA=OFF", "-DEXPECT_EXAMPLES=OFF")
        configure("alias", "-DBUILD_EXA=ON", "-DEXPECT_EXAMPLES=ON")
        configure("alias", "-DMORIS_USE_EXAMPLES=OFF", "-DEXPECT_EXAMPLES=OFF")
        configure("alias", "-DBUILD_ALL=ON", "-DEXPECT_EXAMPLES=ON")
        configure("profiling-alias", "-DUSE_GPERFTOOLS=ON", "-DEXPECT_PROFILING=ON")
        configure("profiling-canonical", "-DUSE_GPERFTOOLS=ON",
                  "-DMORIS_USE_GPERFTOOLS=OFF", "-DEXPECT_PROFILING=OFF")
    print("Module selection, exclusions, cache transitions, and examples alias checks passed.")


if __name__ == "__main__":
    main()
