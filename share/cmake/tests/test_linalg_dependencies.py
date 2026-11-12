#!/usr/bin/env python3
# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.
"""Check the real LINALG dependency selection without installed third parties."""
from pathlib import Path
import subprocess
import tempfile

cmake_dir = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix="moris-linalg-dependencies-") as temporary:
    root = Path(temporary)
    wrappers = root / "wrappers"
    wrappers.mkdir()
    dependencies = root / "dependencies"
    dependencies.mkdir()
    for name in ("eigen", "armadillo", "openblas", "superlu", "arpack", "gperftools"):
        (wrappers / (name + "_new.cmake")).write_text(
            'set_property(GLOBAL APPEND PROPERTY REQUESTED_TPLS "' + name + '")\n')
    (wrappers / "suitesparse_new.cmake").write_text(
        'message(FATAL_ERROR "Dense LINALG requested optional SuiteSparse")\n')
    for name in ("COM", "ALG"):
        (dependencies / (name + "_Depends.cmake")).write_text("# Dependency stand-in.\n")
    (root / "CMakeLists.txt").write_text('''
cmake_minimum_required(VERSION 3.17)
project(LinalgDependencies LANGUAGES NONE)
set(LINALG LINALG)
set(COM COM)
set(ALG ALG)
set(MORIS_TPL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/wrappers")
set(MORIS_DEPENDS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies")
set(ARMADILLO_EIGEN "${BACKEND}")
set(ACML_LAPACK_MKL_OPENBLAS openblas)
if(BACKEND STREQUAL "eigen")
    set(MORIS_USE_EIGEN ON)
    set(MORIS_USE_ARMA OFF)
    set(expected_implementation Eigen_Impl)
else()
    set(MORIS_USE_EIGEN OFF)
    set(MORIS_USE_ARMA ON)
    set(expected_implementation Arma_Impl)
endif()
include("${MORIS_CMAKE_DIR}/dependencies/LINALG_Depends.cmake")
set(expected "${BACKEND};openblas;superlu;arpack")
if(MORIS_USE_GPERFTOOLS)
    list(APPEND expected gperftools)
endif()
get_property(requested GLOBAL PROPERTY REQUESTED_TPLS)
if(NOT "${LINALG_TPL_DEPENDENCIES}" STREQUAL "${expected}" OR
   NOT "${requested}" STREQUAL "${expected}")
    message(FATAL_ERROR "Unexpected LINALG dependencies: ${requested}")
endif()
if(NOT LINALG_IMPLEMENTATION_INCLUDES STREQUAL expected_implementation)
    message(FATAL_ERROR "Incorrect backend include directory")
endif()
# Including dependency metadata twice must not rediscover third parties.
include("${MORIS_CMAKE_DIR}/dependencies/LINALG_Depends.cmake")
get_property(repeated GLOBAL PROPERTY REQUESTED_TPLS)
if(NOT "${repeated}" STREQUAL "${requested}")
    message(FATAL_ERROR "Repeated dependency discovery")
endif()
''')
    for backend, profiling in (("eigen", "OFF"), ("armadillo", "OFF"),
                                ("eigen", "ON"), ("eigen", "OFF")):
        result = subprocess.run(
            ["cmake", "-S", str(root), "-B", str(root / "build"),
             "-DMORIS_CMAKE_DIR=" + str(cmake_dir), "-DBACKEND=" + backend,
             "-DMORIS_USE_GPERFTOOLS=" + profiling],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        assert result.returncode == 0, result.stdout
print("Eigen/Armadillo selection, reconfiguration, and profiling work without SuiteSparse.")
