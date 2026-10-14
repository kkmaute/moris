# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Public cache options. Keep names and defaults compatible with existing build scripts.
include_guard(GLOBAL)

# MPI support.
option(MORIS_HAVE_PARALLEL "Build MORIS with support for MPI communication." ON)

# MPI package.
set(MORIS_MPI_LIBS "MPICH" "OPENMPI")
set(MORIS_USE_MPI "OPENMPI" CACHE STRING "Set MPI package: ${MORIS_MPI_LIBS}")

# MPI execution command
set(MORIS_EXECUTE_COMMAND "mpirun" CACHE STRING "The command used to execute MORIS tests.")

# SHARED lib support.
option(MORIS_HAVE_SHARED "Build MORIS with support for SHARED libraries." OFF)

option(MORIS_USE_ARMA "Use the Armadillo library as the dense linear algebra package." ON)

option(MORIS_USE_EIGEN "Use the Eigen3 library as the dense linear algebra package." OFF)

option(MORIS_USE_32BIT "Use 32bit numbers for identifiers and sizes." ON)

option(MORIS_USE_LTO "Use LTO optimization." ON)

option(MORIS_HAVE_DEBUG "Use the debug version of the code." OFF)

option(MORIS_USE_VALGRIND "Use valgrind to run tests." OFF)

option(MORIS_USE_MATRIX_FILL "Fill matrices by default with NANs or MaxValues." OFF)

option(MORIS_USE_CHECK_MEMORY "Check memory usage for matrices and cells." OFF)

option(MORIS_USE_XTK "Have XTK Library" ON)

option(MORIS_HAVE_SYMBOLIC "Use symbolic information in executable." OFF)

option(MORIS_HAVE_SYMBOLIC_STRONG "Use symbolic information (strong) in executable." OFF)

option(MORIS_HAVE_PEDANTIC "Use strict compiler checking" ON)

option(MORIS_USE_ACML "Use the AMD ACML library as the linear algebra package." OFF)

option(MORIS_USE_OPENBLAS "Use the open blas library as the linear algebra package." ON)

option(MORIS_USE_LAPACK "Use the lapack library as the linear algebra package." OFF)

option(MORIS_USE_MKL "Use the Inel Math kernel library as the linear algebra package." OFF)

option(MORIS_USE_PARDISO "Build MORIS with PARDISO used in Trilinos and Petsc." ON)

option(MORIS_USE_MUMPS "Build MORIS with MUMPS used in Trilinos and Petsc." ON)

option(MORIS_HAVE_PETSC "Build MORIS with support for PETSC." ON)

option(MORIS_HAVE_SLEPC "Build MORIS with support for SLEPC." ON)

option(MORIS_HAVE_GCMMA "Build MORIS with support for GCMMA." ON)

option(MORIS_HAVE_SNOPT "Build MORIS with support for SNOPT." ON)

option(MORIS_HAVE_LBFGS "Build MORIS with support for LBFGS." ON)

option(MORIS_HAVE_ARBORX "Build MORIS with support for ArborX." ON)

option(MORIS_USE_TESTS "Compile unit tests." ON)

option(MORIS_USE_PERF_TESTS " Build the performance tests" OFF)

# Accept the legacy examples selector on a fresh cache.
if(DEFINED BUILD_EXA AND NOT DEFINED MORIS_USE_EXAMPLES)
    set(MORIS_USE_EXAMPLES "${BUILD_EXA}" CACHE BOOL "Build the examples.")
endif()
option(MORIS_USE_EXAMPLES "Build the examples." ON)

option(MORIS_HAVE_PARALLEL_TESTS "Run unit tests in parallel. MORIS_USE_TESTS must be ON." ON)

option(MORIS_PERFORM_CHECK "Use performance check." ON)

option(MORIS_USE_INTEL "USE BLANCA" OFF)

option(MORIS_USE_OPENMP "Use OpenMP" OFF)

option(MORIS_USE_NTOP_LIB " Build with NTOP interface library - experimental" OFF)

option(MORIS_USE_CLANG_TIDY " Use clang-tidy" OFF)

# Package Options

option(BUILD_MAIN "Build main executable." OFF)

option(BUILD_ALL "Build all executables." OFF)

option(BUILD_NONE "Build no executables." OFF)

include("${CMAKE_CURRENT_LIST_DIR}/moris_modules.cmake")
foreach(component IN LISTS MORIS_BUILD_MODULES)
    option(BUILD_${component} "Build the ${component} executable/tests." OFF)
endforeach()
option(BUILD_TUTORIALS "Build the tutorial executables." OFF)

# Optional instrumentation and experimental integrations.
# Accept the old name only when the canonical option has not been supplied.
if(DEFINED USE_GPERFTOOLS AND NOT DEFINED MORIS_USE_GPERFTOOLS)
    set(MORIS_USE_GPERFTOOLS "${USE_GPERFTOOLS}" CACHE BOOL "Use Google Profiling tools.")
endif()
option(MORIS_USE_GPERFTOOLS "Use Google Profiling tools (check paths in CMakeLists.txt)." OFF)
option(MORIS_USE_MEMORY_DEBUG "Use memkind library for memory debugging)." OFF)
option(MORIS_USE_SACADO "UseSacado (check paths in CMakeLists.txt)." OFF)

# Provide choices to cmake-gui without changing the existing cache values.
set_property(CACHE MORIS_USE_MPI PROPERTY STRINGS ${MORIS_MPI_LIBS})
