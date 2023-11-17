#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# MKL Find Module ---------------------------------------------------------
# -------------------------------------------------------------------------

set(MKL_ENV_VARS
    $ENV{MKLDIR}
    $ENV{MKL_DIR}
    $ENV{mkl_DIR}
    $ENV{MKL_ROOT}
    $ENV{mkl_ROOT}
    $ENV{MKL_PATH}
    $ENV{mkl_PATH} )

find_path(MKL_DIR
    NAMES
    include/mkl.h
    HINTS
    ${MKL_ENV_VARS}
    PATHS
    /usr/lib/mkl
    /usr/lib
    /usr
    /usr/local )

find_path(MKL_INCLUDE_DIRS
    NAMES
    mkl.h
    HINTS
    ${MKL_DIR}
    ${MKL_ENV_VARS}
    PATH_SUFFIXES
    include )

find_library(MKL_lp64
    NAMES
    mkl_intel_lp64
    HINTS
    ${MKL_DIR}/lib/intel64
    ${MKL_DIR}/lib/intel64_lin )

find_library(MKL_core
    NAMES
    mkl_core
    HINTS
    ${MKL_DIR}/lib/intel64
    ${MKL_DIR}/lib/intel64_lin )

find_library(MKL_sequential
    NAMES
    mkl_sequential
    HINTS
    ${MKL_DIR}/lib/intel64
    ${MKL_DIR}/lib/intel64_lin )

find_library(MKL_pthread
    NAMES
    pthread
    HINTS
    ${MKL_DIR}/lib/intel64
    ${MKL_DIR}/lib/intel64_lin )

set(MKL_LIBRARIES
    ${MKL_lp64}
    ${MKL_sequential}
    ${MKL_core}
#    ${MKL_pthread}
    CACHE FILEPATH "List of library paths." )


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG
    MKL_DIR
    MKL_INCLUDE_DIRS
    MKL_lp64 MKL_core
    MKL_sequential
    MKL_pthread
    MKL_LIBRARIES )

mark_as_advanced(MKL_DIR
    MKL_INCLUDE_DIRS
    MKL_lp64 MKL_core
    MKL_sequential
    MKL_pthread
    MKL_LIBRARIES )

_import_libraries(MKL_LIBRARY_TARGETS "${MKL_LIBRARIES}")


if(NOT TARGET MKL::mkl)
    add_library(MKL::mkl INTERFACE IMPORTED GLOBAL)
    target_link_libraries(MKL::mkl INTERFACE ${MKL_LIBRARY_TARGETS})
endif()

if(NOT TARGET MKL::all_libs)
    add_library(MKL::all_libs INTERFACE IMPORTED)
    target_link_libraries(MKL::all_libs INTERFACE ${MKL_LIBRARY_TARGETS})
endif()
