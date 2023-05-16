#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# LAPACK Find Module ------------------------------------------------------
# -------------------------------------------------------------------------

set(LAPACK_ENV_VARS
    $ENV{LAPACKDIR}
    $ENV{LAPACK_DIR}
    $ENV{lapack_DIR}
    $ENV{LAPACK_ROOT}
    $ENV{lapack_ROOT}
    $ENV{LAPACK_PATH}
    $ENV{lapack_PATH} )

find_library(LAPACK_lapack
    NAMES
    lapack
    HINTS
    ${LAPACK_ENV_VARS}
    PATH_SUFFIXES
    lib64 )

find_library(LAPACK_blas
    NAMES
    blas
    HINTS
    ${LAPACK_ENV_VARS}
    PATH_SUFFIXES
    lib64 )

set(LAPACK_LIBRARIES
    ${LAPACK_lapack}
    ${LAPACK_blas}
    CACHE FILEPATH "List of library paths." )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LAPACK DEFAULT_MSG LAPACK_LIBRARIES)

mark_as_advanced(LAPACK_lapack LAPACK_blas LAPACK_LIBRARIES)

_import_libraries(LAPACK_LIBRARY_TARGETS ${LAPACK_LIBRARIES})

add_library(LAPACK::lapack INTERFACE IMPORTED GLOBAL)
target_link_libraries(LAPACK::lapack INTERFACE ${LAPACK_LIBRARY_TARGETS})

add_library(LAPACK::all_libs INTERFACE IMPORTED)
target_link_libraries(LAPACK::all_libs
  INTERFACE LAPACK::lapack
  )
