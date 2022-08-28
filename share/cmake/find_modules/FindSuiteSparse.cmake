#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# SuiteSparse Find Module ------------------------------------------------------
# -------------------------------------------------------------------------

set(SUITESPARSE_ENV_VARS
    $ENV{SUITESPARSEDIR}
    $ENV{SUITESPARSE_DIR}
    $ENV{SuiteSparse_DIR}
    $ENV{SUITESPARSE_ROOT}
    $ENV{SuiteSparse_ROOT}
    $ENV{SUITESPARSE_PATH}
    $ENV{SuiteSparse_PATH} )

find_path(SUITESPARSE_DIR
    NAMES
    include/umfpack.h
    HINTS
    ${SUITESPARSE_ENV_VARS}
    )

find_library(SUITESPARSE_LIBRARIES
    NAMES
    umfpack
    HINTS
    ${SUITESPARSE_DIR}
    PATH_SUFFIXES
    lib
    lib64
    )

find_path(SUITESPARSE_INCLUDE_DIRS
    NAMES
    umfpack.h
    HINTS
    ${SUITESPARSE_DIR}/include
    )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUITESPARSE DEFAULT_MSG SUITESPARSE_LIBRARIES SUITESPARSE_INCLUDE_DIRS)

mark_as_advanced(SUITESPARSE_DIR
    SUITESPARSE_LIBRARIES
    SUITESPARSE_INCLUDE_DIRS
    )

