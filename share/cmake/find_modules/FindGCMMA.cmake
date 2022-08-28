#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# GCMMA Find Module ------------------------------------------------------
# -------------------------------------------------------------------------

set(GCMMA_ENV_VARS
    $ENV{GCMMADIR}
    $ENV{GCMMA_DIR}
    $ENV{gcmma_DIR}
    $ENV{GCMMA_ROOT}
    $ENV{gcmma_ROOT}
    $ENV{GCMMA_PATH}
    $ENV{gcmma_PATH} )

find_path(GCMMA_DIR
    NAMES
    lib/libgcmma.a
    lib64/libgcmma.a
    HINTS
    ${GCMMA_ENV_VARS} )

find_path(GCMMA_INCLUDE_DIRS
    NAMES
    optalggcmmacall.hpp
    mma.hpp
    HINTS
    ${GCMMA_DIR}
    PATH_SUFFIXES
    include )

find_path(GCMMA_LIBRARY_DIRS
    NAMES
    libgcmma.a
    HINTS
    ${GCMMA_DIR}
    PATH_SUFFIXES
    lib
    lib64 )

find_library(GCMMA_LIBRARIES
    NAMES
    gcmma
    HINTS
    ${GCMMA_DIR}
    PATH_SUFFIXES
    lib
    lib64 )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GCMMA DEFAULT_MSG
    GCMMA_DIR
    GCMMA_INCLUDE_DIRS
    GCMMA_LIBRARY_DIRS
    GCMMA_LIBRARIES )

mark_as_advanced(GCMMA_DIR
    GCMMA_INCLUDE_DIRS
    GCMMA_LIBRARY_DIRS
    GCMMA_LIBRARIES )

