#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Gperftools Find Module --------------------------------------------------
# -------------------------------------------------------------------------

set(GPERFTOOLS_ENV_VARS
    $ENV{GPERFTOOLSDIR}
    $ENV{GPERFTOOLS_DIR}
    $ENV{Gperftools_DIR}
    $ENV{gperftools_DIR}
    $ENV{GPERFTOOLS_ROOT}
    $ENV{Gperftools_ROOT}
    $ENV{gperftools_ROOT}
    $ENV{GPERFTOOLS_PATH}
    $ENV{Gperftools_PATH}
    $ENV{gperftools_PATH} )

find_path(GPERFTOOLS_DIR
    NAMES
    lib/libtcmalloc_and_profiler.a
    lib64/libtcmalloc_and_profiler.a
    lib/libtcmalloc_and_profiler.so
    lib64/libtcmalloc_and_profiler.so
    HINTS
    ${GPERFTOOLS_ENV_VARS} )

find_path(GPERFTOOLS_INCLUDE_DIRS
    NAMES
    gperftools/profiler.h
    HINTS
    ${GPERFTOOLS_DIR}
    PATH_SUFFIXES
    include )

find_path(GPERFTOOLS_LIBRARY_DIRS
    NAMES
    libtcmalloc_and_profiler.a
    libtcmalloc_and_profiler.so
    HINTS
    ${GPERFTOOLS_DIR}
    PATH_SUFFIXES
    lib
    lib64 )

find_library(GPERFTOOLS_LIBRARIES
    NAMES
    tcmalloc_and_profiler
    HINTS
    ${GPERFTOOLS_DIR}
    PATH_SUFFIXES
    lib
    lib64 )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GPERFTOOLS DEFAULT_MSG
    GPERFTOOLS_DIR
    GPERFTOOLS_INCLUDE_DIRS
    GPERFTOOLS_LIBRARY_DIRS
    GPERFTOOLS_LIBRARIES )

mark_as_advanced(GPERFTOOLS_DIR
    GPERFTOOLS_INCLUDE_DIRS
    GPERFTOOLS_LIBRARY_DIRS
    GPERFTOOLS_LIBRARIES )

