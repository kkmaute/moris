#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# LBFGSB Find Module -------------------------------------------------------
# -------------------------------------------------------------------------

set(LBFGSB_ENV_VARS
    $ENV{LBFGSBDIR}
    $ENV{LBFGSB_DIR}
    $ENV{lbfgsb_DIR}
    $ENV{LBFGSB_ROOT}
    $ENV{lbfgsb_ROOT}
    $ENV{LBFGSB_PATH}
    $ENV{lbfgsb_PATH} )

find_path(LBFGSB_LIBRARY_DIRS
    NAMES
    liblbfgsb.a
    PATHS
    ${LBFGSB_ENV_VARS}
    PATH_SUFFIXES
    lib64
    lib )

find_library(LBFGSB_LIBRARIES
    NAMES
    lbfgsb
    HINTS
    ${LBFGSB_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LBFGSB DEFAULT_MSG 
    LBFGSB_LIBRARY_DIRS
    LBFGSB_LIBRARIES )

mark_as_advanced(LBFGSB_LIBRARY_DIRS LBFGSB_LIBRARIES)


