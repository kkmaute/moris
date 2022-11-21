#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

if (SUPERLU_DIST_LIBRARIES)
    set(SUPERLU_FIND_QUIETLY TRUE)
endif (SUPERLU_DIST_LIBRARIES)

set(SUPERLU_DIST_ENV_VARS
    $ENV{SUPERLUDISTDIR}
    $ENV{SUPERLU_DIST_DIR}
    $ENV{SuperLU_DIST_DIR}
    $ENV{SUPERLU_DIST_ROOT}
    $ENV{SuperLU_DIST_ROOT}
    $ENV{SUPERLU_DIST_PATH}
    $ENV{SuperLU_DIST_PATH} )

find_library(SUPERLU_DIST_LIBRARIES 
    NAMES
    superlu_dist
    HINTS
    ${SUPERLU_DIST_ENV_VARS}
    PATHS
    ${LIB_INSTALL_DIR}
    PATH_SUFFIXES
    lib
    lib64 )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SuperLU_DIST DEFAULT_MSG 
                                  SUPERLU_DIST_LIBRARIES)

mark_as_advanced(SUPERLU_DIST_LIBRARIES)

