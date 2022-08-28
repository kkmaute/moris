#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

if (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
    set(SUPERLU_FIND_QUIETLY TRUE)
endif (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

set(SUPERLU_ENV_VARS
    $ENV{SUPERLUDIR}
    $ENV{SUPERLU_DIR}
    $ENV{SuperLU_DIR}
    $ENV{SUPERLU_ROOT}
    $ENV{SuperLU_ROOT}
    $ENV{SUPERLU_PATH}
    $ENV{SuperLU_PATH} )

find_path(SUPERLU_INCLUDES
    NAMES
    supermatrix.h
    HINTS
    ${SUPERLU_ENV_VARS}
    PATHS
    ${INCLUDE_INSTALL_DIR}
    PATH_SUFFIXES
    superlu
    include
    SRC )

find_library(SUPERLU_LIBRARIES 
    NAMES
    superlu
    HINTS
    ${SUPERLU_ENV_VARS}
    PATHS
    ${LIB_INSTALL_DIR}
    PATH_SUFFIXES
    lib
    lib64 )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUPERLU DEFAULT_MSG
                                  SUPERLU_INCLUDES SUPERLU_LIBRARIES )

mark_as_advanced(SUPERLU_INCLUDES SUPERLU_LIBRARIES)

