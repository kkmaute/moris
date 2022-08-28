#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# EIGEN Find Module ------------------------------------------------------
# -------------------------------------------------------------------------

set(EIGEN_ENV_VARS
    $ENV{EIGENDIR}
    $ENV{EIGEN_DIR}
    $ENV{Eigen_DIR}
    $ENV{Eigen3_DIR}
    $ENV{EIGEN_ROOT}
    $ENV{Eigen_ROOT}
    $ENV{Eigen3_ROOT}
    $ENV{EIGEN_PATH}
    $ENV{Eigen_PATH}
    $ENV{Eigen3_PATH})

find_path(EIGEN_DIRS
    NAMES
    Eigen
    HINTS
    ${EIGEN_ENV_VARS}
    PATH_SUFFIXES
    include/Eigen )


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EIGEN DEFAULT_MSG EIGEN_DIRS)

mark_as_advanced(EIGEN_DIRS)

