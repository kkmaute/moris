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
    PATHS
    ${MKL_ENV_VARS}
    /usr/lib/mkl
    /usr/lib
    /usr
    /usr/local )

find_path(MKL_INCLUDE_DIRS
    NAMES
    mkl.h
    PATHS
    ${MKL_DIR}
    ${MKL_ENV_VARS}
    PATH_SUFFIXES
    include )

find_library(MKL_lp64
    NAMES
    mkl_intel_lp64
    PATHS
    ${MKL_DIR}/lib/intel64
    ${MKL_DIR}/lib/intel64_lin )

find_library(MKL_core
    NAMES
    mkl_core
    PATHS
    ${MKL_DIR}/lib/intel64
    ${MKL_DIR}/lib/intel64_lin )

find_library(MKL_sequential
    NAMES
    mkl_sequential
    PATHS
    ${MKL_DIR}/lib/intel64
    ${MKL_DIR}/lib/intel64_lin )

find_library(MKL_pthread
    NAMES
    pthread
    PATHS
    ${MKL_DIR}/lib/intel64
    ${MKL_DIR}/lib/intel64_lin )

set(MKL_LIBRARIES
    ${MKL_lp64}
    ${MKL_core}
    ${MKL_sequential}
    ${MKL_pthread} )


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARIES)

mark_as_advanced(MKL_LIBRARIES)
