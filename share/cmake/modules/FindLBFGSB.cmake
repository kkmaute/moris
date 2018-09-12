# LBFGSB Find Module ------------------------------------------------------
# -------------------------------------------------------------------------

set(LBFGSB_ENV_VARS
    $ENV{LBFGSBDIR}
    $ENV{LBFGSB_DIR}
    $ENV{lbfgsb_DIR}
    $ENV{LBFGSB_ROOT}
    $ENV{lbfgsb_ROOT}
    $ENV{LBFGSB_PATH}
    $ENV{lbfgsb_PATH}
    $ENV{L-BFGS-B_DIR}
    $ENV{L-BFGS-B_ROOT}
    $ENV{L-BFGS-B_PATH} )

find_path(LBFGSB_LIBRARY_DIRS
    NAMES
    liblbfgsb.a
    HINTS
    ${LBFGSB_ENV_VARS}
    PATH_SUFFIXES
    lib
    lib64 )

find_library(LBFGSB_LIBRARIES
    NAMES
    lbfgsb
    HINTS
    ${LBFGSB_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LBFGSB DEFAULT_MSG
    LBFGSB_LIBRARY_DIRS
    LBFGSB_LIBRARIES)

mark_as_advanced(LBFGSB_LIBRARY_DIRS LBFGSB_LIBRARIES)
