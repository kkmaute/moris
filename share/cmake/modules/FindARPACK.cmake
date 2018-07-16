# ARPACK Find Module ------------------------------------------------------
# -------------------------------------------------------------------------

set(ARPACK_ENV_VARS
    $ENV{ARPACKDIR}
    $ENV{ARPACK_DIR}
    $ENV{arpack_DIR}
    $ENV{ARPACK_ROOT}
    $ENV{arpack_ROOT}
    $ENV{ARPACK_PATH}
    $ENV{arpack_PATH} )

find_library(ARPACK_LIBRARIES
    NAMES
    arpack
    HINTS
    ${ARPACK_ENV_VARS}
    PATH_SUFFIXES
    lib
    lib64 )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK DEFAULT_MSG ARPACK_LIBRARIES)

mark_as_advanced(ARPACK_LIBRARIES)
