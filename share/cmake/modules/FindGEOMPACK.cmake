# GEOMPACK Find Module ------------------------------------------------------
# -------------------------------------------------------------------------

set(GEOMPACK_ENV_VARS
    $ENV{GEOMPACKDIR}
    $ENV{GEOMPACK_DIR}
    $ENV{geompack_DIR}
    $ENV{GEOMPACK_ROOT}
    $ENV{geompack_ROOT}
    $ENV{GEOMPACK_PATH}
    $ENV{geompack_PATH} )

find_library(GEOMPACK_LIBRARIES
    NAMES
    geompack
    geompack2
    geompack3
    HINTS
    ${GEOMPACK_ENV_VARS}
    PATH_SUFFIXES
    lib
    lib64 )

# if(${GEOMPACK_LIBRARIES} NOT STREQUAL "libgeompack3

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GEOMPACK DEFAULT_MSG GEOMPACK_LIBRARIES)

mark_as_advanced(GEOMPACK_LIBRARIES)
