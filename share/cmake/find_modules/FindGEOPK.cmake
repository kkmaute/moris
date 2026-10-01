#
# GEOPK Find Module
# -------------------------------------------------------------------------
set(GEOPK_ENV_VARS
    $ENV{GEOPK_DIR}
    $ENV{GEOPK_PATH} )

find_path(GEOPK_LIBRARY_DIRS
    NAMES
        libgeompack3.gcc.a
    PATHS
        ${GEOPK_ENV_VARS} )

find_library(GEOPK_LIBRARIES
    NAMES
        geompack3.gcc
    HINTS
        ${GEOPK_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
    GEOPK
    DEFAULT_MSG
    GEOPK_LIBRARY_DIRS
    GEOPK_LIBRARIES )

mark_as_advanced(
    GEOPK_LIBRARY_DIRS
    GEOPK_LIBRARIES )