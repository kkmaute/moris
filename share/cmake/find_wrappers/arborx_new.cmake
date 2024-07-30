#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# ArborX library ----------------------------------------------------------
# -------------------------------------------------------------------------

if (NOT ARBORX_FOUND_ONCE)
    set(ARBORX_ENV_VARS
            $ENV{ARBX_DIR}
            $ENV{ARBORXDIR}
            $ENV{ARBORX_DIR}
            $ENV{ArborX_DIR}
            $ENV{ARBORX_ROOT}
            $ENV{ArborX_ROOT}
            $ENV{ARBORX_PATH}
            $ENV{ArborX_PATH}
            $ENV{ARBORX_INSTALL_DIR}
    )
    # add installation directory ${ARBORX_INSTALL_DIR} to the cmake search path
    list(APPEND CMAKE_PREFIX_PATH ${ARBORX_INSTALL_DIR})

    find_package(ArborX REQUIRED HINTS ${ARBORX_ENV_VARS})

    set(MORIS_ARBORX_INCLUDE_DIRS ${ARBORX_DIRS})
    set(MORIS_ARBORX_DEFINITIONS "-DMORIS_HAVE_ARBORX")

    mark_as_advanced(ArborX_DIR
            MORIS_ARBORX_INCLUDE_DIRS)

    if (ArborX_FOUND)
        set(ARBORX_FOUND_ONCE TRUE)
    endif ()

#    message(STATUS "ARBORX_INCLUDE_DIRS: ${ARBORX_DIRS}")
endif ()

if (TARGET ArborX::ArborX AND NOT TARGET ${MORIS}::arborx)
    add_library(${MORIS}::arborx ALIAS ArborX::ArborX)
endif ()

set(ARBORX_FOUND_ONCE TRUE)

