#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# FEM Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if FEM has already been included
if(DEFINED FEM_CONFIGURED_ONCE)
    return()
endif()

set(FEM_CONFIGURED_ONCE "YES")

# FEM is only a wrapper around these moris packages (Aug. 20, 2018)
set(BUILD_INT ON CACHE BOOL "Build the integration and interpolation executable." FORCE)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)

set(BUILD_MDL ON CACHE BOOL "Build the model executable." FORCE)
include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)

set(BUILD_VIS ON CACHE BOOL "Build the model executable." FORCE)
include(${MORIS_DEPENDS_DIR}/VIS_Depends.cmake)

set(BUILD_MSI ON CACHE BOOL "Build the model solver interface executable." FORCE)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)


