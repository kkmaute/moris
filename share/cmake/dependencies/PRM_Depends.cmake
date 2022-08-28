#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Parameters Dependencies ---------------------------------------------
# -------------------------------------------------------------------------

# Check if PRM has already been included
if(DEFINED PRM_CONFIGURED_ONCE)
    return()
endif()

set(PRM_CONFIGURED_ONCE "YES")

# Add PRM to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${PRM})

# Third party libraries needed directly
set(PRM_TPL_DEPENDENCIES
	"mpi" )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/VIS_Depends.cmake)

