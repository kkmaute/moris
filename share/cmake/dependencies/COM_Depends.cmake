#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Communications Dependencies ---------------------------------------------
# -------------------------------------------------------------------------

# Check if COM has already been included
if(DEFINED COM_CONFIGURED_ONCE)
    return()
endif()

set(COM_CONFIGURED_ONCE "YES")

# Add COM to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${COM})

# Third party libraries needed directly
set(COM_TPL_DEPENDENCIES
	"mpi" )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)

