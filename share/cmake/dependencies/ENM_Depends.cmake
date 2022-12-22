#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Enums Dependencies ---------------------------------------------
# -------------------------------------------------------------------------

# Check if ENM has already been included
if(DEFINED ENM_CONFIGURED_ONCE)
    return()
endif()

set(ENM_CONFIGURED_ONCE "YES")

# Add ENM to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${PRM}/${ENM})

# Third party libraries needed directly
set(ENM_TPL_DEPENDENCIES
	"" )

# Make sure needed moris libraries are built
# include(${MORIS_DEPENDS_DIR}/MRS_Depends.cmake)

