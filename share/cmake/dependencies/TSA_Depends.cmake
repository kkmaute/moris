#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Time Solver Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if TSA has already been included
if(DEFINED TSA_CONFIGURED_ONCE)
    return()
endif()

set(TSA_CONFIGURED_ONCE "YES")

# Add TSA to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${SOL}/${TSA})

# Include libraries needed by TSA
set(TSA_TPL_DEPENDENCIES
    ""
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/SOL_CORE_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)

# Includes needed for tests
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

