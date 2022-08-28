#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Time Solver Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if SOL_CORE has already been included
if(DEFINED SOL_CORE_CONFIGURED_ONCE)
    return()
endif()

set(SOL_CORE_CONFIGURED_ONCE "YES")

# Add SOL_CORE to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${SOL}/${SOL_CORE})

# Include libraries needed by SOL_CORE
set(SOL_CORE_TPL_DEPENDENCIES
    ""
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)

