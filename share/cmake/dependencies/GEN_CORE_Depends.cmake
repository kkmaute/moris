#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Time Solver Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if GEN_CORE has already been included
if(DEFINED GEN_CORE_CONFIGURED_ONCE)
    return()
endif()

set(GEN_CORE_CONFIGURED_ONCE "YES")

# Add GEN_CORE to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${GEN}/${GEN_CORE})

# Include libraries needed by GEN_CORE
set(GEN_CORE_TPL_DEPENDENCIES
    ""
    )


# Make sure needed moris libraries are built


