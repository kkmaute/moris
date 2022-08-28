#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Example Dependencies ------------------------------
# ---------------------------------------------------

# Check if MAIN has already been included
if(DEFINED EXA_CONFIGURED_ONCE)
    return()
endif()

set(EXA_CONFIGURED_ONCE "YES")

# Add EXA to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${EXA})

# Include all libraries for EXA
set(MAIN_TPL_DEPENDENCIES
    "trilinos"
    )

# Moris packages needed by main
include(${MORIS_DEPENDS_DIR}/main_includes.cmake)
