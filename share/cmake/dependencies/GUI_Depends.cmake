#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# GUI Dependencies ------------------------------
# -------------------------------------------------------------------------

# Check if GUI has already been included
if(DEFINED GUI_CONFIGURED_ONCE)
    return()
endif()

set(GUI_CONFIGURED_ONCE "YES")

# Add GUI to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${GUI})

# Third party libraries needed by GUI library
set(GUI_TPL_DEPENDENCIES
    ""
    )

# Moris packages needed by GUI
include(${MORIS_DEPENDS_DIR}/ENM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/FEM_Depends.cmake)
