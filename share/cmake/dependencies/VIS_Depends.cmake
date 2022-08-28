#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# VIS Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if VIS has already been included
if(DEFINED VIS_CONFIGURED_ONCE)
    return()
endif()

set(VIS_CONFIGURED_ONCE "YES")

# Add VIS to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${FEM}/${VIS})

# Include libraries needed by VIS
set(VIS_TPL_DEPENDENCIES
    ""
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)

# needed includes for test
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)

