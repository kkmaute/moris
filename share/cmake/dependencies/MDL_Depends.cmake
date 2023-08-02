#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Model Dependencies ------------------------------
# -------------------------------------------------------------------------

# Check if MDL has already been included
if(DEFINED MDL_CONFIGURED_ONCE)
    return()
endif()

set(MDL_CONFIGURED_ONCE "YES")

# Add MDL to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${FEM}/${MDL})

# Include libraries needed by MDL
set(MDL_TPL_DEPENDENCIES
    ""
    )

# Moris packages needed by MDL
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TSA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/SDF_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/VIS_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/PRM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)

# needed for tests
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/SDF_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)

