#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Main Dependencies ---------------------------------
# ---------------------------------------------------

# Check if MAIN has already been included
if(DEFINED MAIN_CONFIGURED_ONCE)
    return()
endif()

set(MAIN_CONFIGURED_ONCE "YES")

# Add MAIN to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${MAIN})

# Include all libraries for MAIN
set(MAIN_TPL_DEPENDENCIES
    "trilinos"
    )

# Moris packages needed by main
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/GEN_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/MOD_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/OPT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/PRM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/SOL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/STK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TOL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TSA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/VIS_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/WRK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MIG_Depends.cmake)

