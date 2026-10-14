# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# MIG Dependencies

# Check if MIG has already been included.
if(DEFINED MIG_CONFIGURED_ONCE)
    return()
endif()

set(MIG_CONFIGURED_ONCE "YES")

# Add MIG to the source directory list.
list(APPEND MORIS_SOURCE_DIRS ${MIG})

# Third-party libraries directly used by MIG.
set(MIG_TPL_DEPENDENCIES "")

# Include required MORIS modules.
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/GEN_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TOL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)

#needed for testing
include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)
