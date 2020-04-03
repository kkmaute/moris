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
include(${MORIS_DEPENDS_DIR}/GEN_CORE_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MAP_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/SDF_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)

# Includes needed for test
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
