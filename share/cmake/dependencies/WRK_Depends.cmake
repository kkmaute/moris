# WRK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if WRK has already been included
if(DEFINED WRK_CONFIGURED_ONCE)
    return()
endif()

set(WRK_CONFIGURED_ONCE "YES")

# Add WRK to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${WRK})

# Include libraries needed by HMR
set(WRK_TPL_DEPENDENCIES
    ""
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/GEN_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MAP_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MIG_Depends.cmake)

# added as temp fix for hmr exe, test, and tutorials
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
