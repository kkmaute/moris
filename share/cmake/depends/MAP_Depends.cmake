# MTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MAP has already been included
if(DEFINED MAP_CONFIGURED_ONCE)
    return()
endif()

set(MAP_CONFIGURED_ONCE "YES")

# Add MTK to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${MTK}/${MAP})

# Include libraries needed by HMR
set(MAP_TPL_DEPENDENCIES
    "superlu"
    "boost"
    )
    
# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/FEM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)

# Include third party libraries indirectly needed by MAP
list(APPEND MAP_TPL_DEPENDENCIES
    ${HMR_TPL_DEPENDENCIES}
    )