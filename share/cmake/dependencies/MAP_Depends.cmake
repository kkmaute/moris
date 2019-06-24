# MTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MAP has already been included
if(DEFINED MAP_CONFIGURED_ONCE)
    return()
endif()

set(MAP_CONFIGURED_ONCE "YES")

# Add MTK to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${MTK}/${MAP})

# Third party libraries directly used by MAP library
set(MAP_TPL_DEPENDENCIES
    ""
    )
    
# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TSA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
