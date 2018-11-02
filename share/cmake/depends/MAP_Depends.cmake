# MTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MAP has already been included
if(DEFINED MAP_CONFIGURED_ONCE)
    return()
endif()

set(MAP_CONFIGURED_ONCE "YES")

# Add MTK to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${MTK}/${MAP})


# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/FEM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)

# Include third party libraries indirectly needed by MTK
list(APPEND MAP_TPL_DEPENDENCIES
    )