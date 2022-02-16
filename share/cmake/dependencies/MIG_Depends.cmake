# MIG Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MAP has already been included
if(DEFINED MIG_CONFIGURED_ONCE)
    return()
endif()

set(MIG_CONFIGURED_ONCE "YES")

# Add MTK to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${MIG})

# Third party libraries directly used by MAP library
set(MAP_TPL_DEPENDENCIES
    ""
    )
    
# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/GEN_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TOL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)

#needed for testing
include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)
