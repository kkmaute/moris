# Geometry Engine Dependencies --------------------------------------------
# -------------------------------------------------------------------------

# Check if GEN has already been included
if(DEFINED GEN_CONFIGURED_ONCE)
    return()
endif()

set(GEN_CONFIGURED_ONCE "YES")

# Add GEN to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${GEN}/${GEN_MAIN})

# Third party libraries used directly by GEN
set(GEN_TPL_DEPENDENCIES
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
include(${MORIS_DEPENDS_DIR}/WRK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)

# Includes needed for test
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)


