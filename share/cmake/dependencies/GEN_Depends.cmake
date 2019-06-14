# Geometry Engine Dependencies --------------------------------------------
# -------------------------------------------------------------------------

# Check if GEN has already been included
if(DEFINED GEN_CONFIGURED_ONCE)
    return()
endif()

set(GEN_CONFIGURED_ONCE "YES")

# Add GEN to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${GEN})

# Third party libraries used directly by GEN
set(GEN_LIB_TPL_DEPENDENCIES
	""
	)

# Additional third party libraries needed by tests
set(GEN_TEST_TPL_DEPENDENCIES
	""
	)

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake) #> headers
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MAP_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/SDF_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake) #> headers
#include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake) #> headers
#include(${MORIS_DEPENDS_DIR}/MAP_Depends.cmake) #> headers
#include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake) #> headers
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake) #> headers
#include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake) #> headers
#include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake) #> headers

# Includes needed for test
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

# Include third party libraries indirectly needed by GEN
#list(APPEND GEN_TPL_DEPENDENCIES
#    ${MTK_TPL_DEPENDENCIES}
#    ${LINALG_TPL_DEPENDENCIES}
#    ${ARMADILLO_EIGEN}
#    )
