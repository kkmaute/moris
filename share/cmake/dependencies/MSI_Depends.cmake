# MSI Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MSI has already been included
if(DEFINED MSI_CONFIGURED_ONCE)
    return()
endif()

set(MSI_CONFIGURED_ONCE "YES")

# Add MSI to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${FEM}/${MSI})

# Include libraries needed by MSI
set(MSI_LIB_TPL_DEPENDENCIES
    #${ACML_LAPACK_MKL}
    #${ARMADILLO_EIGEN} #> used for hierarchical
    )

# Additional third party libraries needed for tests\
set(MSI_TEST_TPL_DEPENDENCIES
	""
	)

# Make sure needed moris libraries are built
#include(${MORIS_DEPENDS_DIR}/MOD_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/TOL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)

# needed includes for test
#include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)

# Include third party libraries indirectly needed by MSI
#list(APPEND MSI_TPL_DEPENDENCIES
#    ${MOD_TPL_DEPENDENCIES}
#    ${TOL_TPL_DEPENDENCIES}
#    ${LINALG_TPL_DEPENDENCIES}
#    ${DLA_TPL_DEPENDENCIES}
#    )

