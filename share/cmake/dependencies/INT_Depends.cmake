# Integration and Interpolation Dependencies ------------------------------
# -------------------------------------------------------------------------

# Check if INT has already been included
if(DEFINED INT_CONFIGURED_ONCE)
    return()
endif()

set(INT_CONFIGURED_ONCE "YES")

# Add INT to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${FEM}/${INT})

# Include libraries needed by INT
set(INT_LIB_TPL_DEPENDENCIES
    #${ARMADILLO_EIGEN}
    #"superlu" #Armadillo
    #${ACML_LAPACK_MKL} #SuperLU
    )

# Additional third party libraries needed for test
set(INT_TEST_TPL_DEPENDENCIES
	#${ARMADILLO_EIGEN} #> remove when possible... i think
	)

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)

# included for tests
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

# Include third party libraries indirectly needed by INT
#list(APPEND INT_TPL_DEPENDENCIES
#     ${LINALG_TPL_DEPENDENCIES}
#     ${MTK_TPL_DEPENDENCIES}
#     ${DLA_TPL_DEPENDENCIES}
#     )
