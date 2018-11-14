# SDF Dependencies ------------------------------
# -------------------------------------------------------------------------

# Check if SDF has already been included
if(DEFINED SDF_CONFIGURED_ONCE)
    return()
endif()

set(SDF_CONFIGURED_ONCE "YES")

# Add SDF to the source directory list
#list(APPEND MORIS_SOURCE_DIRS ${GEN}/${SDF})

# Include libraries needed by SDF
# PETSc and Trilinos; add later
# include(${MORIS_CMAKE_DIR}/PETSc.cmake)
set(SDF_TPL_DEPENDENCIES
    ${ARMADILLO_EIGEN}
    "superlu" #Armadillo
    ${ACML_LAPACK_MKL} #SuperLU
    )
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)

list(APPEND SDF_TPL_DEPENDENCIES
	"boost"
#	${CON_TPL_DEPENDENCIES}
     ${LINALG_TPL_DEPENDENCIES} )