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
# PETSc and Trilinos; add later
# include(${MORIS_CMAKE_DIR}/PETSc.cmake)
set(INT_TPL_DEPENDENCIES
    ${ARMADILLO_EIGEN}
    "superlu" #Armadillo
    ${ACML_LAPACK_MKL} #SuperLU
    )

include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)

list(APPEND INT_TPL_DEPENDENCIES
     ${LNA_TPL_DEPENDENCIES}
     ${MSI_TPL_DEPENDENCIES} )
