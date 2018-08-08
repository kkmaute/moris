# Integration and Interpolation Dependencies ------------------------------
# -------------------------------------------------------------------------

# Check if INT has already been included
if(DEFINED INT_CONFIGURED_ONCE)
    return()
endif()

set(INT_CONFIGURED_ONCE "YES")

# Add INT to the source directory list
list(APPEND MORIS_SRC_DIRS ${FEM}/${INT})

# Include libraries needed by INT
# PETSc and Trilinos; add later
# include(share/cmake/PETSc.cmake)
set(INT_TPL_DEPENDENCIES
    ${MATRIX_LIB}
    "superlu" #Armadillo
    ${MATH_LIB} #SuperLU
    )

include(share/cmake/MatthewCMake/LNA_Depends.cmake)
include(share/cmake/MatthewCMake/INT_Depends.cmake)
include(share/cmake/MatthewCMake/MSI_Depends.cmake)

list(APPEND INT_TPL_DEPENDENCIES
     ${LNA_TPL_DEPENDENCIES}
     ${MSI_TPL_DEPENDENCIES} )
