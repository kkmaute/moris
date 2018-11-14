# Model Dependencies ------------------------------
# -------------------------------------------------------------------------

# Check if MDL has already been included
if(DEFINED MDL_CONFIGURED_ONCE)
    return()
endif()

set(MDL_CONFIGURED_ONCE "YES")

# Add MDL to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${FEM}/${MDL})

# Include libraries needed by MDL
# PETSc and Trilinos; add later
# include(${MORIS_CMAKE_DIR}/PETSc.cmake)
set(MDL_TPL_DEPENDENCIES
    ${ARMADILLO_EIGEN}
    "superlu" #Armadillo
    ${ACML_LAPACK_MKL} #SuperLU
    )
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)

list(APPEND MDL_TPL_DEPENDENCIES
     ${LINALG_TPL_DEPENDENCIES}
     ${MSI_TPL_DEPENDENCIES} 
     ${INT_TPL_DEPENDENCIES} )