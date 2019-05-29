# NonLinear Algebra Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if NLA has already been included
if(DEFINED NLA_CONFIGURED_ONCE)
    return()
endif()

set(NLA_CONFIGURED_ONCE "YES")

# Add NLA to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${SOL}/${NLA})

# Include libraries needed by NLA
set(NLA_TPL_DEPENDENCIES
    "PETSc"
    "trilinos"
    ${ACML_LAPACK_MKL}
    "mpi"
    ${ARMADILLO_EIGEN}
    "superlu"
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/SOL_CORE_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)

# Include third party libraries indirectly needed by NLA
list(APPEND NLA_TPL_DEPENDENCIES
    ${LINALG_TPL_DEPENDENCIES}
    ${DLA_TPL_DEPENDENCIES}
    ${SOL_CORE_TPL_DEPENDENCIES}
    )
