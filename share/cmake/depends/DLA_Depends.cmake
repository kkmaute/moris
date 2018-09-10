# Distributed Linear Algebra Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if DLA has already been included
if(DEFINED DLA_CONFIGURED_ONCE)
    return()
endif()

set(DLA_CONFIGURED_ONCE "YES")

# Add DLA to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${DLA})

# Include libraries needed by DLA
set(DLA_TPL_DEPENDENCIES
    "PETSc"
    "trilinos"
    ${ACML_LAPACK_MKL}
    "mpi"
    ${ARMADILLO_EIGEN}
    "superlu"
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)

# Include third party libraries indirectly needed by DLA
list(APPEND DLA_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES}
    )
