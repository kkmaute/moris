# Distributed Linear Algebra Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if DLA has already been included
if(DEFINED DLA_CONFIGURED_ONCE)
    return()
endif()

set(DLA_CONFIGURED_ONCE "YES")

# Add DLA to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${SOL}/${DLA})

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
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)

# Maybe just needed for tests? #!#
include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MAP_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/SDF_Depends.cmake)

# Include third party libraries indirectly needed by DLA
list(APPEND DLA_TPL_DEPENDENCIES
    ${LINALG_TPL_DEPENDENCIES}
    )
