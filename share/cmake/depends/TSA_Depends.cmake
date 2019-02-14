# Time Solver Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if TSA has already been included
if(DEFINED TSA_CONFIGURED_ONCE)
    return()
endif()

set(TSA_CONFIGURED_ONCE "YES")

# Add TSA to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${SOL}/${TSA})

# Include libraries needed by TSA
set(TSA_TPL_DEPENDENCIES
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
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)

# Include third party libraries indirectly needed by TSA
list(APPEND TSA_TPL_DEPENDENCIES
	${SOL_CORE_TPL_DEPENDENCIES}
    ${LINALG_TPL_DEPENDENCIES}
    ${DLA_TPL_DEPENDENCIES}
    ${NLA_TPL_DEPENDENCIES}
    )
