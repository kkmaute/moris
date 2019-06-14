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
set(DLA_LIB_TPL_DEPENDENCIES
    "PETSc"
    #"trilinos" #
    #${ACML_LAPACK_MKL} #
    "mpi"
    #${ARMADILLO_EIGEN} #
    "superlu"
    )

# Additional third party libraries for tests
set(DLA_TEST_TPL_DEPENDENCIES
	${ARMADILLO_EIGEN}
	)

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)

# Maybe just needed for tests? #!#
#include(${MORIS_DEPENDS_DIR}/SOL_CORE_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/MAP_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/SDF_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

# Include third party libraries indirectly needed by DLA
#list(APPEND DLA_TPL_DEPENDENCIES
#    ${LINALG_TPL_DEPENDENCIES}
#    ${SOL_CORE_TPL_DEPENDENCIES}
#    )
