# Tools Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if TOL has already been included
if(DEFINED TOL_CONFIGURED_ONCE)
    return()
endif()

set(TOL_CONFIGURED_ONCE "YES")

# Add TOL to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${TOL})

# Include third party libraries needed by TOL
set(TOL_TPL_DEPENDENCIES
    ${ACML_LAPACK_MKL}
    ${ARMADILLO_EIGEN}
    "superlu"
    "trilinos"
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)

# Include third party libraries indirectly needed by TOL
list(APPEND TOL_TPL_DEPENDENCIES
	${CON_TPL_DEPENDENCIES}
    ${LINALG_TPL_DEPENDENCIES}
    )
