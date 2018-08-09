# Tools Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if TOL has already been included
if(DEFINED TOL_CONFIGURED_ONCE)
    return()
endif()

set(TOL_CONFIGURED_ONCE "YES")

# Add TOL to the source directory list
list(APPEND MORIS_SRC_DIRS ${TOL})

# Include libraries needed by TOL
# needs some tpls
set(TOL_TPL_DEPENDENCIES
    ${MATH_LIB}
    ${MATRIX_LIB}
    "superlu" )

include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)

list(APPEND TOL_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )
