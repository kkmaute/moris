# Linear Algebra Dependencies ---------------------------------------------
# -------------------------------------------------------------------------

# Check if LNA has already been included
if(DEFINED LNA_CONFIGURED_ONCE)
    return()
endif()

set(LNA_CONFIGURED_ONCE "YES")

# Add LNA to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${LNA})

# Include libraries needed by LNA
# Some tpls
set(LNA_TPL_DEPENDENCIES
    "viennacl"
    ${ACML_LAPACK_MKL}
    ${ARMADILLO_EIGEN}
    "superlu" )

include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

list(APPEND LNA_TPL_DEPENDENCIES
    ${ALG_TPL_DEPENDENCIES} )
