# STK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if STK has already been included
if(DEFINED STK_CONFIGURED_ONCE)
    return()
endif()

set(STK_CONFIGURED_ONCE "YES")

# Add STK to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${STK})

# Include libraries needed by STK
set(STK_TPL_DEPENDENCIES
    #"trilinos"
    "boost" #> used for hierarchical
    ${ACML_LAPACK_MKL}
    ${ARMADILLO_EIGEN} #> used for hierarchical
    "superlu"
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/MOD_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TOL_Depends.cmake)

include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake) #> headers
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)

# Include third party libraries indirectly needed by STK
list(APPEND STK_TPL_DEPENDENCIES
    ${MOD_TPL_DEPENDENCIES}
    ${TOL_TPL_DEPENDENCIES}
    ${LNA_TPL_DEPENDENCIES}
    ${DLA_TPL_DEPENDENCIES}
    )
