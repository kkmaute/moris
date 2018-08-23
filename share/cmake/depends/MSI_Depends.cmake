# MSI Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MSI has already been included
if(DEFINED MSI_CONFIGURED_ONCE)
    return()
endif()

set(MSI_CONFIGURED_ONCE "YES")

# Add MSI to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${FEM}/${MSI})

# Include libraries needed by MSI
set(MSI_TPL_DEPENDENCIES
    #"trilinos"
    ${ACML_LAPACK_MKL}
    ${ARMADILLO_EIGEN} #> used for hierarchical
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/MOD_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TOL_Depends.cmake)

include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake) #> headers
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)

# Include third party libraries indirectly needed by MSI
list(APPEND MSI_TPL_DEPENDENCIES
    ${MOD_TPL_DEPENDENCIES}
    ${TOL_TPL_DEPENDENCIES}
    ${LNA_TPL_DEPENDENCIES}
    ${DLA_TPL_DEPENDENCIES}
    )

