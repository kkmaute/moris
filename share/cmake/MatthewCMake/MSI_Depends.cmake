# MSI Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MSI has already been included
if(DEFINED MSI_CONFIGURED_ONCE)
    return()
endif()

set(MSI_CONFIGURED_ONCE "YES")

# Add MSI to the source directory list
list(APPEND MORIS_SRC_DIRS ${FEM}/${MSI})

# Include libraries needed by MSI
# needs some tpls
# also hierarchical
set(MSI_TPL_DEPENDENCIES
    #"trilinos"
    "boost" #> used for hierarchical
    ${MATH_LIB}
    ${MATRIX_LIB} #> used for hierarchical
     )

include(${SHARE}/${CMAKE}/MatthewCMake/MOD_Depends.cmake)
include(${SHARE}/${CMAKE}/MatthewCMake/TOL_Depends.cmake)

include(share/cmake/MatthewCMake/LNA_Depends.cmake) #> headers
include(share/cmake/MatthewCMake/DLA_Depends.cmake)

list(APPEND MSI_TPL_DEPENDENCIES
    ${MOD_TPL_DEPENDENCIES}
    ${TOL_TPL_DEPENDENCIES}
    ${LNA_TPL_DEPENDENCIES}
    ${DLA_TPL_DEPENDENCIES})

