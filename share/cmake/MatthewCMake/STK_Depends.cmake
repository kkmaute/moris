# STK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if STK has already been included
if(DEFINED STK_CONFIGURED_ONCE)
    return()
endif()

set(STK_CONFIGURED_ONCE "YES")

# Add STK to the source directory list
list(APPEND MORIS_SRC_DIRS ${STK})

# Include libraries needed by STK
# needs some tpls
# also hierarchical
set(STK_TPL_DEPENDENCIES
    #"trilinos"
    "boost" #> used for hierarchical
    ${MATH_LIB}
    ${MATRIX_LIB} #> used for hierarchical
     )

include(${SHARE}/${CMAKE}/MatthewCMake/MOD_Depends.cmake)
include(${SHARE}/${CMAKE}/MatthewCMake/TOL_Depends.cmake)

include(share/cmake/MatthewCMake/LNA_Depends.cmake) #> headers
include(share/cmake/MatthewCMake/DLA_Depends.cmake)

list(APPEND STK_TPL_DEPENDENCIES
    ${MOD_TPL_DEPENDENCIES}
    ${TOL_TPL_DEPENDENCIES}
    ${LNA_TPL_DEPENDENCIES}
    ${DLA_TPL_DEPENDENCIES})

