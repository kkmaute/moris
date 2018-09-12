# XTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if XTK has already been included
if(DEFINED XTK_CONFIGURED_ONCE)
    return()
endif()

set(XTK_CONFIGURED_ONCE "YES")

# Add XTK to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${XTK})


# List moris projects directly needed by PROJ
set(LINALG_MORIS_DEPENDENCIES
    ${LINALG} )

# needs some tpls
set(XTK_TPL_DEPENDENCIES
    ${ARMADILLO_EIGEN}
    "superlu"
    "trilinos"
    )

