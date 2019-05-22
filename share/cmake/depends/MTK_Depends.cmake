# MTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MTK has already been included
if(DEFINED MTK_CONFIGURED_ONCE)
    return()
endif()

set(MTK_CONFIGURED_ONCE "YES")

# Add MTK to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${MTK})

# Include libraries needed by MTK
set(MTK_TPL_DEPENDENCIES
    ${ARMADILLO_EIGEN}
    "superlu"
    "trilinos"
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)

set(MTK_PROJ ${MORIS_PACKAGE_DIR}/${MTK}/src/stk_impl)


include_directories(${MTK_PROJ})

# Include third party libraries indirectly needed by MTK
list(APPEND MTK_TPL_DEPENDENCIES
    )
