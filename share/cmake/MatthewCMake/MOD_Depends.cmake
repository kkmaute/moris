# Model Dependencies ------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MOD has already been included
if(DEFINED MOD_CONFIGURED_ONCE)
    return()
endif()

set(MOD_CONFIGURED_ONCE "YES")

# Add MOD to the source directory list
list(APPEND MORIS_SRC_DIRS ${MOD})

# Include libraries needed by MOD
# N/A
set(MOD_TPL_DEPENDENCIES
    ${MATRIX_LIB} )

include(share/cmake/MatthewCMake/LNA_Depends.cmake)

list(APPEND MOD_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )
