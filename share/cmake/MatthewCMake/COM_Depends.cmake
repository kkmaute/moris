# Communications Dependencies ---------------------------------------------
# -------------------------------------------------------------------------

# Check if COM has already been included
if(DEFINED COM_CONFIGURED_ONCE)
    return()
endif()

set(COM_CONFIGURED_ONCE "YES")

# Add COM to the source directory list
list(APPEND MORIS_SRC_DIRS ${COM})

# Include libraries needed by COM
set(COM_TPL_DEPENDENCIES "")

include(share/cmake/MatthewCMake/LNA_Depends.cmake)

list(APPEND COM_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES}
    ${MATRIX_LIB} )
