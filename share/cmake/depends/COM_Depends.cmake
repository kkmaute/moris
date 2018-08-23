# Communications Dependencies ---------------------------------------------
# -------------------------------------------------------------------------

# Check if COM has already been included
if(DEFINED COM_CONFIGURED_ONCE)
    return()
endif()

set(COM_CONFIGURED_ONCE "YES")

# Add COM to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${COM})

# Include libraries needed by COM
set(COM_TPL_DEPENDENCIES
    ${ARMADILLO_EIGEN}
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)

# Include third party libraries indirectly needed by COM
list(APPEND COM_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES}
    )
