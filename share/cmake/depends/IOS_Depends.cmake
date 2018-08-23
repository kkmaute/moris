# IOS Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if IOS has already been included
if(DEFINED IOS_CONFIGURED_ONCE)
    return()
endif()

set(IOS_CONFIGURED_ONCE "YES")

# Add IOS to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${IOS})

# Include libraries needed by IOS
set(IOS_TPL_DEPENDENCIES
    "boost"
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)

# Include third party libraries indirectly needed by IOS
list(APPEND IOS_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES}
    )
