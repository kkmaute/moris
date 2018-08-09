# IOS Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if IOS has already been included
if(DEFINED IOS_CONFIGURED_ONCE)
    return()
endif()

set(IOS_CONFIGURED_ONCE "YES")

# Add IOS to the source directory list
list(APPEND MORIS_SRC_DIRS ${IOS})

# Include libraries needed by IOS
# N/A
set(IOS_TPL_DEPENDENCIES
    "boost" )

include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)

list(APPEND IOS_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )
