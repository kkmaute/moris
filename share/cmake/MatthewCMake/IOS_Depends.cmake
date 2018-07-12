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
