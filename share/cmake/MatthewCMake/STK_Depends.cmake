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

