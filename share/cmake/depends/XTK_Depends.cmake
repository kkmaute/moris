# XTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if XTK has already been included
if(DEFINED XTK_CONFIGURED_ONCE)
    return()
endif()

set(XTK_CONFIGURED_ONCE "YES")

# Add XTK to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${XTK})

# Include libraries needed by XTK

