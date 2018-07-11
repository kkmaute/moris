# MTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MTK has already been included
if(DEFINED MTK_CONFIGURED_ONCE)
    return()
endif()

set(MTK_CONFIGURED_ONCE "YES")

# Add MTK to the header directory list
list(APPEND MORIS_HEADER_DIRS ${MTK})

# Include libraries needed by MTK
# N/A
