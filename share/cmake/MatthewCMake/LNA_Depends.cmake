# Linear Algebra Dependencies ---------------------------------------------
# -------------------------------------------------------------------------

# Check if LNA has already been included
if(DEFINED LNA_CONFIGURED_ONCE)
    return()
endif()

set(LNA_CONFIGURED_ONCE "YES")

# Add LNA to the header directory list
list(APPEND MORIS_HEADER_DIRS ${LNA})

# Include libraries needed by LNA
# N/A
