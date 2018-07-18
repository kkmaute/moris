# Optimization Dependencies -----------------------------------------------
# -------------------------------------------------------------------------

# Check if OPT has already been included
if(DEFINED OPT_CONFIGURED_ONCE)
    return()
endif()

set(OPT_CONFIGURED_ONCE "YES")

# Add OPT to the source directory list
list(APPEND MORIS_SRC_DIRS ${OPT})

# Include libraries needed by OPT
# some tpls needed
include(share/cmake/MatthewCMake/LNA_Depends.cmake)

set(OPT_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )
