# Time Solver Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if SOL_CORE has already been included
if(DEFINED SOL_CORE_CONFIGURED_ONCE)
    return()
endif()

set(SOL_CORE_CONFIGURED_ONCE "YES")

# Add SOL_CORE to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${SOL}/${SOL_CORE})

# Include libraries needed by SOL_CORE
set(SOL_CORE_TPL_DEPENDENCIES
	# ${ARMADILLO_EIGEN} # for example
    )

# Make sure needed moris libraries are built
#include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake) # for example

# Include third party libraries indirectly needed by SOL_CORE
list(APPEND SOL_CORE_TPL_DEPENDENCIES
	#${LINALG_TPL_DEPENDENCIES} # for example
    )
