# Model Dependencies ------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MOD has already been included
if(DEFINED MOD_CONFIGURED_ONCE)
    return()
endif()

set(MOD_CONFIGURED_ONCE "YES")

# Add MOD to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${MOD})

# Include libraries needed by MOD
# N/A
set(MOD_TPL_DEPENDENCIES
    ${ARMADILLO_EIGEN} )

include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)

list(APPEND MOD_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )
