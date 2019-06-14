# Algorithms Dependencies -------------------------------------------------
# -------------------------------------------------------------------------

# Check if ALG has already been included
if(DEFINED ALG_CONFIGURED_ONCE)
    return()
endif()

set(ALG_CONFIGURED_ONCE "YES")

# Add ALG to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${ALG})

# Include libraries needed by ALG
set(ALG_TPL_DEPENDENCIES
    "boost"
    )

# Make sure needed moris libraries are built
#include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)

# Include third party libraries indirectly needed by ALG
list(APPEND ALG_TPL_DEPENDENCIES
    #${LINALG_TPL_DEPENDENCIES}
    )

list(REMOVE_DUPLICATES ALG_TPL_DEPENDENCIES)    