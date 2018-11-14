# Geometry Engine Dependencies --------------------------------------------
# -------------------------------------------------------------------------

# Check if GEN has already been included
if(DEFINED GEN_CONFIGURED_ONCE)
    return()
endif()

set(GEN_CONFIGURED_ONCE "YES")

# Add GEN to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${GEN})

# Include libraries needed by GEN
set(GEN_TPL_DEPENDENCIES "")

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake) #> headers
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake) #> headers

# Include third party libraries indirectly needed by GEN
list(APPEND GEN_TPL_DEPENDENCIES
    ${MTK_TPL_DEPENDENCIES}
    ${LINALG_TPL_DEPENDENCIES}
    ${ARMADILLO_EIGEN}
    )
