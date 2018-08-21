# FEM Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if FEM has already been included
if(DEFINED FEM_CONFIGURED_ONCE)
    return()
endif()

set(FEM_CONFIGURED_ONCE "YES")

# FEM is only a wrapper around these moris packages (Aug. 20, 2018)
include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
