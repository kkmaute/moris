# Distributed Linear Algebra Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if MAIN has already been included
if(DEFINED MAIN_CONFIGURED_ONCE)
    return()
endif()

set(MAIN_CONFIGURED_ONCE "YES")

# Add MAIN to the source directory list
list(APPEND MORIS_SRC_DIRS ${MAIN})

# Include all libraries for MAIN
set(MAIN_TPL_DEPENDENCIES
    "")

include(share/cmake/MatthewCMake/ALG_Depends.cmake)
include(share/cmake/MatthewCMake/COM_Depends.cmake)
include(share/cmake/MatthewCMake/DLA_Depends.cmake)
include(share/cmake/MatthewCMake/GEN_Depends.cmake)
include(share/cmake/MatthewCMake/HMR_Depends.cmake)
include(share/cmake/MatthewCMake/INT_Depends.cmake)
include(share/cmake/MatthewCMake/IOS_Depends.cmake)
include(share/cmake/MatthewCMake/LNA_Depends.cmake)
include(share/cmake/MatthewCMake/MOD_Depends.cmake)
include(share/cmake/MatthewCMake/MSI_Depends.cmake)
include(share/cmake/MatthewCMake/MTK_Depends.cmake)
include(share/cmake/MatthewCMake/OPT_Depends.cmake)
include(share/cmake/MatthewCMake/STK_Depends.cmake)
include(share/cmake/MatthewCMake/TOL_Depends.cmake)
# include(share/cmake/MatthewCMake/XTK_Depends.cmake)

list(APPEND MAIN_TPL_DEPENDENCIES
    ${LNA_TPL_DEPENDENCIES} )
