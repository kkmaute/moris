# XTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if XTK has already been included
if(DEFINED XTK_CONFIGURED_ONCE)
    return()
endif()

set(XTK_CONFIGURED_ONCE "YES")

# Add XTK to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${XTK})



# Make sure needed moris libraries are built
#include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/FEM_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)

include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/TOL_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

# Moris packages included for exe's
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)

if(BUILD_HMR)
	include(${MORIS_DEPENDS_DIR}/GEN_Depends.cmake)
	include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
	include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
	include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
	include(${MORIS_DEPENDS_DIR}/MAP_Depends.cmake)
endif()

# needs some tpls
set(XTK_LIB_TPL_DEPENDENCIES
    #${ARMADILLO_EIGEN}
    #"superlu" # needed for exe linking order
    #"trilinos"
    #"arpack" # needed for exe linking order
    )

set(XTK_EXE_TPL_DEPENDENCIES
	#"superlu"
	#"trilinos"
	)

