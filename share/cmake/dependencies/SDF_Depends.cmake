# SDF Dependencies ------------------------------
# -------------------------------------------------------------------------

# Check if SDF has already been included
if(DEFINED SDF_CONFIGURED_ONCE)
    return()
endif()

set(SDF_CONFIGURED_ONCE "YES")

# Add SDF to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${GEN}/${SDF})

# Third party libraries needed by SDF library
# PETSc and Trilinos; add later
# include(${MORIS_TPL_DIR}/PETSc.cmake)
set(SDF_LIB_TPL_DEPENDENCIES
    "superlu" # needed for test linking order
    "arpack" # needed for test linking order
    )

# Additional third party libraries needed by SDF executables
set(SDF_EXE_TPL_DPENDENCIES
	""
	)

# Moris packages needed by SDF
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
#include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)

#list(APPEND SDF_TPL_DEPENDENCIES
#	"boost"
#	${CON_TPL_DEPENDENCIES}
#     ${LINALG_TPL_DEPENDENCIES} )