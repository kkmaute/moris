#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Distributed Linear Algebra Dependencies ---------------------------------
# -------------------------------------------------------------------------

# Check if DLA has already been included
if(DEFINED DLA_CONFIGURED_ONCE)
    return()
endif()

set(DLA_CONFIGURED_ONCE "YES")

# Add DLA to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${SOL}/${DLA})

# Include libraries needed by DLA
set(DLA_TPL_DEPENDENCIES
    "mpi"
    )

# Include petsc if used
if(${MORIS_HAVE_PETSC})
    list(APPEND DLA_TPL_DEPENDENCIES "PETSc")
    if(${MORIS_HAVE_SLEPC})
        list(APPEND DLA_TPL_DEPENDENCIES "SLEPc")
    endif()
endif()

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)

# just needed for tests
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

