#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Optimization Dependencies -----------------------------------------------
# -------------------------------------------------------------------------

# Check if OPT has already been included
if(DEFINED OPT_CONFIGURED_ONCE)
    return()
endif()

set(OPT_CONFIGURED_ONCE "YES")

# Add OPT to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${OPT})

# Include libraries needed by OPT
set(OPT_TPL_DEPENDENCIES )

# Include gcmma if used
if(${MORIS_HAVE_GCMMA})
    list(APPEND OPT_TPL_DEPENDENCIES "gcmma")
endif()
    
# Include snopt if used
if(${MORIS_HAVE_SNOPT})
    list(APPEND OPT_TPL_DEPENDENCIES "snopt")
endif()
    
# Include lbfgsb if used
if(${MORIS_HAVE_LBFGS})
    list(APPEND OPT_TPL_DEPENDENCIES "lbfgsb")
endif()
    

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)

# for test
include(${MORIS_DEPENDS_DIR}/PRM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

