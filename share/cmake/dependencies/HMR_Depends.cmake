#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# HMR Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if HMR has already been included
if(DEFINED HMR_CONFIGURED_ONCE)
    return()
endif()

set(HMR_CONFIGURED_ONCE "YES")

# Add HMR to the source directory list
list(APPEND MORIS_SOURCE_DIRS ${HMR})

# Include libraries needed by HMR
set(HMR_TPL_DEPENDENCIES
    ""
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/GEN_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

# added as temp fix for hmr exe, test, and tutorials
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)

# added for tutorials
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/SDF_Depends.cmake)
