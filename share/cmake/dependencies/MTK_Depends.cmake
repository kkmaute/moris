#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# MTK Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if MTK has already been included
if(DEFINED MTK_CONFIGURED_ONCE)
    return()
endif()

set(MTK_CONFIGURED_ONCE "YES")

# Add MTK to the header directory list
list(APPEND MORIS_SOURCE_DIRS ${MTK})

# Third party libraries needed by MTK-lib
set(MTK_TPL_DEPENDENCIES
    "trilinos"
    )

# Make sure needed moris libraries are built
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)

# for test
include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)

