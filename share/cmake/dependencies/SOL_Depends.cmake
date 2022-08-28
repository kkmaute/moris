#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# SOL Dependencies --------------------------------------------------------
# -------------------------------------------------------------------------

# Check if SOL has already been included
if(DEFINED SOL_CONFIGURED_ONCE)
    return()
endif()

set(SOL_CONFIGURED_ONCE "YES")

# SOL is only a wrapper around these moris packages (Aug. 20, 2018)
include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/SOL_CORE_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/TSA_Depends.cmake)

