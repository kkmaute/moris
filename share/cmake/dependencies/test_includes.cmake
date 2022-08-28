#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Unit Test Dependencies --------------------------------------------------
# -------------------------------------------------------------------------

# Make sure to build dependencies
include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)

# Create linkable test target
add_library(test-libs INTERFACE)
target_link_libraries(test-libs INTERFACE ${LINALG}-lib ${COM}-lib)

