#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# LINALG --------------------------------------------------------------------
# -------------------------------------------------------------------------
# Third party library dependencies
#list(REMOVE_DUPLICATES LINALG_TPL_DEPENDENCIES)
#foreach(LINALG_TPL ${LINALG_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${LINALG_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_LINALG)
    add_subdirectory(test)
endif()

