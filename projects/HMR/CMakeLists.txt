#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# HMR ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES HMR_TPL_DEPENDENCIES)
#foreach(HMR_TPL ${HMR_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${HMR_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

if( BUILD_TUTORIALS AND BUILD_HMR )
    add_subdirectory(tutorials)
endif()

# Add test directory
if(MORIS_USE_TESTS AND BUILD_HMR)
    add_subdirectory(test)
endif()

