#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# IOS ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES IOS_TPL_DEPENDENCIES)
#foreach(IOS_TPL ${IOS_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${IOS_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_IOS)
    add_subdirectory(test)
endif()

