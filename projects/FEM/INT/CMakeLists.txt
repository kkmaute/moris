#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Interpolation and Integration -------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES INT_TPL_DEPENDENCIES)
#foreach(INT_TPL ${INT_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${INT_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_INT)
    add_subdirectory(test)
endif()

