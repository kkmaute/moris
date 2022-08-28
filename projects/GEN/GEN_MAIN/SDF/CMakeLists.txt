#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Interpolation and Integration -------------------------------------------
# -------------------------------------------------------------------------

# Third party library dependencies
#list(REMOVE_DUPLICATES SDF_TPL_DEPENDENCIES) <-- Matt, please find out why this does not work on facon.
#foreach(SDF_TPL ${SDF_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${SDF_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_SDF)
    add_subdirectory(test)
endif()
