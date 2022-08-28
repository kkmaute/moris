#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Third party library dependencies
#list(REMOVE_DUPLICATES MAP_TPL_DEPENDENCIES)
#foreach(MAP_TPL ${MAP_TPL_DEPENDENCIES})
#    include(${MORIS_TPL_DIR}/${MAP_TPL}.cmake)
#endforeach()

# Add source directory
add_subdirectory(src)
