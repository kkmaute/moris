#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Interpolation and Integration -------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

if( BUILD_TUTORIALS AND BUILD_MDL )
    add_subdirectory(tutorials)
endif()

# Add test directory
if(MORIS_USE_TESTS AND BUILD_MDL)
    add_subdirectory(test)
endif()
