#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Parameters ----------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)
add_subdirectory(ENM)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_PRM)
    add_subdirectory(test)
endif()
