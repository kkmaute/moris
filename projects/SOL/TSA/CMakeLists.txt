#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Time Solver Algebra ----------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_TSA)
    add_subdirectory(test)
endif()

# Not sure yet what else should be here

