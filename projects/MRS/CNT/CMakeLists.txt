#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Containers --------------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_CNT)
    add_subdirectory(test)
endif()

