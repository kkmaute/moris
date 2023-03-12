#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# MTK ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

# Add test directory
if(MORIS_USE_TESTS AND BUILD_MTK)
     add_subdirectory(test) #> no tests
endif()

# Add test directory
if(BUILD_TUTORIALS AND BUILD_MTK)
     add_subdirectory(tutorial) #> no tests
endif()
