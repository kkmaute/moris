#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# GEN_CORE Geometry Engine ------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

## Add test directory # maybe add core tests when in debug
#if(MORIS_USE_TESTS AND BUILD_SOL_CORE)
#    add_subdirectory(test)
#endif()

