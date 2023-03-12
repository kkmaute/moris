#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Geometry Engine ---------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directory
add_subdirectory(src)

#if(BUILD_SDF)
#	add_subdirectory(SDF)
#endif()

# Add test directory
if(MORIS_USE_TESTS AND BUILD_GEN)
    add_subdirectory(test)
endif()

