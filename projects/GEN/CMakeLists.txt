# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Geometry Engine

# Build the module library.
add_subdirectory(src)

#if(BUILD_SDF)
#	add_subdirectory(SDF)
#endif()

# Build tests when this module is selected.
if(MORIS_USE_TESTS AND BUILD_GEN)
    add_subdirectory(test)
endif()
