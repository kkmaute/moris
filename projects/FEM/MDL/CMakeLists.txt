# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Interpolation and Integration

# Build the module library.
add_subdirectory(src)

if(BUILD_TUTORIALS AND BUILD_MDL)
    add_subdirectory(tutorials)
endif()

# Build tests when this module is selected.
if(MORIS_USE_TESTS AND BUILD_MDL)
    add_subdirectory(test)
endif()
