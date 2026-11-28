# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

include(${CMAKE_CURRENT_SOURCE_DIR}/example_test_helpers.cmake)

add_subdirectory(thermal)
add_subdirectory(structure)
add_subdirectory(fluid)
add_subdirectory(optimization)
add_subdirectory(mesh_generation)
