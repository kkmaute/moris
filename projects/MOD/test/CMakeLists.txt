# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Model Tests

# Sources.
set(TEST_SOURCES test_main.cpp cl_Bspline.cpp EgnManager.cpp)

# Test dependencies.
set(TEST_DEPENDENCIES)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MOD_TEST ${MOD}-test)

# Define the executable and its dependencies.
add_executable(${MOD_TEST} ${TEST_SOURCES})
target_link_libraries(${MOD_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${MOD_TEST} PROPERTIES OUTPUT_NAME ${MTK_TEST}${EXE_EXT})

# Mark as test

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${MOD_TEST} ${MORIS_TEST_PROCS})
