# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Exceptions Tests

# Sources.
set(TEST_SOURCES test_main.cpp cl_invalid_size.cpp cl_out_of_range.cpp cl_tic_toc_error.cpp)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(EXC_TEST ${EXC}-test)

# Define the executable and its dependencies.
add_executable(${EXC_TEST} ${TEST_SOURCES})
target_link_libraries(${EXC_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${EXC_TEST} PROPERTIES OUTPUT_NAME ${EXC_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${EXC_TEST} ${MORIS_TEST_PROCS})
