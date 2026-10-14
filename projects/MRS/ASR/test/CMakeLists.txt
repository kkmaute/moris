# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Assert Tests

# Sources.
set(TEST_SOURCES test_main.cpp assert.cpp)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(ASR_TEST ${ASR}-test)

# Define the executable and its dependencies.
add_executable(${ASR_TEST} ${TEST_SOURCES})
target_link_libraries(${ASR_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${ASR_TEST} PROPERTIES OUTPUT_NAME ${ASR_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${ASR_TEST} ${MORIS_TEST_PROCS})
