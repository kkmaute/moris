# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Chronos Tests

# Sources.
set(TEST_SOURCES test_main.cpp cl_Stopwatch.cpp cl_Waitbar.cpp)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${CHR}-lib ${MORIS_BASE_LIBS})

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(CHR_TEST ${CHR}-test)

# Define the executable and its dependencies.
add_executable(${CHR_TEST} ${TEST_SOURCES})
target_link_libraries(${CHR_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${CHR_TEST} PROPERTIES OUTPUT_NAME ${CHR_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${CHR_TEST} ${MORIS_TEST_PROCS})
