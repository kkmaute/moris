# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Optimization Tests

# Sources.
set(TEST_SOURCES test_main.cpp opt_test.cpp)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${OPT}-lib ${COM}-lib ${PRM}-lib ${MORIS_BASE_LIBS})

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(OPT_TEST ${OPT}-test)

# Define the executable and its dependencies.
add_executable(${OPT_TEST} ${TEST_SOURCES})
target_link_libraries(${OPT_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${OPT_TEST} PROPERTIES OUTPUT_NAME ${OPT_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${OPT_TEST} ${MORIS_TEST_PROCS})
