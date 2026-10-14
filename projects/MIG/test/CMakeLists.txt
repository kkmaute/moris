# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# MIG Tests

# Sources.
set(TEST_SOURCES test_main.cpp UT_MIG_Coords.cpp)
# )

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${MIG}-lib ${ALG}-lib ${LINALG}-lib ${MORIS_BASE_LIBS})

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MIG_TEST ${MIG}-test)

# Define the executable and its dependencies.
add_executable(${MIG_TEST} ${TEST_SOURCES})
target_link_libraries(${MIG_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${MIG_TEST} PROPERTIES OUTPUT_NAME ${MIG_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${MIG_TEST} ${MORIS_TEST_PROCS})
