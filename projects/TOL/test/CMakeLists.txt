# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Tools Tests

# Sources.
set(TEST_SOURCES test_main.cpp cl_Debug.cpp cl_Geometry.cpp)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${TOL}-lib ${MORIS_BASE_LIBS})

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(TOL_TEST ${TOL}-test)

# Define the executable and its dependencies.
add_executable(${TOL_TEST} ${TEST_SOURCES})
target_link_libraries(${TOL_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${TOL_TEST} PROPERTIES OUTPUT_NAME ${TOL}-test${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${TOL_TEST} ${MORIS_TEST_PROCS})
