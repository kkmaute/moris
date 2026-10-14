# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Parameters Tests

# Sources.
set(TEST_SOURCES test_main.cpp UT_FEM_Parameters.cpp)

# List snippet includes
set(SNIPPET_INCLUDES)

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${PRM}-lib)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(PRM_TEST ${PRM}-test)

# Define the executable and its dependencies.
add_executable(${PRM_TEST} ${TEST_SOURCES})
target_link_libraries(${PRM_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${PRM_TEST} PROPERTIES OUTPUT_NAME ${PRM_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${PRM_TEST} ${MORIS_TEST_PROCS})
