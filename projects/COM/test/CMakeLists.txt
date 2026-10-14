# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Communications Tests

# Sources.
set(TEST_SOURCES test_main.cpp cl_Communication_Tools.cpp)

# List snippet includes
set(SNIPPET_INCLUDES)

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${COM}-lib)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(COM_TEST ${COM}-test)

# Define the executable and its dependencies.
add_executable(${COM_TEST} ${TEST_SOURCES})
target_link_libraries(${COM_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${COM_TEST} PROPERTIES OUTPUT_NAME ${COM_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${COM_TEST} ${MORIS_TEST_PROCS})
