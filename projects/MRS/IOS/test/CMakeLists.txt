# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# IOS Tests

# Sources.
set(TEST_SOURCES test_main.cpp cl_Logger.cpp fn_to_stdio.cpp UT_IOS_Parsing_Tool.cpp
                 UT_IOS_File_To_Array.cpp
)

# List additional includes
include_directories(${MORIS_DIR}/snippets/ios)
include_directories(${CMAKE_BINARY_DIR}/generated/)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name for general IOS test
set(IOS_TEST ${IOS}-test)

# Define the executable and its dependencies.
add_executable(${IOS_TEST} ${TEST_SOURCES})
target_link_libraries(${IOS_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${IOS_TEST} PROPERTIES OUTPUT_NAME ${IOS_TEST}${EXE_EXT})

# List source files for Query
set(QUERY_SOURCES query_main.cpp cl_query_test.cpp)

# Set test name for Query
set(IOS_QUERY ${IOS}-query)

# Create executable for Query
add_executable(${IOS_QUERY} ${QUERY_SOURCES})
target_link_libraries(${IOS_QUERY} ${TEST_DEPENDENCIES})
set_target_properties(${IOS_QUERY} PROPERTIES OUTPUT_NAME ${IOS_QUERY}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${IOS_TEST} ${MORIS_TEST_PROCS})
