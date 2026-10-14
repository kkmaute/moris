# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Algorithms Tests

# Sources.
set(TEST_SOURCES test_main.cpp fn_equal_to.cpp fn_num_digits.cpp)

# List snippet includes
set(SNIPPET_INCLUDES algorithms)

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${ALG}-lib)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(ALG_TEST ${ALG}-test)

# Define the executable and its dependencies.
add_executable(${ALG_TEST} ${TEST_SOURCES})
target_link_libraries(${ALG_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${ALG_TEST} PROPERTIES OUTPUT_NAME ${ALG_TEST}${EXE_EXT})

# Add tests
moris_register_unit_tests(${ALG_TEST} ${MORIS_TEST_PROCS})
