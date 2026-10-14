# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Time Solver Algebra Tests

include_directories(TSA_Test_Proxy)

# Sources.
set(TEST_SOURCES
    cl_TSA_Time_Solver_Tests.cpp cl_TSA_Time_Solver_Staggered_Tests.cpp
    TSA_Test_Proxy/cl_TSA_Solver_Interface_Proxy.cpp
    TSA_Test_Proxy/cl_TSA_Solver_Interface_Proxy2.cpp test_main.cpp
)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${ALG}-lib ${NLA}-lib ${MSI}-lib)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(TSA_TEST ${TSA}-test)

# Define the executable and its dependencies.
add_executable(${TSA_TEST} ${TEST_SOURCES})
target_link_libraries(${TSA_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${TSA_TEST} PROPERTIES OUTPUT_NAME ${TSA_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${TSA_TEST} ${MORIS_TEST_PROCS})
