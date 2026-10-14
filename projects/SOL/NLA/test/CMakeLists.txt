# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Non Linear Algebra Tests
include_directories(${MORIS_PACKAGE_DIR}/SOL/NLA/test/NLA_Test_Proxy)
include_directories(${MORIS_PACKAGE_DIR}/SOL/TSA/test/TSA_Test_Proxy)

# Sources.
set(TEST_SOURCES
    test_main.cpp
    cl_NLA_Newton_Solver_Test.cpp
    cl_NLA_NonlinearDatabase.cpp
    NLA_Test_Proxy/cl_NLA_Solver_Interface_Proxy.cpp
    NLA_Test_Proxy/cl_NLA_Solver_Interface_Proxy2.cpp
    ${MORIS_PACKAGE_DIR}/SOL/TSA/test/TSA_Test_Proxy/cl_TSA_Solver_Interface_Proxy2.cpp
)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${ALG}-lib ${MSI}-lib ${DLA}-lib ${NLA}-lib)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(NLA_TEST ${NLA}-test)

# Define the executable and its dependencies.
add_executable(${NLA_TEST} ${TEST_SOURCES})
target_link_libraries(${NLA_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${NLA_TEST} PROPERTIES OUTPUT_NAME ${NLA_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${NLA_TEST} ${MORIS_TEST_PROCS})
