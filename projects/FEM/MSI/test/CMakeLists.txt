# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Model Solver Interface Tests

include_directories(MSI_Test_Proxy)

# Sources.
set(TEST_SOURCES
    test_main.cpp
    cl_MSI_Dof_Manager_Test.cpp
    cl_MSI_Equation_Object_Test.cpp
    cl_MSI_Pdof_Host_Test.cpp
    cl_MSI_Model_Solver_Interface_Test.cpp
    cl_MSI_Multigrid.cpp
    UT_Equation_Object_Pdv.cpp
    #UT_MSI_Element_Time_Sideset.cpp
    MSI_Test_Proxy/cl_MSI_Design_Variable_Interface_Proxy.cpp
    MSI_Test_Proxy/cl_MSI_Solver_Interface_Proxy.cpp
    cl_MSI_SpaceTime_Test.cpp
    UT_Sparsity_Pattern.cpp
)

# Test dependencies.
set(TEST_DEPENDENCIES
    test-libs
    ${ALG}-lib
    ${NLA}-lib
    ${MSI}-lib
    ${PRM}-lib
    ${INT}-lib
    ${MORIS_BASE_LIBS}
)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MSI_TEST ${MSI}-test)

# Define the executable and its dependencies.
add_executable(${MSI_TEST} ${TEST_SOURCES})
target_link_libraries(${MSI_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${MSI_TEST} PROPERTIES OUTPUT_NAME ${MSI_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${MSI_TEST} ${MORIS_TEST_PROCS})
