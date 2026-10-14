# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Distributed Linear Algebra Tests

include_directories(DLA_Test_Proxy)

# Sources.
set(TEST_SOURCES
    test_main.cpp
    cl_Dist_Vector_Test.cpp #> VectorType error
    cl_Linear_Solver_Test.cpp
    cl_Map_Test.cpp
    UT_convert_epetra_to_matrix.cpp
    cl_Sparse_Matrix_Test.cpp
    cl_DLA_Multigrid_Test.cpp
    DLA_Test_Proxy/cl_Solver_Interface_Proxy.cpp
    map.mtx
    matrix.mtx
    vector.mtx
)

if(${MORIS_HAVE_SLEPC})
    list(APPEND TEST_SOURCES cl_Eigen_Solver_Slepc.cpp)
endif()

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${ALG}-lib ${DLA}-lib ${COM}-lib ${MTK}-lib ${PRM}-lib)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(DLA_TEST ${DLA}-test)

# Define the executable and its dependencies.
add_executable(${DLA_TEST} ${TEST_SOURCES})
target_link_libraries(${DLA_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${DLA_TEST} PROPERTIES OUTPUT_NAME ${DLA_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${DLA_TEST} ${MORIS_TEST_PROCS})
