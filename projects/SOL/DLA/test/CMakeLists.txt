#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Distributed Linear Algebra Tests ----------------------------------------
# -------------------------------------------------------------------------

    include_directories(DLA_Test_Proxy)
    
# List source files
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
    vector.mtx )

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
	${ALG}-lib
    ${DLA}-lib
    ${COM}-lib
    ${MTK}-lib
    ${PRM}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(DLA_TEST ${DLA}-test)

# Create executable
add_executable(${DLA_TEST} ${TEST_SOURCES})
target_link_libraries(${DLA_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${DLA_TEST} PROPERTIES OUTPUT_NAME ${DLA_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(DLA_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${DLA_TEST_PROCS})
		set(DLA_PARALLEL_TEST ${DLA_TEST}-${PROCS}-procs)
		add_test(NAME ${DLA_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${DLA_TEST}>)
	endforeach()
endif()


