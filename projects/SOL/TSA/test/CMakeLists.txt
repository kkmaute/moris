#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Time Solver Algebra Tests ----------------------------------------
# -------------------------------------------------------------------------

    include_directories(TSA_Test_Proxy)
    
# List source files
set(TEST_SOURCES
    cl_TSA_Time_Solver_Tests.cpp
    cl_TSA_Time_Solver_Staggered_Tests.cpp
    TSA_Test_Proxy/cl_TSA_Solver_Interface_Proxy.cpp
    TSA_Test_Proxy/cl_TSA_Solver_Interface_Proxy2.cpp
    test_main.cpp
    )

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
	${ALG}-lib
	${NLA}-lib
	${MSI}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(TSA_TEST ${TSA}-test)

# Create executable
add_executable(${TSA_TEST} ${TEST_SOURCES})
target_link_libraries(${TSA_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${TSA_TEST} PROPERTIES OUTPUT_NAME ${TSA_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(TSA_TEST_PROCS ${MORIS_TEST_PROCS})
    foreach(PROCS ${TSA_TEST_PROCS})
        set(TSA_PARALLEL_TEST ${TSA_TEST}-${PROCS}-procs)
        add_test(NAME ${TSA_PARALLEL_TEST} 
            COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${TSA_TEST}>)
    endforeach()
endif()


