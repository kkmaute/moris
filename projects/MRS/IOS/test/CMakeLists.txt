#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# IOS Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Logger.cpp
    fn_to_stdio.cpp
    UT_IOS_Parsing_Tool.cpp
    UT_IOS_File_To_Array.cpp )

# List additional includes
include_directories(${MORIS_DIR}/snippets/ios)
include_directories(${CMAKE_BINARY_DIR}/generated/)

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    )


# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name for general IOS test
set(IOS_TEST ${IOS}-test)

# Create executable
add_executable(${IOS_TEST} ${TEST_SOURCES})
target_link_libraries(${IOS_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${IOS_TEST} PROPERTIES OUTPUT_NAME ${IOS_TEST}${EXE_EXT})


# List source files for Query
set(QUERY_SOURCES
    query_main.cpp
    cl_query_test.cpp )

# Set test name for Query
set(IOS_QUERY ${IOS}-query)

# Create executable for Query
add_executable(${IOS_QUERY} ${QUERY_SOURCES})
target_link_libraries(${IOS_QUERY} ${TEST_DEPENDENCIES})
set_target_properties(${IOS_QUERY} PROPERTIES OUTPUT_NAME ${IOS_QUERY}${EXE_EXT})


# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(IOS_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${IOS_TEST_PROCS})
		set(IOS_PARALLEL_TEST ${IOS_TEST}-${PROCS}-procs)
		add_test(NAME ${IOS_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${IOS_TEST}>)
	endforeach()
endif()

