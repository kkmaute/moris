#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Containers Tests --------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Array.cpp
    cl_Bi_Map.cpp
	ut_Vector.cpp
    cl_Dist_Map.cpp
    cl_Map.cpp
    cl_Param_List.cpp
    cl_Tuple.cpp )

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
	${CNT}-lib
    ${COM}-lib
    ${ALG}-lib
    ${MORIS_BASE_LIBS})

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(CNT_TEST ${CNT}-test)

# Create executable
add_executable(${CNT_TEST} ${TEST_SOURCES})
target_link_libraries(${CNT_TEST} PUBLIC ${TEST_DEPENDENCIES})
set_target_properties(${CNT_TEST} PROPERTIES OUTPUT_NAME ${CNT_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(CNT_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${CNT_TEST_PROCS})
		set(CNT_PARALLEL_TEST ${CNT_TEST}-${PROCS}-procs)
		add_test(NAME ${CNT_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${CNT_TEST}>)
	endforeach()
endif()

