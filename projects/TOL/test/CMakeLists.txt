#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Tools Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Debug.cpp
    cl_Geometry.cpp)

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
    ${TOL}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(TOL_TEST ${TOL}-test)

# Create executable
add_executable(${TOL_TEST} ${TEST_SOURCES})
target_link_libraries(${TOL_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${TOL_TEST} PROPERTIES OUTPUT_NAME ${TOL}-test${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(TOL_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${TOL_TEST_PROCS})
		set(TOL_PARALLEL_TEST ${TOL_TEST}-${PROCS}-procs)
		add_test(NAME ${TOL_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${TOL_TEST}>)
	endforeach()
endif()

