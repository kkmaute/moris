#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Exceptions Tests --------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_invalid_size.cpp
    cl_out_of_range.cpp
    cl_tic_toc_error.cpp )

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(EXC_TEST ${EXC}-test)

# Create executable
add_executable(${EXC_TEST} ${TEST_SOURCES})
target_link_libraries(${EXC_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${EXC_TEST} PROPERTIES OUTPUT_NAME ${EXC_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(EXC_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${EXC_TEST_PROCS})
		set(EXC_PARALLEL_TEST ${EXC_TEST}-${PROCS}-procs)
		add_test(NAME ${EXC_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${EXC_TEST}>)
	endforeach()
endif()

