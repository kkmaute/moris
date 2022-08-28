#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Chronos Tests -----------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Stopwatch.cpp
    cl_Waitbar.cpp )

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
	${CHR}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(CHR_TEST ${CHR}-test)

# Create executable
add_executable(${CHR_TEST} ${TEST_SOURCES})
target_link_libraries(${CHR_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${CHR_TEST} PROPERTIES OUTPUT_NAME ${CHR_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(CHR_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${CHR_TEST_PROCS})
		set(CHR_PARALLEL_TEST ${CHR_TEST}-${PROCS}-procs)
		add_test(NAME ${CHR_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${CHR_TEST}>)
	endforeach()
endif()

