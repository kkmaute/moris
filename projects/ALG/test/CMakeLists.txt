#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Algorithms Tests --------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    fn_equal_to.cpp
    fn_num_digits.cpp)

# List snippet includes
set(SNIPPET_INCLUDES
    algorithms )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
	${ALG}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(ALG_TEST ${ALG}-test)

# Create executable
add_executable(${ALG_TEST} ${TEST_SOURCES})
target_link_libraries(${ALG_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${ALG_TEST} PROPERTIES OUTPUT_NAME ${ALG_TEST}${EXE_EXT})

# Add tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(ALG_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${ALG_TEST_PROCS})
		set(ALG_PARALLEL_TEST ${ALG_TEST}-${PROCS}-procs)
		add_test(NAME ${ALG_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${ALG_TEST}>)
	endforeach()
endif()

