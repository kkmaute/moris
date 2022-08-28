#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Communications Tests ----------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Communication_Tools.cpp)

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
    ${COM}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(COM_TEST ${COM}-test)

# Create executable
add_executable(${COM_TEST} ${TEST_SOURCES})
target_link_libraries(${COM_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${COM_TEST} PROPERTIES OUTPUT_NAME ${COM_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(COM_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${COM_TEST_PROCS})
		set(COM_PARALLEL_TEST ${COM_TEST}-${PROCS}-procs)
		add_test(NAME ${COM_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${COM_TEST}>)
	endforeach()
endif()


