#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Optimization Tests ------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    opt_test.cpp )

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${OPT}-lib
    ${COM}-lib
    ${PRM}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(OPT_TEST ${OPT}-test)

# Create executable
add_executable(${OPT_TEST} ${TEST_SOURCES})
target_link_libraries(${OPT_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${OPT_TEST} PROPERTIES OUTPUT_NAME ${OPT_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(OPT_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${OPT_TEST_PROCS})
		set(OPT_PARALLEL_TEST ${OPT_TEST}-${PROCS}-procs)
		add_test(NAME ${OPT_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${OPT_TEST}>)
	endforeach()
endif()


