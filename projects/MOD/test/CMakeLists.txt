#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Model Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Bspline.cpp
    EgnManager.cpp )

# List test dependencies
set(TEST_DEPENDENCIES
    )
    
# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MOD_TEST ${MOD}-test)

# Create executable
add_executable(${MOD_TEST} ${TEST_SOURCES})
target_link_libraries(${MOD_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${MOD_TEST} PROPERTIES OUTPUT_NAME ${MTK_TEST}${EXE_EXT})

# Mark as test
add_test(NAME ${MOD_TEST}
    COMMAND ${MORIS_EXECUTE_COMMAND} $<TARGET_FILE:${MOD_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(MOD_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${MOD_TEST_PROCS})
		set(MOD_PARALLEL_TEST ${MOD_TEST}-${PROCS}-procs)
		add_test(NAME ${MOD_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${MOD_TEST}>)
	endforeach()
endif()


