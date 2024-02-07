#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# SDF Tests ---------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    ut_SDF_Line.cpp
    ut_SDF_Triangle.cpp 
    ut_SDF_Raycast.cpp
    ut_SDF_Core.cpp
    ut_SDF_Object.cpp
    )

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
    ${SDF}-lib
    ${MTK}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(SDF_TEST ${SDF}-test)

# Create executable
add_executable(${SDF_TEST} ${TEST_SOURCES})
target_link_libraries(${SDF_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${SDF_TEST} PROPERTIES OUTPUT_NAME ${SDF_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(SDF_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${SDF_TEST_PROCS})
		set(SDF_PARALLEL_TEST ${SDF_TEST}-${PROCS}-procs)
		add_test(NAME ${SDF_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${SDF_TEST}>)
	endforeach()
endif()


