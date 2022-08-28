#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Parameters Tests ----------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES 
    test_main.cpp
    UT_FEM_Parameters.cpp )

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${PRM}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(PRM_TEST ${PRM}-test)

# Create executable
add_executable(${PRM_TEST} ${TEST_SOURCES})
target_link_libraries(${PRM_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${PRM_TEST} PROPERTIES OUTPUT_NAME ${PRM_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(PRM_TEST_PROCS ${MORIS_TEST_PROCS})
    foreach(PROCS ${PRM_TEST_PROCS})
        set(PRM_PARALLEL_TEST ${PRM_TEST}-${PROCS}-procs)
        add_test(NAME ${PRM_PARALLEL_TEST} 
            COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${PRM_TEST}>)
    endforeach()
endif()


