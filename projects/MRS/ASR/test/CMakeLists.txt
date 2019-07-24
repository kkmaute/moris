# Assert Tests ------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    assert.cpp)

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(ASR_TEST ${ASR}-test)

# Create executable
add_executable(${ASR_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${ASR_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${ASR_TEST}
    COMMAND ${ASR_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(ASR_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${ASR_TEST_PROCS})
		set(ASR_PARALLEL_TEST ${ASR_TEST}-${PROCS}-procs)
		add_test(NAME ${ASR_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${ASR_TEST}${EXE_EXT}>)
	endforeach()
endif()
