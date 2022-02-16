# MIG Tests ----------------------------------------
# -------------------------------------------------------------------------

    
# List source files
set(TEST_SOURCES
    test_main.cpp
    UT_MIG_Coords.cpp )
   # )

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
	${MIG}-lib
    ${ALG}-lib
	${LINALG}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MIG_TEST ${MIG}-test)

# Create executable
add_executable(${MIG_TEST} ${TEST_SOURCES})
target_link_libraries(${MIG_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${MIG_TEST} PROPERTIES OUTPUT_NAME ${MIG_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(MIG_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${MIG_TEST_PROCS})
		set(MIG_PARALLEL_TEST ${MIG_TEST}-${PROCS}-procs)
		add_test(NAME ${MIG_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${MIG_TEST}>)
	endforeach()
endif()

