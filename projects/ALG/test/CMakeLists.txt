# Algorithms Tests --------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    fn_equal_to.cpp
    fn_num_digits.cpp)

# List test includes
set(TEST_INCLUDES
    #${LINALG}
    #${COM}
    #${ALG}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List snippet includes
set(SNIPPET_INCLUDES
    algorithms )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
	${ALG}-lib
    ${COM}-lib )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(ALG_TEST ${ALG}-test)

# Create executable
add_executable(${ALG_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${ALG_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${ALG_TEST}
    COMMAND ${ALG_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(ALG_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${ALG_TEST_PROCS})
		set(ALG_PARALLEL_TEST ${ALG_TEST}-${PROCS}-procs)
		add_test(NAME ${ALG_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${ALG_TEST}${EXE_EXT}>)
	endforeach()
endif()

