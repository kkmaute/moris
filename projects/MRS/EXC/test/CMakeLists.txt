# Exceptions Tests --------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_invalid_size.cpp
    cl_out_of_range.cpp
    cl_tic_toc_error.cpp )

# List test includes
set(TEST_INCLUDES
    ${COM} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(EXC_TEST ${EXC}-test)

# Create executable
add_executable(${EXC_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${EXC_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${EXC_TEST}
    COMMAND ${EXC_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(EXC_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${EXC_TEST_PROCS})
		set(EXC_PARALLEL_TEST ${EXC_TEST}-${PROCS}-procs)
		add_test(NAME ${EXC_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${EXC_TEST}${EXE_EXT}>)
	endforeach()
endif()
