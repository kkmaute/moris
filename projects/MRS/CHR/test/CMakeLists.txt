# Chronos Tests -----------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Stopwatch.cpp
    cl_Waitbar.cpp )

# List test includes
set(TEST_INCLUDES
    ${COM} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additonal includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_BOOST_LIBS})

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(CHR_TEST ${CHR}-test)

# Create executable
add_executable(${CHR_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${CHR_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${CHR_TEST}
    COMMAND ${CHR_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(CHR_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${CHR_TEST_PROCS})
		set(CHR_PARALLEL_TEST ${CHR_TEST}-${PROCS}-procs)
		add_test(NAME ${CHR_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${CHR_TEST}${EXE_EXT}>)
	endforeach()
endif()
