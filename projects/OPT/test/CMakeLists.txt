# Optimization Tests ------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    opt_test.cpp )

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${OPT} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${OPT}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_GCMMA_LIBS}
    ${MORIS_SNOPT_LIBS} )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(OPT_TEST ${OPT}-test)

# Create executable
add_executable(${OPT_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${OPT_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${OPT_TEST}
    COMMAND ${OPT_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(OPT_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${OPT_TEST_PROCS})
		set(OPT_PARALLEL_TEST ${OPT_TEST}-${PROCS}-procs)
		add_test(NAME ${OPT_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${OPT_TEST}${EXE_EXT}>)
	endforeach()
endif()

