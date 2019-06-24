# Tools Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Debug.cpp
    cl_Geometry.cpp)

# List test includes
set(TEST_INCLUDES
    #${LINALG}
    #${COM}
    #${ALG}
    #${TOL}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
    ${TOL}-lib
    #${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    #${MORIS_ACML_LAPACK_MKL_LIBS}
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_ARPACK_LIBS}
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(TOL_TEST ${TOL}-test)

# Create executable
add_executable(${TOL_TEST} ${TEST_SOURCES})
target_link_libraries(${TOL_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${TOL_TEST} PROPERTIES OUTPUT_NAME ${TOL}-test${EXE_EXT})

# Mark as test
add_test(NAME ${TOL_TEST}
    COMMAND $<TARGET_FILE:${TOL_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(TOL_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${TOL_TEST_PROCS})
		set(TOL_PARALLEL_TEST ${TOL_TEST}-${PROCS}-procs)
		add_test(NAME ${TOL_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${TOL_TEST}>)
	endforeach()
endif()
