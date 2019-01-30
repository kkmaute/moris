# Non Linear Algebra Tests ----------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_NLA_Newton_Solver_Test.cpp
    cl_NLA_NonlinearDatabase.cpp)

# List test includes
set(TEST_INCLUDES
    ${SOL}/${DLA}
    ${LINALG}
    ${COM}
    ${SOL}/${NLA}
    ${ALG}
    ${FEM}/${MSI} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${DLA}-lib
    ${COM}-lib
    ${NLA}-lib
    ${MSI}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(NLA_TEST ${NLA}-test)

# Create executable
add_executable(${NLA_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${NLA_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${NLA_TEST}
    COMMAND ${NLA_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(NLA_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${NLA_TEST_PROCS})
		set(NLA_PARALLEL_TEST ${NLA_TEST}-${PROCS}-procs)
		add_test(NAME ${NLA_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${NLA_TEST}${EXE_EXT}>)
	endforeach()
endif()

