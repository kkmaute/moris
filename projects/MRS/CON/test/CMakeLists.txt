# Containers Tests --------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Array.cpp
    cl_Bi_Map.cpp
    cl_Cell.cpp
    cl_Dist_Map.cpp
    cl_Map.cpp
    cl_Param_List.cpp
    cl_Tuple.cpp )

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${ALG})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_BOOST_LIBS})

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(CON_TEST ${CON}-test)

# Create executable
add_executable(${CON_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${CON_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${CON_TEST}
    COMMAND ${CON_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(CON_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${CON_TEST_PROCS})
		set(CON_PARALLEL_TEST ${CON_TEST}-${PROCS}-procs)
		add_test(NAME ${CON_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${CON_TEST}${EXE_EXT}>)
	endforeach()
endif()
