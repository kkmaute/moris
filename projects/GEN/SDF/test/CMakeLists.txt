# SDF Tests ---------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_SDF_Triangle.cpp 
    cl_SDF_Core.cpp)

# List test includes
set(TEST_INCLUDES
    #${GEN}/${SDF}
    #${MTK}
    #${LINALG}
    #${COM}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List snippet includes
#set(SNIPPET_INCLUDES
#    )

#foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
#    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
#endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
    ${SDF}-lib
    ${MTK}-lib
    #${COM}-lib
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_TRILINOS_LIBS}
    #${MORIS_ARPACK_LIBS}
    #${MORIS_BOOST_LIBS}
    ${MORIS_BASE_LIBS}
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(SDF_TEST ${SDF}-test)

# Create executable
add_executable(${SDF_TEST} ${TEST_SOURCES})
target_link_libraries(${SDF_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${SDF_TEST} PROPERTIES OUTPUT_NAME ${SDF_TEST}${EXE_EXT})

# Mark as test
add_test(NAME ${SDF_TEST}
    COMMAND $<TARGET_FILE:${SDF_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(SDF_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${SDF_TEST_PROCS})
		set(SDF_PARALLEL_TEST ${SDF_TEST}-${PROCS}-procs)
		add_test(NAME ${SDF_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${SDF_TEST}>)
	endforeach()
endif()

