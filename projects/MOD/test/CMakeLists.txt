# Model Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Bspline.cpp
    EgnManager.cpp )

# List test includes
set(TEST_INCLUDES
    ${COM}
    ${MOD}
    ${ALG})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MOD}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MOD_TEST ${MOD}-test)

# Create executable
add_executable(${MOD_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MOD_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MOD_TEST}
    COMMAND ${MOD_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(MOD_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${MOD_TEST_PROCS})
		set(MOD_PARALLEL_TEST ${MOD_TEST}-${PROCS}-procs)
		add_test(NAME ${MOD_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${MOD_TEST}${EXE_EXT}>)
	endforeach()
endif()

