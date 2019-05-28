# IOS Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Logger.cpp
    fn_to_stdio.cpp )

# List test includes
set(TEST_INCLUDES
    ${IOS}
    ${COM})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes
include_directories(${MORIS_DIR}/snippets/ios)

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib)

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(IOS_TEST ${IOS}-test)

# Create executable
add_executable(${IOS_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${IOS_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${IOS_TEST}
    COMMAND ${IOS_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(IOS_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${IOS_TEST_PROCS})
		set(IOS_PARALLEL_TEST ${IOS_TEST}-${PROCS}-procs)
		add_test(NAME ${IOS_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${IOS_TEST}${EXE_EXT}>)
	endforeach()
endif()
