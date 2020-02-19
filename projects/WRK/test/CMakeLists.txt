# Parameters Tests ----------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES 
    test_main.cpp
    UT_WRK_Dummy.cpp )

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${WRK}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(WRK_TEST ${WRK}-test)

# Create executable
add_executable(${WRK_TEST} ${TEST_SOURCES})
target_link_libraries(${WRK_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${WRK_TEST} PROPERTIES OUTPUT_NAME ${WRK_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(WRK_TEST_PROCS ${MORIS_TEST_PROCS})
    foreach(PROCS ${WRK_TEST_PROCS})
        set(WRK_PARALLEL_TEST ${WRK_TEST}-${PROCS}-procs)
        add_test(NAME ${WRK_PARALLEL_TEST} 
            COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} $<TARGET_FILE:${WRK_TEST}>)
    endforeach()
endif()

