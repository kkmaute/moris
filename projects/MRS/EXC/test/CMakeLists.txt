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
    ${LNA}
    ${IOS}
    ${COM} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../../${TEST_INCLUDE}/src)
endforeach()

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(../../../snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${EXC}-test ${TEST_SOURCES})
target_link_libraries(${EXC}-test ${TEST_DEPENDENCIES})
