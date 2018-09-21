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
    ${COM} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${EXC}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${EXC}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${EXC}-test
    COMMAND ${EXC}-test${EXE_EXT} )
