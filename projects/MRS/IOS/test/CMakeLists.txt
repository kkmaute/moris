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

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${IOS}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${IOS}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${IOS}-test
    COMMAND ${IOS}-test${EXE_EXT} )
