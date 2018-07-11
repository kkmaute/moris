# Algorithms Tests --------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    faux_main.cpp
    fn_equal_to.cpp
    fn_num_digits.cpp)

# List test includes
set(TEST_INCLUDES
    ${ALG}
    ${LNA}
    ${IOS}
    ${COM} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../${TEST_INCLUDE}/src)
endforeach()

# List snippet includes
set(SNIPPET_INCLUDES
    algorithms )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(../../../snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${ALG}-test ${TEST_SOURCES})
target_link_libraries(${ALG}-test ${TEST_DEPENDENCIES})
