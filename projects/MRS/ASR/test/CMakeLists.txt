# Assert Tests ------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    assert.cpp)

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../../${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${ASR}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${ASR}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${ASR}-test
    COMMAND ${ASR}-test${EXE_EXT} )
