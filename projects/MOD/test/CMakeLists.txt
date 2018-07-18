# Model Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Bspline.cpp
    EgnManager.cpp )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM}
    ${MOD}
    ${ALG})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MOD}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${MOD}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MOD}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MOD}-test
    COMMAND ${MOD}-test${EXE_EXT} )

