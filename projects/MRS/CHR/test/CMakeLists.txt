# Chronos Tests -----------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Stopwatch.cpp
    cl_Waitbar.cpp )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../../${TEST_INCLUDE}/src)
endforeach()

# List additonal includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MORIS_LDLIBS})

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${CHR}-test ${TEST_SOURCES})
target_link_libraries(${CHR}-test ${TEST_DEPENDENCIES})
