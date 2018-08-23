# Linear Algebra Tests --------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES )

# List test includes
set(TEST_INCLUDES
    ${LINALG} )
    
    
foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
#add_executable(${LINALG}-test${EXE_EXT} ${TEST_SOURCES})
#target_link_libraries(${LINALG}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${LINALG}-test
    COMMAND ${LINALG}-test${EXE_EXT} )
