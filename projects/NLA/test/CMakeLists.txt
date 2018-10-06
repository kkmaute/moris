# Non Linear Algebra Tests ----------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_NLA_Newton_Solver_Test.cpp)

# List test includes
set(TEST_INCLUDES
    ${DLA}
    ${LINALG}
    ${COM}
    ${NLA}
    ${ALG} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${DLA}-lib
    ${COM}-lib
    ${NLA}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${NLA}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${NLA}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${NLA}-test
    COMMAND ${NLA}-test${EXE_EXT} )

