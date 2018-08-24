# Distributed Linear Algebra Tests ----------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Dist_Vector_Test.cpp #> VectorType error
    cl_Linear_Solver_Test.cpp
    cl_Map_Epetra.cpp
    cl_Sparse_Matrix_Test.cpp
    map.mtx
    matrix.mtx
    vector.mtx )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${COM}
    ${DLA}
    ${ALG} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${DLA}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${DLA}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${DLA}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${DLA}-test
    COMMAND ${DLA}-test${EXE_EXT} )

