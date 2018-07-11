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
    ${DLA}
    ${LNA}
    ${IOS}
    ${COM}
    ${ALG} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../${TEST_INCLUDE}/src)
endforeach()

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(../../../snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${DLA}-lib
    ${MORIS_LDLIBS})

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${DLA}-test ${TEST_SOURCES})
target_link_libraries(${DLA}-test ${TEST_DEPENDENCIES})
