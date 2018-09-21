# MTK Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
	cl_Mesh.cpp
	cl_MTK_Vertex.cpp
	cl_MTK_Cell.cpp)

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${ALG}
    ${DLA}
    ${MTK}
    ${TOL} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes
include_directories(${MORIS_PACKAGE_DIR}/MTK/src/stk_impl)

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MTK}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${MTK}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MTK}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MTK}-test
    COMMAND ${MTK}-test${EXE_EXT} )