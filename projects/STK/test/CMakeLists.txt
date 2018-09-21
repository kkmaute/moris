# STK Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Mesh_From_Data.cpp
    cl_Mesh_From_String.cpp
    cl_Mesh_Input_Buckets.cpp )

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${ALG}
    ${DLA}
    ${MOD}
    ${MTK}
    ${TOL}
    ${XTK}
    ${STK} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# Additional includes
#include_directories(${MORIS_PACKAGE_DIR}/${STK}/src/Hierarchical)

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${STK}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${TOL}-lib
    ${MORIS_SUPERLU_LIBS}
    ${DLA}-lib
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${STK}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${STK}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${STK}-test
    COMMAND ${STK}-test${EXE_EXT} )
