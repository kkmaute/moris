# Tools Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Debug.cpp
    cl_Geometry.cpp
    cl_Interpolation.cpp)

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM}
    ${ALG}
    ${TOL})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${TOL}-lib
    ${MORIS_LDLIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS} )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${TOL}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${TOL}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${TOL}-test
    COMMAND ${TOL}-test${EXE_EXT} )
