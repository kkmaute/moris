# HMR Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_HMR_Background_Mesh.cpp
    cl_HMR_Background_Mesh_Private.cpp
    cl_HMR_BSpline_Mesh.cpp
    cl_HMR_Lagrange_Mesh.cpp
    cl_HMR_T_Matrix_Private.cpp
    cl_HMR_L2.cpp )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${COM}
    ${HMR}
    ${MTK}
    ${ALG} #MTK
    ${TOL} #MTK
    ${XTK} #MTK
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    ${COM}-lib
    ${HMR}-lib
    ${MDL}-lib
    ${MSI}-lib
    ${INT}-lib
    ${DLA}-lib
    ${MTK}-lib
    ${STK}-lib
    ${TOL}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${HMR}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${HMR}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${HMR}-test
    COMMAND ${HMR}-test${EXE_EXT} )

