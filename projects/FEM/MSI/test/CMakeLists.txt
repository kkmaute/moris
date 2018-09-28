# MSI ----------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_MSI_Dof_Manager_Test.cpp
    cl_MSI_Equation_Object_Test.cpp
    cl_MSI_Pdof_Host_Test.cpp
    cl_MSI_Model_Solver_Interface_Test.cpp)

# List test includes
set(TEST_INCLUDES
    ${DLA}
    ${NLA}
    ${COM}
    ${LINALG}
    ${ALG}
    ${MTK}
    ${INT}
)

include_directories(${MORIS_INCDIRS})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

    include_directories(../../INT/src)

# List additional includes
include_directories(../src)
# include_directories(${MORIS_PACKAGE_DIR}/IOS/src)
# include_directories(${MORIS_PACKAGE_DIR}/DLA/src)
# include_directories(${MORIS_PACKAGE_DIR}/FEM/INT/src)

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${NLA}-lib
    ${MSI}-lib
    ${INT}-lib
    ${MORIS_BASE_LIBS}
#    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${MSI}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MSI}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MSI}-test
    COMMAND ${MSI}-test${EXE_EXT} )




