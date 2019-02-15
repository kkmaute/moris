# Integration and Interpolation Tests -------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_FEM_Interpolation_Function_Lagrange_Quad4.cpp
    cl_FEM_Interpolation_Function_Lagrange_Quad8.cpp
    cl_FEM_Interpolation_Function_Lagrange_Quad9.cpp
    cl_FEM_Interpolation_Function_Lagrange_Quad16.cpp
    cl_FEM_Interpolation_Function_Lagrange_Hex8.cpp
    cl_FEM_Interpolation_Function_Lagrange_Hex20.cpp
    cl_FEM_Interpolation_Function_Lagrange_Hex27.cpp
    cl_FEM_Interpolation_Function_Lagrange_Hex64.cpp
    cl_FEM_Geometry_Interpolator_Test.cpp
    cl_FEM_Field_Interpolator_Test.cpp
    cl_FEM_Field_Integrator_Test.cpp
    cl_FEM_Field_IWG_Helmholtz.cpp
    cl_FEM_Field_IWG_Hamilton_Jacobi.cpp
    cl_FEM_Field_IWG_Olsson_CLS.cpp
    cl_FEM_Field_IWG_Test.cpp
    )

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${FEM}/${INT}
    ${MTK}
    ${ALG}
    ${FEM}/${MSI})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MSI}-lib
    ${INT}-lib
    ${DLA}-lib
    ${MTK}-lib
    ${TOL}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(INT_TEST ${INT}-test)

# Create executable
add_executable(${INT_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${INT_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${INT_TEST}
    COMMAND ${INT_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(INT_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${INT_TEST_PROCS})
		set(INT_PARALLEL_TEST ${INT_TEST}-${PROCS}-procs)
		add_test(NAME ${INT_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${INT_TEST}${EXE_EXT}>)
	endforeach()
endif()
