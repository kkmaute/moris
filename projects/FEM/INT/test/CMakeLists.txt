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
    cl_FEM_Interpolation_Integration_Lagrange_Tri_Tet.cpp
    cl_FEM_Geometry_Interpolator_Test.cpp
    cl_FEM_Field_Interpolator_Test.cpp
    cl_FEM_Integrator_Test.cpp
    cl_FEM_Element_L2.cpp
    cl_FEM_Element_Diffusion_Test.cpp
    cl_FEM_Side_Geometry_Interpolation_Test.cpp
    cl_FEM_Integration_Mesh_Test.cpp
    cl_FEM_Integration_Mesh_Test2.cpp
    cl_FEM_Double_Sided_Sideset_Test.cpp
    cl_FEM_Element_HJ_Test.cpp
    cl_FEM_Element_Block_Test.cpp
    cl_FEM_Model_Param_Interpolator_Test.cpp
    )

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${MTK}
    ${ALG}
	${HMR}
    ${FEM}/${INT}
    ${FEM}/${MSI}
    ${FEM}/${MDL}
    ${SOL}/${DLA}
    ${SOL}/${NLA}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MSI}-lib
    ${MDL}-lib
    ${INT}-lib
    ${DLA}-lib
    ${NLA}-lib
    ${MTK}-lib
    ${TOL}-lib
	${HMR}-lib
	${MAP}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

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
