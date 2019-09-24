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
    cl_FEM_Geometry_Interpolator_3rd_Order_Derivs_Test.cpp
    cl_FEM_Field_Interpolator_Test.cpp
    cl_FEM_Integrator_Test.cpp
    cl_FEM_Element_L2.cpp
    cl_FEM_Element_Diffusion_Test.cpp
    cl_FEM_Side_Geometry_Interpolation_Test.cpp
    cl_FEM_Time_Side_Geometry_Interpolation_Test.cpp
    cl_FEM_Integration_Mesh_Test.cpp
    cl_FEM_Integration_Mesh_Test2.cpp
    cl_FEM_Double_Sided_Sideset_Test.cpp
    cl_FEM_Element_HJ_Test.cpp
    UT_FEM_Set.cpp
    UT_FEM_Property.cpp
    UT_FEM_Property_User_Defined_Info.cpp
    UT_FEM_IWG.cpp
    UT_FEM_IWG_User_Defined_Info.cpp
    UT_FEM_IWG_Spatial_Diffusion.cpp
    UT_FEM_CM.cpp
    UT_FEM_IWG_Spatial_Diffusion_Interface.cpp
    )

# List test dependencies
set(TEST_DEPENDENCIES
	test-libs
    ${INT}-lib
    ${ALG}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(INT_TEST ${INT}-test)

# Create executable
add_executable(${INT_TEST} ${TEST_SOURCES})
target_link_libraries(${INT_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${INT_TEST} PROPERTIES OUTPUT_NAME ${INT_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(INT_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${INT_TEST_PROCS})
		set(INT_PARALLEL_TEST ${INT_TEST}-${PROCS}-procs)
		add_test(NAME ${INT_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} $<TARGET_FILE:${INT_TEST}>)
	endforeach()
endif()
