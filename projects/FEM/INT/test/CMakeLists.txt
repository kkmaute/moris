#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

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
    #UT_FEM_Element.cpp
    UT_FEM_Property.cpp
    UT_FEM_IWG.cpp
    
    UT_FEM_CM.cpp
    UT_FEM_CM_Diff_Lin_Iso.cpp
    UT_FEM_CM_Diff_PC.cpp
    UT_FEM_CM_Diff_Lin_Iso_Turbulence.cpp
    UT_FEM_CM_Fluid.cpp
    UT_FEM_CM_Compressible_Fluid_Ideal.cpp
    UT_FEM_CM_Compressible_Newtonian_Fluid.cpp
    UT_FEM_CM_Compressible_Fluid_VdW.cpp
    UT_FEM_CM_Struc_Linear_Isotropic_Damage.cpp
    UT_FEM_CM_Struc_Linear_Isotropic.cpp
    UT_FEM_CM_Struc_Linear_Isotropic_Axisymmetric.cpp
    UT_FEM_CM_Spalart_Allmaras_Turbulence.cpp
    UT_FEM_CM_Struct_Linear_MT.cpp
    UT_FEM_CM_Struc_NonLinear_Saint_Venant_Kirchhoff.cpp
    UT_FEM_CM_Struc_NonLinear_Neo_Hookean.cpp
    
    UT_FEM_SP_Fluid.cpp
    UT_FEM_SP_YZBeta_Advection.cpp
    UT_FEM_SP_Crosswind.cpp
    
    UT_FEM_IWG_Advection_Bulk.cpp
    UT_FEM_IWG_Diffusion_Bulk.cpp
    UT_FEM_IWG_Diffusion_Phase_Change_Bulk.cpp
    UT_FEM_IWG_GGLS_Diffusion_Phase_Change.cpp
    UT_FEM_IWG_Diffusion_Dirichlet.cpp
    UT_FEM_IWG_Diffusion_Robin_Nitsche.cpp
    UT_FEM_IWG_Diffusion_Interface.cpp
    UT_FEM_IWG_Diffusion_Ghost.cpp
    UT_FEM_IWG_Diffusion_Virtual_Work_Ghost.cpp
    UT_FEM_IWG_Ghost_Analytic.cpp
    UT_FEM_IWG_Struc_Linear_Bulk.cpp
    UT_FEM_IWG_Struc_Linear_Pressure_Bulk.cpp
    UT_FEM_IWG_Struc_Linear_Dirichlet.cpp
    UT_FEM_IWG_Struc_Linear_Pressure_Dirichlet.cpp
    UT_FEM_IWG_Struc_Linear_Interface.cpp
    UT_FEM_IWG_Struc_Linear_Contact_Nitsche.cpp
    
    UT_FEM_IWG_Struc_Linear_Interface_SLM_Constraint_and_Mixed.cpp
    UT_FEM_IWG_Struc_Linear_Interface_SLM_L2.cpp
    
    UT_FEM_IWG_Nonlocal_Bulk.cpp
    UT_FEM_IWG_Nonlocal_Interface.cpp
    UT_FEM_IWG_History_Bulk.cpp
    
    #UT_FEM_IWG_Struc_Linear_Ghost.cpp
    UT_FEM_IWG_Struc_Linear_Virtual_Work_Ghost.cpp

    UT_FEM_IWG_Struc_Linear_Axisymmetric_Bulk.cpp
    UT_FEM_IWG_Struc_Linear_Axisymmetric_Pressure_Bulk.cpp
    UT_FEM_IWG_Struc_Linear_Axisymmetric_Dirichlet.cpp
    UT_FEM_IWG_Struc_Linear_Axisymmetric_Pressure_Dirichlet.cpp
    UT_FEM_IWG_Struc_Linear_Axisymmetric_Interface.cpp
    UT_FEM_IWG_Struc_Linear_Axisymmetric_Virtual_Work_Ghost.cpp
    
    UT_FEM_IWG_Struc_Nonlinear_Bulk.cpp
    UT_FEM_IWG_Struc_Nonlinear_Dirichlet.cpp
    UT_FEM_IWG_Struc_Nonlinear_Interface.cpp
    
    UT_FEM_IWG_Incompressible_NS_Velocity_Bulk.cpp
    UT_FEM_IWG_Incompressible_NS_Pressure_Bulk.cpp
    UT_FEM_IWG_Incompressible_NS_Dirichlet_Nitsche.cpp
    UT_FEM_IWG_Incompressible_NS_SlipBoundary_Nitsche.cpp
    UT_FEM_IWG_Incompressible_NS_Interface_Nitsche.cpp
    UT_FEM_IWG_Incompressible_NS_Pressure_Neumann.cpp
    
    UT_FEM_IWG_Compressible_NS_Bulk.cpp
    UT_FEM_IWG_Compressible_NS_Bulk_Analytical.cpp
    UT_FEM_IWG_Compressible_NS_Dirichlet_Analytical.cpp
    UT_FEM_IWG_Compressible_NS_Bulk_Flux_Matrices.cpp
    UT_FEM_IWG_Compressible_NS_Boundary.cpp
    UT_FEM_IWG_Compressible_NS_Dirichlet_Nitsche.cpp

    UT_FEM_IWG_Compressible_NS_Density_Bulk.cpp
    UT_FEM_IWG_Compressible_NS_Density_Bulk_Analytical.cpp
    UT_FEM_IWG_Compressible_NS_Velocity_Bulk.cpp
    UT_FEM_IWG_Compressible_NS_Velocity_Bulk_Analytical.cpp
    UT_FEM_IWG_Compressible_NS_Temperature_Bulk.cpp
    UT_FEM_IWG_Compressible_NS_Temperature_Bulk_Analytical.cpp
    UT_FEM_IWG_Compressible_NS_Advective_Momentum_Boundary.cpp
    UT_FEM_IWG_Compressible_NS_Advective_Energy_Boundary.cpp
    UT_FEM_IWG_Compressible_NS_Neumann_Boundaries.cpp
    UT_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.cpp
    UT_FEM_IWG_Compressible_NS_Temperature_Dirichlet_Nitsche.cpp
    
    #UT_FEM_IWG_FS_Struc_Interface.cpp
    UT_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.cpp
    UT_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.cpp
    UT_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.cpp
    
    UT_FEM_IQI_Dof.cpp
    UT_FEM_IQI_Strain_Energy.cpp
    UT_FEM_IQI_Power_Dissipation.cpp
    UT_FEM_IQI_Power_Dissipation_Bulk.cpp
  
    UT_FEM_Input.cpp
    UT_FEM_Geometry_Interpolator.cpp
    UT_FEM_Integration_Rule.cpp
    
    FEM_Test_Proxy/cl_FEM_Design_Variable_Interface_Proxy.cpp
    FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp
    FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Compressible_UT.cpp
    FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Compressible_1D.cpp
    FEM_Test_Proxy/cl_FEM_Fields_for_NS_Compressible_UT.cpp
    FEM_Test_Proxy/cl_FEM_Fields_for_NS_Compressible_UT_for_BD.cpp
    FEM_Test_Proxy/cl_FEM_Flux_Matrix_Reference_Values.cpp
    FEM_Test_Proxy/cl_FEM_Inputs_for_Diffusion_UT.cpp
    FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp
    FEM_Test_Proxy/cl_FEM_Inputs_for_Ghost_UT.cpp
    FEM_Test_Proxy/fn_FEM_Convert_Dimensions.cpp
    )

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${INT}-lib
    ${ALG}-lib
    ${PRM}-lib
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
            COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${INT_TEST}>)
    endforeach()
endif()

