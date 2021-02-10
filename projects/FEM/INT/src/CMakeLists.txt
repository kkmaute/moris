# Integration and Interpolation Source ------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${INT}_VERSION ${MORIS_VERSION})

# List header dependencies                                                                                                                                                                                                                                                     
set(INT_INCLUDES
	"Advection"
	"Compressible_Fluid"
	"Compressible_Fluid_Unified"
	"Diffusion"
	"Elements"
    "Hamilton_Jacobi"
    "Helmholtz"
    "Incompressible_Fluid"
    "Integration_Coeffs"
    "Interpolation_Functions"
    "IQIs"
    "Linear_Elasticity"
    "Turbulence")

set(INT_INTERNAL_INCLUDES)

# Include depedency directories
foreach(INT_INCLUDES ${INT_INCLUDES})
    list(APPEND INT_INTERNAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/${INT_INCLUDES}")
	include_directories(${MORIS_PACKAGE_DIR}/FEM/INT/src/${INT_INCLUDES})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# INT Library

# List header files
set(HEADERS
    cl_FEM_Element_Factory.hpp
    cl_FEM_ElementProxy.hpp
    cl_FEM_NodeProxy.hpp
    cl_FEM_Cluster.hpp
    cl_FEM_Element.hpp
    cl_FEM_Set.hpp
    Elements/cl_FEM_Element_Bulk.hpp
    Elements/cl_FEM_Element_Sideset.hpp
    Elements/cl_FEM_Element_Double_Sideset.hpp
    Elements/cl_FEM_Element_Time_Sideset.hpp
    Elements/cl_FEM_Element_Time_Boundary.hpp
    cl_FEM_IWG_Factory.hpp
    cl_FEM_IWG.hpp
    cl_FEM_IWG_L2.hpp
    
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.hpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Bulk.hpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Dirichlet.hpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.hpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.hpp    
    
    Advection/cl_FEM_IWG_Advection_Bulk.hpp

    Diffusion/cl_FEM_IWG_Diffusion_Bulk.hpp
    Diffusion/cl_FEM_IWG_Diffusion_Dirichlet_Nitsche.hpp
    Diffusion/cl_FEM_IWG_Diffusion_Neumann.hpp
    Diffusion/cl_FEM_IWG_Diffusion_Robin.hpp
    Diffusion/cl_FEM_IWG_Diffusion_Radiation.hpp
    Diffusion/cl_FEM_IWG_Diffusion_Interface.hpp
    Diffusion/cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.hpp
    
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Velocity_Bulk.hpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.hpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.hpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche.hpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche.hpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Pressure_Neumann.hpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Velocity_Interface.hpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Pressure_Interface.hpp
    
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Density_Bulk.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Velocity_Bulk.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Temperature_Bulk.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Advective_Momentum_Flux_Boundary.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Advective_Energy_Flux_Boundary.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Mass_Flux_Neumann.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Traction_Neumann.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Heat_Flux_Neumann.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.hpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Temperature_Dirichlet_Nitsche.hpp
    
    Turbulence/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp
    Turbulence/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.hpp
    Turbulence/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.hpp

    #Helmholtz/cl_FEM_IWG_Helmholtz_Bulk.hpp
    #Helmholtz/cl_FEM_IWG_Helmholtz_Bulk2.hpp
    #Helmholtz/cl_FEM_IWG_Helmholtz_Interface.hpp

    #Hamilton_Jacobi/cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp
    #Hamilton_Jacobi/cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp
    Hamilton_Jacobi/cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp

    cl_FEM_IWG_FS_Struc_Interface.hpp
    cl_FEM_IWG_Time_Continuity_Dof.hpp
    cl_FEM_IWG_Ghost_Normal_Field.hpp
    #cl_FEM_IWG_LSNormal_Bulk.hpp
    #cl_FEM_IWG_Olsson_CLS_Bulk.hpp
    #cl_FEM_IWG_Olsson_CLS_Interface.hpp
    cl_FEM_Enums.hpp
    cl_FEM_Geometry_Interpolator.hpp
    
    cl_FEM_Integration_Coeffs.hpp    
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Bar_1.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Bar_2.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Bar_3.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Bar_4.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Bar_5.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Bar_6.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Base.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Hex_2x2x2.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Hex_3x3x3.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Hex_4x4x4.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Hex_5x5x5.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Quad_2x2.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Quad_3x3.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Quad_4x4.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Quad_5x5.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_1.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_3.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_4.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_6.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_7.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_12.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_13.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_16.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_19.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tri_25.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tet_1.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tet_4.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tet_5.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tet_11.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tet_15.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tet_20.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tet_35.hpp
    Integration_Coeffs/cl_FEM_Integration_Coeffs_Tet_56.hpp

    cl_FEM_Integration_Rule.hpp
    cl_FEM_Integrator.hpp
    cl_FEM_Integrator_Test_Polynomial.hpp
    cl_FEM_Interpolation_Function_Base.hpp
    cl_FEM_Interpolation_Function_Factory.hpp
    cl_FEM_Interpolation_Function.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Bar1.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Bar2.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Bar3.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Bar4.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Hex27.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Hex64.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Hex8.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Quad16.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Quad4.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Quad8.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Quad9.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Tri3.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Tri6.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Tri10.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Tet4.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Tet10.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Lagrange_Tet20.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Constant_Bar2.hpp
    Interpolation_Functions/cl_FEM_Interpolation_Function_Constant_Point.hpp
    cl_FEM_Interpolation_Rule.hpp
    cl_FEM_Node_Base.hpp
    cl_FEM_Node.hpp
    cl_FEM_Field_Interpolator.hpp
    cl_FEM_Field_Interpolator_Manager.hpp
    cl_FEM_Field.hpp        
    fn_FEM_Rotation_Matrix.hpp
    cl_FEM_Property.hpp
    
    cl_FEM_Material_Model.hpp
    cl_FEM_MM_Perfect_Gas.hpp
    #cl_FEM_MM_Van_Der_Waals_Fluid.hpp
    cl_FEM_MM_Factory.hpp

    cl_FEM_Constitutive_Model.hpp
    cl_FEM_CM_Factory.hpp
    Diffusion/cl_FEM_CM_Diffusion_Linear_Isotropic.hpp
    Diffusion/cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.hpp
    Linear_Elasticity/cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp
    Linear_Elasticity/cl_FEM_CM_Struc_Linear_Isotropic.hpp
    Incompressible_Fluid/cl_FEM_CM_Fluid_Incompressible.hpp
    Turbulence/cl_FEM_CM_Fluid_Turbulence.hpp
    Compressible_Fluid/cl_FEM_CM_Fluid_Compressible_Ideal.hpp
    Compressible_Fluid_Unified/cl_FEM_CM_Compressible_Newtonian_Fluid.hpp
    Compressible_Fluid/cl_FEM_CM_Fluid_Compressible_Van_der_Waals.hpp
    
    cl_FEM_IQI.hpp
    cl_FEM_IQI_Factory.hpp
    IQIs/cl_FEM_IQI_Volume.hpp
    IQIs/cl_FEM_IQI_Strain_Energy.hpp
    IQIs/cl_FEM_IQI_Stress.hpp
    IQIs/cl_FEM_IQI_Dof.hpp
    IQIs/cl_FEM_IQI_Max_Dof.hpp
    IQIs/cl_FEM_IQI_Max_Stress.hpp
    IQIs/cl_FEM_IQI_Property.hpp
    IQIs/cl_FEM_IQI_L2_Error_Analytic.hpp
    IQIs/cl_FEM_IQI_H1_Error_Analytic.hpp
    IQIs/cl_FEM_IQI_H1_Semi_Error.hpp
    IQIs/cl_FEM_IQI_J_Integral.hpp
    IQIs/cl_FEM_IQI_Volume_Fraction.hpp
    IQIs/cl_FEM_IQI_2D_Drag_Lift_Coefficient.hpp
    IQIs/cl_FEM_IQI_Latent_Heat_Absorption.hpp
    IQIs/cl_FEM_IQI_Turbulent_Kinematic_Viscosity.hpp
    IQIs/cl_FEM_IQI_Total_Pressure.hpp
    IQIs/cl_FEM_IQI_Mass_Flow.hpp
    IQIs/cl_FEM_IQI_Thermal_Energy.hpp

    cl_FEM_Stabilization_Parameter.hpp
    cl_FEM_SP_Factory.hpp
    cl_FEM_SP_Dirichlet_Nitsche.hpp
    cl_FEM_SP_Nitsche_Interface.hpp
    cl_FEM_SP_Reciprocal_Total_Volume.hpp
    cl_FEM_SP_Convective_Ghost.hpp
    cl_FEM_SP_Time_Velocity_Ghost.hpp   
    Advection/cl_FEM_SP_SUPG_Advection.hpp
    Diffusion/cl_FEM_SP_GGLS_Diffusion.hpp
    Incompressible_Fluid/cl_FEM_SP_Velocity_Dirichlet_Nitsche.hpp
    Incompressible_Fluid/cl_FEM_SP_Pressure_Ghost.hpp
    Incompressible_Fluid/cl_FEM_SP_Incompressible_Flow.hpp
    Incompressible_Fluid/cl_FEM_SP_Viscous_Ghost.hpp
    Compressible_Fluid/cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.hpp
    Turbulence/cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.hpp   
    Turbulence/cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp  
    Turbulence/cl_FEM_SP_Turbulence_Dirichlet_Nitsche.hpp   
    Linear_Elasticity/cl_FEM_SP_Ghost_Displacement.hpp
    Linear_Elasticity/cl_FEM_SP_Ghost_Virtual_Work.hpp 
    Linear_Elasticity/cl_FEM_SP_Stab_Penalty_Contact.hpp
    Linear_Elasticity/cl_FEM_SP_Penalty_Contact.hpp 
    
    cl_FEM_Set_User_Info.hpp
    cl_FEM_Phase_User_Info.hpp
    cl_FEM_Interpolation_Element.hpp
    cl_FEM_Model.hpp
    fn_FEM_CM_Phase_State_Functions.hpp
    fn_FEM_Check.hpp
    fn_FEM_FD_Scheme.hpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Library

# List library source files
set(LIB_SOURCES
    cl_FEM_Set.cpp
    cl_FEM_Cluster.cpp
    cl_FEM_Element_Factory.cpp
    Elements/cl_FEM_Element_Bulk.cpp
    Elements/cl_FEM_Element_Sideset.cpp
    Elements/cl_FEM_Element_Double_Sideset.cpp
    Elements/cl_FEM_Element_Time_Sideset.cpp
    Elements/cl_FEM_Element_Time_Boundary.cpp
    
    cl_FEM_IWG_Factory.cpp
    cl_FEM_IWG_L2.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Interface.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Bulk.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Dirichlet.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.cpp
    Linear_Elasticity/cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.cpp
    
    Advection/cl_FEM_IWG_Advection_Bulk.cpp

    Diffusion/cl_FEM_IWG_Diffusion_Bulk.cpp
    Diffusion/cl_FEM_IWG_Diffusion_Dirichlet_Nitsche.cpp
    Diffusion/cl_FEM_IWG_Diffusion_Neumann.cpp
    Diffusion/cl_FEM_IWG_Diffusion_Robin.cpp
    Diffusion/cl_FEM_IWG_Diffusion_Radiation.cpp
    Diffusion/cl_FEM_IWG_Diffusion_Interface.cpp
    Diffusion/cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.cpp
    
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Velocity_Bulk.cpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.cpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.cpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche.cpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche.cpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Pressure_Neumann.cpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Velocity_Interface.cpp
    Incompressible_Fluid/cl_FEM_IWG_Incompressible_NS_Pressure_Interface.cpp
    
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Density_Bulk.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Velocity_Bulk.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Temperature_Bulk.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Advective_Momentum_Flux_Boundary.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Advective_Energy_Flux_Boundary.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Mass_Flux_Neumann.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Traction_Neumann.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Heat_Flux_Neumann.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.cpp
    Compressible_Fluid/cl_FEM_IWG_Compressible_NS_Temperature_Dirichlet_Nitsche.cpp
    
    Turbulence/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.cpp
    Turbulence/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.cpp
    Turbulence/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.cpp

    #Helmholtz/cl_FEM_IWG_Helmholtz_Bulk.cpp
    #Helmholtz/cl_FEM_IWG_Helmholtz_Bulk2.cpp
    #Helmholtz/cl_FEM_IWG_Helmholtz_Interface.cpp

    #Hamilton_Jacobi/cl_FEM_IWG_Hamilton_Jacobi_Bulk.cpp
    #Hamilton_Jacobi/cl_FEM_IWG_Hamilton_Jacobi_Bulk2.cpp
    Hamilton_Jacobi/cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.cpp

    cl_FEM_IWG_FS_Struc_Interface.cpp
    cl_FEM_IWG_Time_Continuity_Dof.cpp
    cl_FEM_IWG_Ghost_Normal_Field.cpp
    #cl_FEM_IWG_LSNormal_Bulk.cpp
    #cl_FEM_IWG_Olsson_CLS_Bulk.cpp
    #cl_FEM_IWG_Olsson_CLS_Interface.cpp
    cl_FEM_IWG.cpp
    
    cl_FEM_Interpolation_Function_Factory.cpp
    cl_FEM_Interpolation_Rule.cpp
    cl_FEM_Geometry_Interpolator.cpp
    cl_FEM_Field_Interpolator.cpp
    cl_FEM_Field_Interpolator_Manager.cpp
    cl_FEM_Field.cpp  
    cl_FEM_Integration_Rule.cpp
    cl_FEM_Integrator.cpp
    cl_FEM_Integrator_Test_Polynomial.cpp
    cl_FEM_Property.cpp
    
    cl_FEM_Material_Model.cpp
    cl_FEM_MM_Factory.cpp
    cl_FEM_Material_Model_Density_Functions.cpp
    cl_FEM_Material_Model_Pressure_Functions.cpp
    cl_FEM_Material_Model_Temperature_Functions.cpp
    Compressible_Fluid_Unified/cl_FEM_MM_Perfect_Gas.cpp
    #Compressible_Fluid_Unified/cl_FEM_MM_Van_Der_Waals_Fluid.cpp

    cl_FEM_Constitutive_Model.cpp
    cl_FEM_CM_Factory.cpp
    Diffusion/cl_FEM_CM_Diffusion_Linear_Isotropic.cpp
    Diffusion/cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.cpp
    Linear_Elasticity/cl_FEM_CM_Struc_Nonlinear_Isotropic.cpp
    Linear_Elasticity/cl_FEM_CM_Struc_Linear_Isotropic.cpp
    Incompressible_Fluid/cl_FEM_CM_Fluid_Incompressible.cpp
    Turbulence/cl_FEM_CM_Fluid_Turbulence.cpp
    Compressible_Fluid/cl_FEM_CM_Fluid_Compressible_Ideal.cpp
    Compressible_Fluid/cl_FEM_CM_Fluid_Compressible_Ideal_Test_Tractions.cpp
    Compressible_Fluid/cl_FEM_CM_Fluid_Compressible_Van_der_Waals.cpp
    Compressible_Fluid/cl_FEM_CM_Fluid_Compressible_Van_der_Waals_Test_Tractions.cpp
    Compressible_Fluid_Unified/cl_FEM_CM_Compressible_Newtonian_Fluid.cpp
    Compressible_Fluid_Unified/cl_FEM_CM_Compressible_Newtonian_Fluid_Test_Tractions.cpp    
    
    cl_FEM_Stabilization_Parameter.cpp
    cl_FEM_SP_Factory.cpp
    cl_FEM_SP_Dirichlet_Nitsche.cpp
    cl_FEM_SP_Nitsche_Interface.cpp
    cl_FEM_SP_Reciprocal_Total_Volume.cpp
    cl_FEM_SP_Convective_Ghost.cpp
    cl_FEM_SP_Time_Velocity_Ghost.cpp   
    Advection/cl_FEM_SP_SUPG_Advection.cpp
    Diffusion/cl_FEM_SP_GGLS_Diffusion.cpp
    Incompressible_Fluid/cl_FEM_SP_Velocity_Dirichlet_Nitsche.cpp
    Incompressible_Fluid/cl_FEM_SP_Pressure_Ghost.cpp
    Incompressible_Fluid/cl_FEM_SP_Incompressible_Flow.cpp
    Incompressible_Fluid/cl_FEM_SP_Viscous_Ghost.cpp
    Compressible_Fluid/cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.cpp
    Turbulence/cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.cpp   
    Turbulence/cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.cpp  
    Turbulence/cl_FEM_SP_Turbulence_Dirichlet_Nitsche.cpp   
    Linear_Elasticity/cl_FEM_SP_Ghost_Displacement.cpp
    Linear_Elasticity/cl_FEM_SP_Ghost_Virtual_Work.cpp 
    Linear_Elasticity/cl_FEM_SP_Stab_Penalty_Contact.cpp
    Linear_Elasticity/cl_FEM_SP_Penalty_Contact.cpp       
    
    cl_FEM_IQI.cpp
    cl_FEM_IQI_Factory.cpp
    IQIs/cl_FEM_IQI_Volume.cpp
    IQIs/cl_FEM_IQI_Strain_Energy.cpp
    IQIs/cl_FEM_IQI_Stress.cpp
    IQIs/cl_FEM_IQI_Dof.cpp
    IQIs/cl_FEM_IQI_Max_Dof.cpp
    IQIs/cl_FEM_IQI_Max_Stress.cpp
    IQIs/cl_FEM_IQI_Property.cpp
    IQIs/cl_FEM_IQI_L2_Error_Analytic.cpp
    IQIs/cl_FEM_IQI_H1_Error_Analytic.cpp
    IQIs/cl_FEM_IQI_H1_Semi_Error.cpp
    IQIs/cl_FEM_IQI_J_Integral.cpp
    IQIs/cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp
    IQIs/cl_FEM_IQI_Volume_Fraction.cpp
    IQIs/cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp
    IQIs/cl_FEM_IQI_Latent_Heat_Absorption.cpp
    IQIs/cl_FEM_IQI_Turbulent_Kinematic_Viscosity.cpp
    IQIs/cl_FEM_IQI_Total_Pressure.cpp
    IQIs/cl_FEM_IQI_Mass_Flow.cpp
    IQIs/cl_FEM_IQI_Thermal_Energy.cpp
    
    cl_FEM_Interpolation_Element.cpp
    cl_FEM_Model.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    ${GEN_CORE}-lib
    ${LINALG}-lib
    ${MTK}-lib
    ${MSI}-lib
    ${VIS}-lib
    ${PRM}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${INT}-lib STATIC ${LIB_SOURCES})
target_include_directories(${INT}-lib PUBLIC
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${INT_INTERNAL_INCLUDES};>"
	$<INSTALL_INTERFACE:${${INT}_HEADER_INSTALL_DIR}> )
target_link_libraries(${INT}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${INT}-lib PROPERTIES OUTPUT_NAME ${INT})

foreach(TPL ${INT_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${INT}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${INT}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${INT}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${INT_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND INT_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND INT_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(INT_COMPILE_DEFINITIONS ${INT}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${INT}Config.cmake.in
    ${${INT}_CONFIG_BUILD_DIR}/${INT}Config.cmake
    INSTALL_DESTINATION 
    ${${INT}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${INT}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${INT}_CONFIG_BUILD_DIR}/${INT}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Installation

# Install *Config and *ConfigVersion files
install(FILES ${${INT}_CONFIG_BUILD_DIR}/${INT}Config.cmake 
    ${${INT}_CONFIG_BUILD_DIR}/${INT}ConfigVersion.cmake
    DESTINATION ${${INT}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${INT}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${INT}-lib
    DESTINATION ${${INT}_LIB_INSTALL_DIR}
    EXPORT ${INT}Targets )

# Export library into *Targets file
install(EXPORT ${INT}Targets
    DESTINATION ${${INT}_CONFIG_INSTALL_DIR} )
