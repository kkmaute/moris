# Integration and Interpolation Source ------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${INT}_VERSION ${MORIS_VERSION})

# List header dependencies                                                                                                                                                                                                                                                     
set(INT_INCLUDES
	"CORE"
    "ELEM"
    "IP"
	"MM"
	"CM"
    "SP"
	"IWG"
    "IQI"
    "MEA")

set(INT_INTERNAL_INCLUDES)

# Include depedency directories
foreach(INT_INCLUDES ${INT_INCLUDES})
    list(APPEND INT_INTERNAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/${INT_INCLUDES}")
	include_directories(${MORIS_PACKAGE_DIR}/FEM/INT/src/${INT_INCLUDES})
endforeach()

# create global list for .../INT/src/ includes
set_property(GLOBAL PROPERTY INT_SRC_LIST "${INT_INTERNAL_INCLUDES}")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# INT Library

# List header files
set(HEADERS

    CORE/cl_FEM_Set.hpp
    CORE/cl_FEM_NodeProxy.hpp
    CORE/cl_FEM_Enums.hpp
    CORE/cl_FEM_Node_Base.hpp
    CORE/cl_FEM_Node.hpp
    CORE/cl_FEM_Field.hpp        
    CORE/fn_FEM_Rotation_Matrix.hpp
    CORE/cl_FEM_Property.hpp
    CORE/cl_FEM_Set_User_Info.hpp
    CORE/cl_FEM_Phase_User_Info.hpp
    CORE/cl_FEM_Model.hpp
    CORE/fn_FEM_Check.hpp
    CORE/fn_FEM_FD_Scheme.hpp

    ELEM/cl_FEM_Cluster.hpp
    ELEM/cl_FEM_Element_Factory.hpp
    ELEM/cl_FEM_ElementProxy.hpp
    ELEM/cl_FEM_Element.hpp

    ELEM/cl_FEM_Element_Bulk.hpp
    ELEM/cl_FEM_Element_Sideset.hpp
    ELEM/cl_FEM_Element_Double_Sideset.hpp
    ELEM/cl_FEM_Element_Time_Sideset.hpp
    ELEM/cl_FEM_Element_Time_Boundary.hpp


    IP/cl_FEM_Geometry_Interpolator.hpp
    IP/cl_FEM_Field_Interpolator.hpp
    IP/cl_FEM_Field_Interpolator_Manager.hpp


    MM/cl_FEM_Material_Model.hpp
    MM/cl_FEM_MM_Factory.hpp

    MM/cl_FEM_MM_Perfect_Gas.hpp
    #MM/cl_FEM_MM_Van_Der_Waals_Fluid.hpp


    CM/cl_FEM_Constitutive_Model.hpp
    CM/cl_FEM_CM_Factory.hpp

    CM/cl_FEM_CM_Diffusion_Linear_Isotropic.hpp
    CM/cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.hpp
    CM/fn_FEM_CM_Phase_State_Functions.hpp
    CM/cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp

    CM/cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp
    CM/cl_FEM_CM_Struc_Linear_Isotropic.hpp
    CM/cl_FEM_CM_Struc_Linear_MoriTanaka.hpp
    CM/cl_FEM_CM_Struc_Linear.hpp

    CM/cl_FEM_CM_Fluid_Incompressible.hpp
    CM/cl_FEM_CM_Fluid_Turbulence.hpp
    CM/cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp

    CM/cl_FEM_CM_Fluid_Compressible_Ideal.hpp
    CM/cl_FEM_CM_Fluid_Compressible_Van_der_Waals.hpp
    CM/cl_FEM_CM_Compressible_Newtonian_Fluid.hpp


    SP/cl_FEM_Stabilization_Parameter.hpp
    SP/cl_FEM_SP_Factory.hpp

    SP/cl_FEM_SP_Dirichlet_Nitsche.hpp
    SP/cl_FEM_SP_Nitsche_Interface.hpp
    SP/cl_FEM_SP_Reciprocal_Total_Volume.hpp
    SP/cl_FEM_SP_Time_Velocity_Ghost.hpp 
    SP/cl_FEM_SP_Convective_Ghost.hpp  
    SP/cl_FEM_SP_SUPG_Advection.hpp
    SP/cl_FEM_SP_YZBeta_Advection.hpp
    SP/cl_FEM_SP_GGLS_Diffusion.hpp
    SP/cl_FEM_SP_Velocity_Dirichlet_Nitsche.hpp
    SP/cl_FEM_SP_Velocity_SlipBoundary_Nitsche.hpp
    SP/cl_FEM_SP_Pressure_Ghost.hpp
    SP/cl_FEM_SP_Incompressible_Flow.hpp
    SP/cl_FEM_SP_Viscous_Ghost.hpp
    SP/cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.hpp
    SP/cl_FEM_SP_Compressible_Dirichlet_Nitsche.hpp
    SP/cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.hpp   
    SP/cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp  
    SP/cl_FEM_SP_Turbulence_Dirichlet_Nitsche.hpp   
    SP/cl_FEM_SP_Ghost_Displacement.hpp
    SP/cl_FEM_SP_Ghost_Virtual_Work.hpp 
    SP/cl_FEM_SP_Stab_Penalty_Contact.hpp
    SP/cl_FEM_SP_Penalty_Contact.hpp 
    SP/cl_FEM_SP_Measure.hpp


    IWG/cl_FEM_IWG_Factory.hpp
    IWG/cl_FEM_IWG.hpp

    IWG/cl_FEM_IWG_FS_Struc_Interface.hpp
    IWG/cl_FEM_IWG_Time_Continuity_Dof.hpp
    IWG/cl_FEM_IWG_Ghost_Normal_Field.hpp
    IWG/cl_FEM_IWG_L2.hpp
    #IWG/cl_FEM_IWG_LSNormal_Bulk.hpp
    #IWG/cl_FEM_IWG_Olsson_CLS_Bulk.hpp
    #IWG/cl_FEM_IWG_Olsson_CLS_Interface.hpp

    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.hpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Fluid_Interface.hpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Bulk.hpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Dirichlet.hpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.hpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.hpp 
    IWG/cl_FEM_IWG_Struc_Stress.hpp

    IWG/cl_FEM_IWG_Advection_Bulk.hpp
    IWG/cl_FEM_IWG_Diffusion_Bulk.hpp
    IWG/cl_FEM_IWG_Diffusion_Dirichlet_Nitsche.hpp
    IWG/cl_FEM_IWG_Diffusion_Neumann.hpp
    IWG/cl_FEM_IWG_Diffusion_Robin.hpp
    IWG/cl_FEM_IWG_Diffusion_Radiation.hpp
    IWG/cl_FEM_IWG_Diffusion_Interface.hpp
    IWG/cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.hpp

    IWG/cl_FEM_IWG_Incompressible_NS_Velocity_Bulk.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_SlipBoundary_Nitsche.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_Neumann.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Velocity_Interface.hpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_Interface.hpp
    
    IWG/cl_FEM_IWG_Compressible_NS_Density_Bulk.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Velocity_Bulk.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Temperature_Bulk.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Advective_Momentum_Flux_Boundary.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Advective_Energy_Flux_Boundary.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Mass_Flux_Neumann.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Traction_Neumann.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Heat_Flux_Neumann.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Temperature_Dirichlet_Nitsche.hpp

    IWG/cl_FEM_IWG_Compressible_NS_Base.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Bulk.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Boundary.hpp
    IWG/cl_FEM_IWG_Compressible_NS_Dirichlet_Nitsche.hpp
    IWG/fn_FEM_IWG_Compressible_NS.hpp
    
    IWG/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp
    IWG/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.hpp
    IWG/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.hpp
    IWG/fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp

    #IWG/cl_FEM_IWG_Helmholtz_Bulk.hpp
    #IWG/cl_FEM_IWG_Helmholtz_Bulk2.hpp
    #IWG/cl_FEM_IWG_Helmholtz_Interface.hpp

    #IWG/cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp
    #IWG/cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp
    IWG/cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp


    IQI/cl_FEM_IQI.hpp
    IQI/cl_FEM_IQI_Factory.hpp
    
    IQI/cl_FEM_IQI_Volume.hpp
    IQI/cl_FEM_IQI_Strain_Energy.hpp
    IQI/cl_FEM_IQI_Stress.hpp
    IQI/cl_FEM_IQI_Dof.hpp
    IQI/cl_FEM_IQI_Max_Dof.hpp
    IQI/cl_FEM_IQI_Max_Stress.hpp
    IQI/cl_FEM_IQI_Property.hpp
    IQI/cl_FEM_IQI_L2_Error_Analytic.hpp
    IQI/cl_FEM_IQI_H1_Error_Analytic.hpp
    IQI/cl_FEM_IQI_H1_Error.hpp
    IQI/cl_FEM_IQI_J_Integral.hpp
    IQI/cl_FEM_IQI_Volume_Fraction.hpp
    IQI/cl_FEM_IQI_2D_Drag_Lift_Coefficient.hpp
    IQI/cl_FEM_IQI_Latent_Heat_Absorption.hpp
    IQI/cl_FEM_IQI_Turbulent_Dynamic_Viscosity.hpp
    IQI/cl_FEM_IQI_Effective_Dynamic_Viscosity.hpp
    IQI/cl_FEM_IQI_Effective_Conductivity.hpp
    IQI/cl_FEM_IQI_Spalart_Allmaras_Coefficient.hpp
    IQI/cl_FEM_IQI_Power_Dissipation.hpp
    IQI/cl_FEM_IQI_Total_Pressure.hpp
    IQI/cl_FEM_IQI_Mass_Flow.hpp
    IQI/cl_FEM_IQI_Thermal_Energy_Convective_Flux.hpp
    IQI/cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.hpp
    IQI/cl_FEM_IQI_Stabilization.hpp
    IQI/cl_FEM_IQI_Homogenized_Constitutive.hpp
    IQI/cl_FEM_IQI_Heat_Method_Penalty.hpp
    IQI/cl_FEM_IQI_Advection_Strong_Residual.hpp
    IQI/cl_FEM_IQI_Strong_Residual_SA.hpp
    IQI/cl_FEM_IQI_Strong_Residual_Incompressible_NS.hpp
    IQI/cl_FEM_IQI_Zienkiewicz_Zhu.hpp
     
    MEA/cl_FEM_Cluster_Measure.hpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Library

# List library source files
set(LIB_SOURCES


    CORE/cl_FEM_Set.cpp
    CORE/cl_FEM_Property.cpp
    CORE/cl_FEM_Model.cpp
    CORE/cl_FEM_Field.cpp 

    ELEM/cl_FEM_Cluster.cpp
    ELEM/cl_FEM_Element_Factory.cpp
    ELEM/cl_FEM_Element_Bulk.cpp
    ELEM/cl_FEM_Element_Sideset.cpp
    ELEM/cl_FEM_Element_Double_Sideset.cpp
    ELEM/cl_FEM_Element_Time_Sideset.cpp
    ELEM/cl_FEM_Element_Time_Boundary.cpp
    ELEM/cl_FEM_Interpolation_Element.cpp

    IP/cl_FEM_Geometry_Interpolator.cpp
    IP/cl_FEM_Field_Interpolator.cpp
    IP/cl_FEM_Field_Interpolator_Manager.cpp

    MM/cl_FEM_Material_Model.cpp
    MM/cl_FEM_MM_Factory.cpp
    MM/cl_FEM_Material_Model_Density_Functions.cpp
    MM/cl_FEM_Material_Model_Pressure_Functions.cpp
    MM/cl_FEM_Material_Model_Temperature_Functions.cpp
    MM/cl_FEM_Material_Model_Thermodynamic_Quantities.cpp
    MM/cl_FEM_Material_Model_FD_Functions.cpp

    MM/cl_FEM_MM_Perfect_Gas.cpp
    #MM/cl_FEM_MM_Van_Der_Waals_Fluid.cpp
    
    CM/cl_FEM_Constitutive_Model.cpp
    CM/cl_FEM_CM_Factory.cpp

    CM/cl_FEM_CM_Diffusion_Linear_Isotropic.cpp
    CM/cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.cpp
    CM/cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.cpp

    CM/cl_FEM_CM_Struc_Nonlinear_Isotropic.cpp
    CM/cl_FEM_CM_Struc_Linear_Isotropic.cpp
     CM/cl_FEM_CM_Struc_Linear_MoriTanaka.cpp
    CM/cl_FEM_CM_Struc_Linear.cpp

    CM/cl_FEM_CM_Fluid_Incompressible.cpp
    CM/cl_FEM_CM_Fluid_Turbulence.cpp
    CM/cl_FEM_CM_Spalart_Allmaras_Turbulence.cpp

    CM/cl_FEM_CM_Fluid_Compressible_Ideal.cpp
    CM/cl_FEM_CM_Fluid_Compressible_Ideal_Test_Tractions.cpp
    CM/cl_FEM_CM_Fluid_Compressible_Van_der_Waals.cpp
    CM/cl_FEM_CM_Fluid_Compressible_Van_der_Waals_Test_Tractions.cpp
    CM/cl_FEM_CM_Compressible_Newtonian_Fluid.cpp
    CM/cl_FEM_CM_Compressible_Newtonian_Fluid_Test_Tractions.cpp  

    SP/cl_FEM_Stabilization_Parameter.cpp
    SP/cl_FEM_SP_Factory.cpp

    SP/cl_FEM_SP_Dirichlet_Nitsche.cpp
    SP/cl_FEM_SP_Nitsche_Interface.cpp
    SP/cl_FEM_SP_Reciprocal_Total_Volume.cpp
    SP/cl_FEM_SP_Convective_Ghost.cpp
    SP/cl_FEM_SP_Time_Velocity_Ghost.cpp
    SP/cl_FEM_SP_SUPG_Advection.cpp
    SP/cl_FEM_SP_YZBeta_Advection.cpp
    SP/cl_FEM_SP_GGLS_Diffusion.cpp
    SP/cl_FEM_SP_Velocity_Dirichlet_Nitsche.cpp
    SP/cl_FEM_SP_Velocity_SlipBoundary_Nitsche.cpp
    SP/cl_FEM_SP_Pressure_Ghost.cpp
    SP/cl_FEM_SP_Incompressible_Flow.cpp
    SP/cl_FEM_SP_Viscous_Ghost.cpp
    SP/cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.cpp
    SP/cl_FEM_SP_Compressible_Dirichlet_Nitsche.cpp
    SP/cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.cpp   
    SP/cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.cpp  
    SP/cl_FEM_SP_Turbulence_Dirichlet_Nitsche.cpp   
    SP/cl_FEM_SP_Ghost_Displacement.cpp
    SP/cl_FEM_SP_Ghost_Virtual_Work.cpp 
    SP/cl_FEM_SP_Stab_Penalty_Contact.cpp
    SP/cl_FEM_SP_Penalty_Contact.cpp  
    SP/cl_FEM_SP_Measure.cpp


    IWG/cl_FEM_IWG.cpp
    IWG/cl_FEM_IWG_Factory.cpp

    IWG/cl_FEM_IWG_L2.cpp
    IWG/cl_FEM_IWG_Time_Continuity_Dof.cpp
    IWG/cl_FEM_IWG_Ghost_Normal_Field.cpp
    #IWG/cl_FEM_IWG_LSNormal_Bulk.cpp
    #IWG/cl_FEM_IWG_Olsson_CLS_Bulk.cpp
    #IWG/cl_FEM_IWG_Olsson_CLS_Interface.cpp

    IWG/cl_FEM_IWG_FS_Struc_Interface.cpp

    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Interface.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Fluid_Interface.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Bulk.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Dirichlet.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.cpp
    IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.cpp
    IWG/cl_FEM_IWG_Struc_Stress.cpp
    
    IWG/cl_FEM_IWG_Advection_Bulk.cpp
    IWG/cl_FEM_IWG_Diffusion_Bulk.cpp
    IWG/cl_FEM_IWG_Diffusion_Dirichlet_Nitsche.cpp
    IWG/cl_FEM_IWG_Diffusion_Neumann.cpp
    IWG/cl_FEM_IWG_Diffusion_Robin.cpp
    IWG/cl_FEM_IWG_Diffusion_Radiation.cpp
    IWG/cl_FEM_IWG_Diffusion_Interface.cpp
    IWG/cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.cpp
    
    IWG/cl_FEM_IWG_Incompressible_NS_Velocity_Bulk.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_Bulk.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Convective_Velocity_Ghost.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_Dirichlet_Nitsche.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Velocity_SlipBoundary_Nitsche.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_SlipBoundary_Nitsche.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_Neumann.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Velocity_Interface.cpp
    IWG/cl_FEM_IWG_Incompressible_NS_Pressure_Interface.cpp
    
    IWG/cl_FEM_IWG_Compressible_NS_Density_Bulk.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Velocity_Bulk.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Temperature_Bulk.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Advective_Momentum_Flux_Boundary.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Advective_Energy_Flux_Boundary.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Mass_Flux_Neumann.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Traction_Neumann.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Heat_Flux_Neumann.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Temperature_Dirichlet_Nitsche.cpp

    IWG/cl_FEM_IWG_Compressible_NS_Base.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Bulk.cpp
    #IWG/cl_FEM_IWG_Compressible_NS_Bulk_Flux_Matrices.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Bulk_Stabilization.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Boundary.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Dirichlet_Nitsche.cpp
    IWG/cl_FEM_IWG_Compressible_NS_Dirichlet_Nitsche_Flux_Matrices.cpp
    IWG/fn_FEM_IWG_Compressible_NS_Dof_Check.cpp
    IWG/fn_FEM_IWG_Compressible_NS_Tools.cpp
    IWG/fn_FEM_IWG_Compressible_NS_Flux_Matrices.cpp
    IWG/fn_FEM_IWG_Compressible_NS_Flux_Matrices_Dof_Derivs.cpp
    IWG/fn_FEM_IWG_Compressible_NS_Flux_Matrices_Var_Derivs.cpp
    
    IWG/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.cpp
    IWG/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Dirichlet.cpp
    IWG/cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.cpp
    IWG/fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.cpp

    #IWG/cl_FEM_IWG_Helmholtz_Bulk.cpp
    #IWG/cl_FEM_IWG_Helmholtz_Bulk2.cpp
    #IWG/cl_FEM_IWG_Helmholtz_Interface.cpp

    #IWG/cl_FEM_IWG_Hamilton_Jacobi_Bulk.cpp
    #IWG/cl_FEM_IWG_Hamilton_Jacobi_Bulk2.cpp
    IWG/cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.cpp
    
    IQI/cl_FEM_IQI.cpp
    IQI/cl_FEM_IQI_Factory.cpp
    IQI/cl_FEM_IQI_Volume.cpp
    IQI/cl_FEM_IQI_Strain_Energy.cpp
    IQI/cl_FEM_IQI_Stress.cpp
    IQI/cl_FEM_IQI_Dof.cpp
    IQI/cl_FEM_IQI_Max_Dof.cpp
    IQI/cl_FEM_IQI_Max_Stress.cpp
    IQI/cl_FEM_IQI_Property.cpp
    IQI/cl_FEM_IQI_L2_Error_Analytic.cpp
    IQI/cl_FEM_IQI_H1_Error_Analytic.cpp
    IQI/cl_FEM_IQI_H1_Error.cpp
    IQI/cl_FEM_IQI_J_Integral.cpp
    IQI/cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp
    IQI/cl_FEM_IQI_Volume_Fraction.cpp
    IQI/cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp
    IQI/cl_FEM_IQI_Latent_Heat_Absorption.cpp
    IQI/cl_FEM_IQI_Turbulent_Dynamic_Viscosity.cpp
    IQI/cl_FEM_IQI_Effective_Dynamic_Viscosity.cpp
    IQI/cl_FEM_IQI_Effective_Conductivity.cpp
    IQI/cl_FEM_IQI_Spalart_Allmaras_Coefficient.cpp
    IQI/cl_FEM_IQI_Power_Dissipation.cpp
    IQI/cl_FEM_IQI_Total_Pressure.cpp
    IQI/cl_FEM_IQI_Mass_Flow.cpp
    IQI/cl_FEM_IQI_Thermal_Energy_Convective_Flux.cpp
    IQI/cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.cpp
    IQI/cl_FEM_IQI_Stabilization.cpp
    IQI/cl_FEM_IQI_Homogenized_Constitutive.cpp
    IQI/cl_FEM_IQI_Heat_Method_Penalty.cpp
    IQI/cl_FEM_IQI_Advection_Strong_Residual.cpp
    IQI/cl_FEM_IQI_Strong_Residual_SA.cpp
    IQI/cl_FEM_IQI_Strong_Residual_Incompressible_NS.cpp
    IQI/cl_FEM_IQI_Zienkiewicz_Zhu.cpp
    
    MEA/cl_FEM_Cluster_Measure.cpp
    
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

set(MTK_INCLUDES
	"IG/"
	"IP/"
	)

# Include depedency directories
foreach(MTK_INCLUDE ${MTK_INCLUDES})
	include_directories(${MORIS_PACKAGE_DIR}/MTK/src/${MTK_INCLUDE})
endforeach()

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
