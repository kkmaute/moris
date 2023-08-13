# List header files
set(HEADERS
    cl_FEM_IQI.hpp
    cl_FEM_IQI_Factory.hpp
    
    cl_FEM_IQI_Volume.hpp
    cl_FEM_IQI_Strain_Energy.hpp
    cl_FEM_IQI_Stress.hpp
    cl_FEM_IQI_Dof.hpp
    cl_FEM_IQI_Eigen_Vector.hpp
    cl_FEM_IQI_ALM_Dof.hpp
    cl_FEM_IQI_Max_Dof.hpp
    cl_FEM_IQI_Max_Stress.hpp
    cl_FEM_IQI_Property.hpp
    cl_FEM_IQI_L2_Error_Analytic.hpp
    cl_FEM_IQI_H1_Error_Analytic.hpp
    cl_FEM_IQI_H1_Error.hpp
    cl_FEM_IQI_J_Integral.hpp
    cl_FEM_IQI_Volume_Fraction.hpp
    cl_FEM_IQI_2D_Drag_Lift_Coefficient.hpp
    cl_FEM_IQI_Latent_Heat_Absorption.hpp
    cl_FEM_IQI_Turbulent_Dynamic_Viscosity.hpp
    cl_FEM_IQI_Effective_Dynamic_Viscosity.hpp
    cl_FEM_IQI_Effective_Conductivity.hpp
    cl_FEM_IQI_Spalart_Allmaras_Coefficient.hpp
    cl_FEM_IQI_Power_Dissipation.hpp
    cl_FEM_IQI_Total_Pressure.hpp
    cl_FEM_IQI_Mass_Flow.hpp
    cl_FEM_IQI_Thermal_Energy_Convective_Flux.hpp
    cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.hpp
    cl_FEM_IQI_Stabilization.hpp
    cl_FEM_IQI_Homogenized_Constitutive.hpp
    cl_FEM_IQI_Heat_Method_Penalty.hpp
    cl_FEM_IQI_Advection_Strong_Residual.hpp
    cl_FEM_IQI_Strong_Residual_SA.hpp
    cl_FEM_IQI_Strong_Residual_Incompressible_NS.hpp
    cl_FEM_IQI_Zienkiewicz_Zhu.hpp
    cl_FEM_IQI_Jump_Dof.hpp
    cl_FEM_IQI_Jump_Traction.hpp
    )

# List library source files
set(LIB_SOURCES
cl_FEM_IQI.cpp
cl_FEM_IQI_Factory.cpp
cl_FEM_IQI_Volume.cpp
cl_FEM_IQI_Strain_Energy.cpp
cl_FEM_IQI_Stress.cpp
cl_FEM_IQI_Dof.cpp
cl_FEM_IQI_Eigen_Vector.cpp
cl_FEM_IQI_ALM_Dof.cpp
cl_FEM_IQI_Max_Dof.cpp
cl_FEM_IQI_Max_Stress.cpp
cl_FEM_IQI_Property.cpp
cl_FEM_IQI_L2_Error_Analytic.cpp
cl_FEM_IQI_H1_Error_Analytic.cpp
cl_FEM_IQI_H1_Error.cpp
cl_FEM_IQI_J_Integral.cpp
cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp
cl_FEM_IQI_Volume_Fraction.cpp
cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp
cl_FEM_IQI_Latent_Heat_Absorption.cpp
cl_FEM_IQI_Turbulent_Dynamic_Viscosity.cpp
cl_FEM_IQI_Effective_Dynamic_Viscosity.cpp
cl_FEM_IQI_Effective_Conductivity.cpp
cl_FEM_IQI_Spalart_Allmaras_Coefficient.cpp
cl_FEM_IQI_Power_Dissipation.cpp
cl_FEM_IQI_Total_Pressure.cpp
cl_FEM_IQI_Mass_Flow.cpp
cl_FEM_IQI_Thermal_Energy_Convective_Flux.cpp
cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.cpp
cl_FEM_IQI_Stabilization.cpp
cl_FEM_IQI_Homogenized_Constitutive.cpp
cl_FEM_IQI_Heat_Method_Penalty.cpp
cl_FEM_IQI_Advection_Strong_Residual.cpp
cl_FEM_IQI_Strong_Residual_SA.cpp
cl_FEM_IQI_Strong_Residual_Incompressible_NS.cpp
cl_FEM_IQI_Zienkiewicz_Zhu.cpp
cl_FEM_IQI_Jump_Dof.cpp
cl_FEM_IQI_Jump_Traction.cpp
    )


# Create library

#include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)
list(TRANSFORM LIB_SOURCES PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
list(TRANSFORM HEADERS PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

#target_sources(morisCombinedLib PUBLIC ${LIB_SOURCES} ${HEADERS})
target_sources(morisCombinedLib PUBLIC ${HEADERS})

target_include_directories(morisCombinedLib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include> )

