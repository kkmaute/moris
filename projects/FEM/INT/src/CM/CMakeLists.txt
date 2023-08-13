# List header files
set(HEADERS
cl_FEM_Constitutive_Model.hpp
cl_FEM_CM_Factory.hpp

cl_FEM_CM_Diffusion_Linear_Isotropic.hpp
cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.hpp
fn_FEM_CM_Phase_State_Functions.hpp
cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp

cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp
cl_FEM_CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff.hpp
cl_FEM_CM_Struc_Nonlinear_Isotropic_Neo_Hookean.hpp
cl_FEM_CM_Struc_Linear_Isotropic.hpp
cl_FEM_CM_Struc_Linear_MoriTanaka.hpp
cl_FEM_CM_Struc_Linear.hpp

cl_FEM_CM_Fluid_Incompressible.hpp
cl_FEM_CM_Fluid_Turbulence.hpp
cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp

cl_FEM_CM_Fluid_Compressible_Ideal.hpp
cl_FEM_CM_Fluid_Compressible_Van_der_Waals.hpp
cl_FEM_CM_Compressible_Newtonian_Fluid.hpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/int Library

# List library source files
set(LIB_SOURCES
cl_FEM_Constitutive_Model.cpp
cl_FEM_CM_Factory.cpp

cl_FEM_CM_Diffusion_Linear_Isotropic.cpp
cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.cpp
cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.cpp

cl_FEM_CM_Struc_Nonlinear_Isotropic.cpp
cl_FEM_CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff.cpp
cl_FEM_CM_Struc_Nonlinear_Isotropic_Neo_Hookean.cpp
cl_FEM_CM_Struc_Linear_Isotropic.cpp
cl_FEM_CM_Struc_Linear_MoriTanaka.cpp
cl_FEM_CM_Struc_Linear.cpp

cl_FEM_CM_Fluid_Incompressible.cpp
cl_FEM_CM_Fluid_Turbulence.cpp
cl_FEM_CM_Spalart_Allmaras_Turbulence.cpp

cl_FEM_CM_Fluid_Compressible_Ideal.cpp
cl_FEM_CM_Fluid_Compressible_Ideal_Test_Tractions.cpp
cl_FEM_CM_Fluid_Compressible_Van_der_Waals.cpp
cl_FEM_CM_Fluid_Compressible_Van_der_Waals_Test_Tractions.cpp
cl_FEM_CM_Compressible_Newtonian_Fluid.cpp
cl_FEM_CM_Compressible_Newtonian_Fluid_Test_Tractions.cpp  
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

