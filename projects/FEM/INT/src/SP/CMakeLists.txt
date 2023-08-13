# List header files
set(HEADERS
cl_FEM_Stabilization_Parameter.hpp
cl_FEM_SP_Factory.hpp

cl_FEM_SP_Dirichlet_Nitsche.hpp
cl_FEM_SP_Robin_Nitsche.hpp
cl_FEM_SP_Nitsche_Interface.hpp
cl_FEM_SP_Reciprocal_Total_Volume.hpp
cl_FEM_SP_Time_Velocity_Ghost.hpp 
cl_FEM_SP_Convective_Ghost.hpp  
cl_FEM_SP_SUPG_Advection.hpp
cl_FEM_SP_YZBeta_Advection.hpp
cl_FEM_SP_GGLS_Diffusion.hpp
cl_FEM_SP_Velocity_Dirichlet_Nitsche.hpp
cl_FEM_SP_Velocity_SlipBoundary_Nitsche.hpp
cl_FEM_SP_Pressure_Ghost.hpp
cl_FEM_SP_Incompressible_Flow.hpp
cl_FEM_SP_Viscous_Ghost.hpp
cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.hpp
cl_FEM_SP_Compressible_Dirichlet_Nitsche.hpp
cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.hpp   
cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.hpp  
cl_FEM_SP_Turbulence_Dirichlet_Nitsche.hpp   
cl_FEM_SP_Ghost_Displacement.hpp
cl_FEM_SP_Ghost_Virtual_Work.hpp 
cl_FEM_SP_Stab_Penalty_Contact.hpp
cl_FEM_SP_Penalty_Contact.hpp 
cl_FEM_SP_Measure.hpp
cl_FEM_SP_Lagrange_Multiplier_L2.hpp
    )

# List library source files
set(LIB_SOURCES

cl_FEM_Stabilization_Parameter.cpp
cl_FEM_SP_Factory.cpp

cl_FEM_SP_Dirichlet_Nitsche.cpp
cl_FEM_SP_Robin_Nitsche.cpp
cl_FEM_SP_Nitsche_Interface.cpp
cl_FEM_SP_Reciprocal_Total_Volume.cpp
cl_FEM_SP_Convective_Ghost.cpp
cl_FEM_SP_Time_Velocity_Ghost.cpp
cl_FEM_SP_SUPG_Advection.cpp
cl_FEM_SP_YZBeta_Advection.cpp
cl_FEM_SP_GGLS_Diffusion.cpp
cl_FEM_SP_Velocity_Dirichlet_Nitsche.cpp
cl_FEM_SP_Velocity_SlipBoundary_Nitsche.cpp
cl_FEM_SP_Pressure_Ghost.cpp
cl_FEM_SP_Incompressible_Flow.cpp
cl_FEM_SP_Viscous_Ghost.cpp
cl_FEM_SP_Compressible_Velocity_Dirichlet_Nitsche.cpp
cl_FEM_SP_Compressible_Dirichlet_Nitsche.cpp
cl_FEM_SP_SUPG_Spalart_Allmaras_Turbulence.cpp   
cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.cpp  
cl_FEM_SP_Turbulence_Dirichlet_Nitsche.cpp   
cl_FEM_SP_Ghost_Displacement.cpp
cl_FEM_SP_Ghost_Virtual_Work.cpp 
cl_FEM_SP_Stab_Penalty_Contact.cpp
cl_FEM_SP_Penalty_Contact.cpp  
cl_FEM_SP_Measure.cpp
cl_FEM_SP_Lagrange_Multiplier_L2.cpp
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

