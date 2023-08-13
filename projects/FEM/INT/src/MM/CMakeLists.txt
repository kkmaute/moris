# List header files
set(HEADERS
cl_FEM_Material_Model.hpp
cl_FEM_MM_Factory.hpp

cl_FEM_MM_Perfect_Gas.hpp
    )

# List library source files
set(LIB_SOURCES
cl_FEM_Material_Model.cpp
cl_FEM_MM_Factory.cpp
cl_FEM_Material_Model_Density_Functions.cpp
cl_FEM_Material_Model_Pressure_Functions.cpp
cl_FEM_Material_Model_Temperature_Functions.cpp
cl_FEM_Material_Model_Thermodynamic_Quantities.cpp
cl_FEM_Material_Model_FD_Functions.cpp

cl_FEM_MM_Perfect_Gas.cpp
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

