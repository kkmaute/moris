# List header files
set(HEADERS
cl_FEM_Cluster.hpp
cl_FEM_Element_Factory.hpp
cl_FEM_ElementProxy.hpp
cl_FEM_Element.hpp

cl_FEM_Element_Bulk.hpp
cl_FEM_Element_Sideset.hpp
cl_FEM_Element_Double_Sideset.hpp
cl_FEM_Element_Time_Sideset.hpp
cl_FEM_Element_Time_Boundary.hpp
    )

# List library source files
set(LIB_SOURCES
cl_FEM_Cluster.cpp
cl_FEM_Element_Factory.cpp
cl_FEM_Element_Bulk.cpp
cl_FEM_Element_Sideset.cpp
cl_FEM_Element_Double_Sideset.cpp
cl_FEM_Element_Time_Sideset.cpp
cl_FEM_Element_Time_Boundary.cpp
cl_FEM_Interpolation_Element.cpp
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

