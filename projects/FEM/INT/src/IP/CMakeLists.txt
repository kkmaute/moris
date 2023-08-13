# List header files
set(HEADERS
cl_FEM_Geometry_Interpolator.hpp
cl_FEM_Field_Interpolator.hpp
cl_FEM_Field_Interpolator_Manager.hpp
    )

# List library source files
set(LIB_SOURCES
cl_FEM_Geometry_Interpolator.cpp
cl_FEM_Field_Interpolator.cpp
cl_FEM_Field_Interpolator_Manager.cpp
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

