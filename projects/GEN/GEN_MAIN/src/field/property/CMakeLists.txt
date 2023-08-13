# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Library
set(HEADERS
cl_GEN_Property.hpp
cl_GEN_Constant_Property.hpp
cl_GEN_Scaled_Field.hpp
cl_GEN_User_Defined_Property.hpp
cl_GEN_BSpline_Property.hpp
fn_GEN_create_properties.hpp
st_GEN_Property_Parameters.hpp
)

# List library source files
set(LIB_SOURCES
cl_GEN_Property.cpp
cl_GEN_Constant_Property.cpp
cl_GEN_Scaled_Field.cpp
cl_GEN_User_Defined_Property.cpp
cl_GEN_BSpline_Property.cpp
fn_GEN_create_properties.cpp
)

include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)
