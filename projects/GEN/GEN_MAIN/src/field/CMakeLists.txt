set(HEADERS
cl_GEN_Field.hpp
cl_GEN_Field_Analytic.hpp
cl_GEN_Field_Discrete_Integration.hpp
cl_GEN_Field_Discrete_Interpolation.hpp
cl_GEN_User_Defined_Field.hpp
cl_GEN_BSpline_Field.hpp
)

# List library source files
set(LIB_SOURCES
cl_GEN_Field.cpp
cl_GEN_Field_Analytic.cpp
cl_GEN_Field_Discrete_Integration.cpp
cl_GEN_Field_Discrete_Interpolation.cpp
cl_GEN_User_Defined_Field.cpp
cl_GEN_BSpline_Field.cpp
st_GEN_Field_Parameters.hpp
)


include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)

add_subdirectory(geometry)
add_subdirectory(property)
