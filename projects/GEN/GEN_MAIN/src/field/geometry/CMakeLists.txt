# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Library
set(HEADERS
cl_GEN_Geometry.hpp
cl_GEN_Circle.hpp
cl_GEN_Superellipse.hpp
cl_GEN_Sphere.hpp
cl_GEN_Superellipsoid.hpp
cl_GEN_Plane.hpp
cl_GEN_Mesh_Field_Geometry.hpp
cl_GEN_Voxel_Input.hpp
cl_GEN_Single_Grain.hpp
cl_GEN_User_Defined_Geometry.hpp
cl_GEN_BSpline_Geometry.hpp
cl_GEN_Stored_Geometry.hpp
cl_GEN_Geometry_SDF.hpp
cl_GEN_Image_SDF_Geometry.hpp
cl_GEN_Multigeometry.hpp
cl_GEN_Swiss_Cheese_Slice.hpp
fn_GEN_create_geometries.hpp
st_GEN_Geometry_Parameters.hpp
)

# List library source files
set(LIB_SOURCES
cl_GEN_Geometry.cpp
cl_GEN_Circle.cpp
cl_GEN_Superellipse.cpp
cl_GEN_Sphere.cpp
cl_GEN_Superellipsoid.cpp
cl_GEN_Plane.cpp
cl_GEN_Mesh_Field_Geometry.cpp
cl_GEN_Voxel_Input.cpp
cl_GEN_Single_Grain.cpp
cl_GEN_User_Defined_Geometry.cpp
cl_GEN_BSpline_Geometry.cpp
cl_GEN_Stored_Geometry.cpp
cl_GEN_Geometry_SDF.cpp
cl_GEN_Image_SDF_Geometry.cpp
cl_GEN_Multigeometry.cpp
cl_GEN_Swiss_Cheese_Slice.cpp
fn_GEN_create_geometries.cpp
)

include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)
