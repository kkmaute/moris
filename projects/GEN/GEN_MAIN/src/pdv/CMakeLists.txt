# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Library
set(HEADERS
cl_GEN_Pdv_Host_Manager.hpp
cl_GEN_Child_Node.hpp
cl_GEN_Interpolation_Pdv_Host.hpp
cl_GEN_Pdv.hpp
cl_GEN_Pdv_Value.hpp
cl_GEN_Pdv_Property.hpp
cl_GEN_Intersection_Node.hpp
cl_GEN_Intersection_Node_Linear.hpp
cl_GEN_Intersection_Node_Bilinear.hpp
)

# List library source files
set(LIB_SOURCES
cl_GEN_Pdv_Host_Manager.cpp
cl_GEN_Interpolation_Pdv_Host.cpp
cl_GEN_Pdv.cpp
cl_GEN_Pdv_Value.cpp
cl_GEN_Pdv_Property.cpp
cl_GEN_Child_Node.cpp
cl_GEN_Intersection_Node.cpp
cl_GEN_Intersection_Node_Linear.cpp
cl_GEN_Intersection_Node_Bilinear.cpp
)

include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)
