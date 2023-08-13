# List header files
set(HEADERS
cl_MTK_Cell_STK.hpp
cl_MTK_Integration_Mesh_STK.hpp
cl_MTK_Interpolation_Mesh_STK.hpp
cl_MTK_Mesh_Core_STK.hpp
cl_MTK_Mesh_Data_STK.hpp
cl_MTK_Vertex_Interpolation_STK.hpp
cl_MTK_Vertex_STK.hpp
cl_MTK_Cell_Cluster_STK.hpp
cl_MTK_Side_Cluster_STK.hpp
)

# List library source files
set( LIB_SOURCES
cl_MTK_Mesh_Core_STK.cpp
cl_MTK_Integration_Mesh_STK.cpp
cl_MTK_Interpolation_Mesh_STK.cpp
cl_MTK_Vertex_Interpolation_STK.cpp
cl_MTK_Cell_Cluster_STK.cpp
cl_MTK_Side_Cluster_STK.cpp
)


include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)
