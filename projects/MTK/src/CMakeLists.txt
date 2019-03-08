# MTK Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MTK}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
cl_Mesh_Enums.hpp
cl_Mesh_Factory.hpp
cl_MTK_Block_Sets_Info.hpp
cl_MTK_Block.cpp
cl_MTK_Block.hpp
cl_MTK_Cell.hpp
cl_MTK_Edge.hpp
cl_MTK_Enums.hpp
cl_MTK_Exodus_IO_Helper.hpp
cl_MTK_Facet_Cluster.hpp
cl_MTK_Facet.hpp
cl_MTK_Fields_Info.hpp
cl_MTK_Matrix_Field_Info.hpp
cl_MTK_Mesh_Data_Input.hpp
cl_MTK_Mesh_Tools.hpp
cl_MTK_Mesh.hpp
cl_MTK_Node_Sets_Info.hpp
cl_MTK_Scalar_Field_Info.hpp
cl_MTK_Sets_Info.hpp
cl_MTK_Side_Sets_Info.hpp
cl_MTK_Topology.hpp
cl_MTK_Vertex_Interpolation.hpp
cl_MTK_Vertex.hpp
cl_MTK_Vertexset.hpp
MTK_Tools.hpp )

# List header dependencies
set(MTK_INCLUDES
    ${LINALG}
    ${COM}
    ${TOL}
    ${ALG}
    )

# Include depedency directories
foreach(MTK_INCLUDE ${MTK_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${MTK_INCLUDE}/src)
endforeach()

# additional includes
include_directories(${MORIS_PACKAGE_DIR}/MTK/src/stk_impl)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Library

# List library source files
set( LIB_SOURCES
    cl_MTK_Block.cpp
    stk_impl/cl_MTK_Mesh_STK.cpp
    stk_impl/cl_MTK_Vertex_Interpolation_STK.cpp)

# List library dependencies
set(LIB_DEPENDENCIES
	${COM}-lib
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MTK}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${MTK}-lib ${LIB_DEPENDENCIES})
set_target_properties(${MTK}-lib PROPERTIES OUTPUT_NAME ${MTK})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Config Files

# Build *Config and *ConfigVersion files
 include(CMakePackageConfigHelpers)
 
 configure_package_config_file(
     ${MTK}Config.cmake.in
     ${${MTK}_CONFIG_BUILD_DIR}/${MTK}Config.cmake
     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MTK}
     PATH_VARS LIB_INSTALL_DIR )
 
 write_basic_package_version_file(
     ${${MTK}_CONFIG_BUILD_DIR}/${MTK}ConfigVersion.cmake
     COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Installation

 # Install *Config and *ConfigVersion files
 install(FILES ${${MTK}_CONFIG_BUILD_DIR}/${MTK}Config.cmake 
     ${${MTK}_CONFIG_BUILD_DIR}/${MTK}ConfigVersion.cmake
     DESTINATION ${CONFIG_INSTALL_DIR}/${MTK} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MTK} )

# Install library and executable
install(TARGETS ${MTK}-lib
    EXPORT ${MTK}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${MTK}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${MTK}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${MTK} )

# Export library in to *Targets file
install(EXPORT ${MTK}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${MTK} )
