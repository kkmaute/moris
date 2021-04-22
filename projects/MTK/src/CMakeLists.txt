# MTK Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MTK}_VERSION ${MORIS_VERSION})

# List header dependencies                                                                                                                                                                                                                                                     
set(MTK_INCLUDES
    "IP"
    "IG")

set(MTK_INTERNAL_INCLUDES)

# Include depedency directories
foreach(MTK_INCLUDES ${MTK_INCLUDES})
    list(APPEND MTK_INTERNAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/${MTK_INCLUDES}")
	include_directories(${MORIS_PACKAGE_DIR}/MTK/src/${MTK_INCLUDES})
endforeach()

# create global list for .../MTK/src/ includes
set_property(GLOBAL PROPERTY MTK_SRC_LIST "${MTK_INTERNAL_INCLUDES}")

# List header files
set(HEADERS
cl_Mesh_Enums.hpp
cl_MTK_Mesh_Factory.hpp
cl_MTK_Block_Sets_Info.hpp
cl_MTK_Block.hpp
cl_MTK_Cell_Cluster_Input.hpp
cl_MTK_Cell_Cluster_Proxy.hpp
cl_MTK_Cell_Cluster.hpp
cl_MTK_Cell_Proxy.hpp
cl_MTK_Cell.hpp
cl_MTK_Cell_Info_Factory.hpp
cl_MTK_Cluster.hpp
cl_MTK_Double_Side_Cluster_Input.hpp
cl_MTK_Double_Side_Cluster.hpp
cl_MTK_Double_Side_Set.hpp
cl_MTK_Edge.hpp
cl_MTK_Enums.hpp
cl_MTK_Exodus_IO_Helper.hpp
cl_MTK_Facet.hpp
cl_MTK_Field.hpp
cl_MTK_Field_Analytic.hpp
cl_MTK_Field_Discrete.hpp
cl_MTK_Fields_Info.hpp
cl_MTK_Cell_Info_Hex8.hpp
cl_MTK_Cell_Info_Hex27.hpp
cl_MTK_Cell_Info_Hex64.hpp
cl_MTK_Integration_Mesh.hpp
cl_MTK_Interpolation_Mesh.hpp
cl_MTK_Matrix_Field_Info.hpp
cl_MTK_Mesh_Core.hpp
cl_MTK_Mesh_Data_Input.hpp
cl_MTK_Mesh_Manager.hpp
cl_MTK_Mesh_Checker.hpp
cl_MTK_Mesh_Tools.hpp
cl_MTK_Mesh.hpp
cl_MTK_Node_Sets_Info.hpp
cl_MTK_Periodic_Boundary_Condition_Helper.hpp
cl_MTK_Scalar_Field_Info.hpp
cl_MTK_Set.hpp
cl_MTK_Sets_Info.hpp
cl_MTK_Side_Cluster_Input.hpp
cl_MTK_Side_Cluster.hpp
cl_MTK_Side_Set.hpp
cl_MTK_Side_Sets_Info.hpp
cl_MTK_Cell_Info.hpp
cl_MTK_Cell_Info_Tet4.hpp
cl_MTK_Cell_Info_Quad9.hpp
cl_MTK_Cell_Info_Quad16.hpp
cl_MTK_Topology.hpp
cl_MTK_Vertex_Interpolation.hpp
cl_MTK_Vertex_Proxy.hpp
cl_MTK_Vertex.hpp
cl_MTK_Writer_Exodus.hpp
cl_MTK_Reader_Exodus.hpp
MTK_Tools.hpp

stk_impl/cl_MTK_Cell_STK.hpp
stk_impl/cl_MTK_Integration_Mesh_STK.hpp
stk_impl/cl_MTK_Interpolation_Mesh_STK.hpp
stk_impl/cl_MTK_Mesh_Core_STK.hpp
stk_impl/cl_MTK_Mesh_Data_STK.hpp
stk_impl/cl_MTK_Vertex_Interpolation_STK.hpp
stk_impl/cl_MTK_Vertex_STK.hpp

IP/cl_MTK_Interpolation_Function_Base.hpp
IP/cl_MTK_Interpolation_Function_Factory.hpp
IP/cl_MTK_Interpolation_Function.hpp
IP/cl_MTK_Interpolation_Rule.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Bar1.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Bar2.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Bar3.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Bar4.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Hex27.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Hex64.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Hex8.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Quad16.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Quad4.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Quad8.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Quad9.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Tri3.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Tri6.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Tri10.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Tet4.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Tet10.hpp
IP/cl_MTK_Interpolation_Function_Lagrange_Tet20.hpp
IP/cl_MTK_Interpolation_Function_Constant_Bar2.hpp
IP/cl_MTK_Interpolation_Function_Constant_Point.hpp
IP/cl_MTK_Space_Interpolator.hpp

IG/cl_MTK_Integration_Coeffs.hpp  
IG/cl_MTK_Integration_Rule.hpp
IG/cl_MTK_Integrator.hpp
IG/cl_MTK_Integrator_Test_Polynomial.hpp
IG/cl_MTK_Integration_Coeffs_Bar_1.hpp
IG/cl_MTK_Integration_Coeffs_Bar_2.hpp
IG/cl_MTK_Integration_Coeffs_Bar_3.hpp
IG/cl_MTK_Integration_Coeffs_Bar_4.hpp
IG/cl_MTK_Integration_Coeffs_Bar_5.hpp
IG/cl_MTK_Integration_Coeffs_Bar_6.hpp
IG/cl_MTK_Integration_Coeffs_Base.hpp
IG/cl_MTK_Integration_Coeffs_Hex_2x2x2.hpp
IG/cl_MTK_Integration_Coeffs_Hex_3x3x3.hpp
IG/cl_MTK_Integration_Coeffs_Hex_4x4x4.hpp
IG/cl_MTK_Integration_Coeffs_Hex_5x5x5.hpp
IG/cl_MTK_Integration_Coeffs_Quad_2x2.hpp
IG/cl_MTK_Integration_Coeffs_Quad_3x3.hpp
IG/cl_MTK_Integration_Coeffs_Quad_4x4.hpp
IG/cl_MTK_Integration_Coeffs_Quad_5x5.hpp
IG/cl_MTK_Integration_Coeffs_Tri_1.hpp
IG/cl_MTK_Integration_Coeffs_Tri_3.hpp
IG/cl_MTK_Integration_Coeffs_Tri_4.hpp
IG/cl_MTK_Integration_Coeffs_Tri_6.hpp
IG/cl_MTK_Integration_Coeffs_Tri_7.hpp
IG/cl_MTK_Integration_Coeffs_Tri_12.hpp
IG/cl_MTK_Integration_Coeffs_Tri_13.hpp
IG/cl_MTK_Integration_Coeffs_Tri_16.hpp
IG/cl_MTK_Integration_Coeffs_Tri_19.hpp
IG/cl_MTK_Integration_Coeffs_Tri_25.hpp
IG/cl_MTK_Integration_Coeffs_Tet_1.hpp
IG/cl_MTK_Integration_Coeffs_Tet_4.hpp
IG/cl_MTK_Integration_Coeffs_Tet_5.hpp
IG/cl_MTK_Integration_Coeffs_Tet_11.hpp
IG/cl_MTK_Integration_Coeffs_Tet_15.hpp
IG/cl_MTK_Integration_Coeffs_Tet_20.hpp
IG/cl_MTK_Integration_Coeffs_Tet_35.hpp
IG/cl_MTK_Integration_Coeffs_Tet_56.hpp)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Library

# List library source files
set( LIB_SOURCES
cl_MTK_Mesh_Factory.cpp
stk_impl/cl_MTK_Mesh_Core_STK.cpp
stk_impl/cl_MTK_Integration_Mesh_STK.cpp
stk_impl/cl_MTK_Interpolation_Mesh_STK.cpp
stk_impl/cl_MTK_Vertex_Interpolation_STK.cpp
stk_impl/cl_MTK_Cell_Cluster_STK.cpp
stk_impl/cl_MTK_Side_Cluster_STK.cpp
cl_MTK_Integration_Mesh.cpp
cl_MTK_Block.cpp
cl_MTK_Mesh_Checker.cpp
cl_MTK_Exodus_IO_Helper.cpp
cl_MTK_Field.cpp
cl_MTK_Field_Analytic.cpp
cl_MTK_Field_Discrete.cpp
cl_MTK_Writer_Exodus.cpp
cl_MTK_Reader_Exodus.cpp
cl_MTK_Cell_Info_Factory.cpp
cl_MTK_Block.cpp
cl_MTK_Cell.cpp
cl_MTK_Cell_Info.cpp
cl_MTK_Cell_Info_Hex8.cpp
cl_MTK_Cell_Info_Hex27.cpp
cl_MTK_Cell_Info_Hex64.cpp
cl_MTK_Cell_Info_Tet4.cpp
cl_MTK_Cell_Info_Tri3.cpp
cl_MTK_Cell_Info_Quad4.cpp
cl_MTK_Cell_Info_Quad9.cpp
cl_MTK_Cell_Info_Quad16.cpp
cl_MTK_Mesh_Core.cpp
cl_MTK_Mesh_Manager.cpp
cl_MTK_Periodic_Boundary_Condition_Helper.cpp
cl_MTK_Double_Side_Cluster.cpp
cl_MTK_Side_Cluster.cpp
    cl_MTK_Mesh_Pair.cpp

IP/cl_MTK_Interpolation_Function_Factory.cpp
IP/cl_MTK_Interpolation_Rule.cpp
IP/cl_MTK_Space_Interpolator.cpp

IG/cl_MTK_Integration_Rule.cpp
IG/cl_MTK_Integrator.cpp
IG/cl_MTK_Integrator_Test_Polynomial.cpp)

# List library dependencies
set(LIB_DEPENDENCIES
	${COM}-lib
	${LINALG}-lib
	${HMR}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MTK}-lib STATIC ${LIB_SOURCES})
target_include_directories(${MTK}-lib PUBLIC
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${MTK_INTERNAL_INCLUDES};>"
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/stk_impl>
	$<INSTALL_INTERFACE:${${MTK}_HEADER_INSTALL_DIR}> )
target_link_libraries(${MTK}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MTK}-lib PROPERTIES OUTPUT_NAME ${MTK})

foreach(TPL ${MTK_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${MTK}-lib
    	PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MTK}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${MTK}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${MTK_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND MTK_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND MTK_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(MTK_COMPILE_DEFINITIONS ${MTK}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${MTK}Config.cmake.in
    ${${MTK}_CONFIG_BUILD_DIR}/${MTK}Config.cmake
    INSTALL_DESTINATION 
    ${${MTK}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${MTK}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${MTK}_CONFIG_BUILD_DIR}/${MTK}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Installation

# Install *Config and *ConfigVersion files
install(FILES ${${MTK}_CONFIG_BUILD_DIR}/${MTK}Config.cmake 
    ${${MTK}_CONFIG_BUILD_DIR}/${MTK}ConfigVersion.cmake
    DESTINATION ${${MTK}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${MTK}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${MTK}-lib
    DESTINATION ${${MTK}_LIB_INSTALL_DIR}
    EXPORT ${MTK}Targets )

# Export library into *Targets file
install(EXPORT ${MTK}Targets
    DESTINATION ${${MTK}_CONFIG_INSTALL_DIR} )
