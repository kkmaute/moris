# XTK Source --------------------------------------------------------------                                                                                                                                                                                                    
# -------------------------------------------------------------------------                                                                                                                                                                                                    

# Set version                                                                                                                                                                                                                                                                  
set(${XTK}_VERSION ${MORIS_VERSION})                                                                                                                                                                                                                                           

# List header files                                                                                                                                                                                                                                                            
set(HEADERS
    fn_num_digits.hpp
    cl_XTK_Matrix_Base_Utilities.hpp
    fn_verify_tet_topology.hpp
    cl_Interpolaton.hpp
    fn_approximate.hpp
    fn_bubble_sort.hpp
    fn_hex_8_volume.hpp
    fn_matrix_to_offset.hpp
    fn_Pairing.hpp
    fn_tet_volume.hpp
    cl_XTK_Basis_Function.hpp
    cl_XTK_Topology.hpp
    cl_XTK_Edge_Topology.hpp
    cl_XTK_Hexahedron_8_Basis_Function.hpp
    cl_XTK_Hexahedron_8_Topology.hpp
    cl_XTK_Linear_Basis_Functions.hpp
    cl_XTK_Quad_4_Basis_Function.hpp
    cl_XTK_Quad_4_Topology.hpp
    cl_XTK_Cell_CM.hpp
   	cl_XTK_Cell_No_CM.hpp
    cl_MTK_Vertex_Interpolation_XTK_Impl.hpp
    cl_MTK_Vertex_XTK_Impl.hpp
    cl_XTK_Background_Mesh.hpp
    cl_XTK_Cell_Cluster.hpp
    cl_XTK_Child_Mesh_Modification_Template.hpp
    cl_XTK_Child_Mesh.hpp
    cl_XTK_Cut_Mesh.hpp
    cl_XTK_Decomposition_Data.hpp
    cl_XTK_Downward_Inheritance.hpp
    cl_XTK_Enrichment.hpp
    cl_XTK_Entity.hpp
    cl_XTK_Enums.hpp
    cl_XTK_Integration_Mesh_Generator.hpp
    cl_XTK_Decomposition_Algorithm.hpp
    cl_XTK_Regular_Subdivision_Interface.hpp
    cl_XTK_Node_Hierarchy_Interface.hpp
    cl_XTK_Elevate_Order_Interface.hpp
    cl_XTK_Cut_Integration_Mesh.hpp
    cl_XTK_External_Mesh_Data.hpp
    cl_XTK_Face_Registry.hpp
    cl_XTK_Ghost_Stabilization.hpp
    cl_XTK_Hole_Seeder.hpp
    cl_XTK_Input_Generator.hpp
    cl_XTK_Interface_Element.hpp
    cl_XTK_Interpolation_Cell.hpp
    cl_XTK_Interpolation_Cell_Unzipped.hpp
    cl_XTK_Interpolation_Vertex_Unzipped.hpp
    cl_XTK_Enriched_Integration_Mesh.hpp
    cl_XTK_Enriched_Interpolation_Mesh.hpp
    cl_XTK_Mesh_Cleanup.hpp
    cl_XTK_Multigrid.hpp
    cl_XTK_Model.hpp
#    cl_XTK_Node.hpp
    cl_XTK_Output_Options.hpp
    cl_XTK_Paramfile.hpp
    cl_XTK_Pending_Node.hpp
    cl_XTK_Side_Cluster.hpp
    cl_XTK_Vertex_Enrichment.hpp
    cl_XTK_Field.hpp
    fn_compute_interface_surface_area.hpp
    fn_compute_xtk_model_areas.hpp
    fn_compute_xtk_model_volumes.hpp
    fn_create_edges_from_element_to_node.hpp
    fn_create_faces_from_element_to_node.hpp
    fn_generate_cross_boundary_pairs.hpp
    fn_generate_element_to_element.hpp
    fn_generate_shared_face_element_graph.hpp
    fn_local_child_mesh_flood_fill.hpp
    fn_mesh_flood_fill.hpp
    fn_prune_element_to_element.hpp
    fn_write_element_ownership_as_field.hpp
    fn_generate_element_to_element_2D.hpp
    xtk_typedefs.hpp
    cl_XTK_Decomposition_Algorithm_Factory.hpp
    cl_XTK_Diagnostics.hpp
    )                                                                                                                                                                                                                                                                   


include_directories(${MORIS_PACKAGE_DIR}/include)

set(GEN_INCLUDES
	"field"
	"field/geometry"
	"field/property"
	"pdv"
	)

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})
	include_directories(${MORIS_PACKAGE_DIR}/GEN/GEN_MAIN/src/${GEN_INCLUDES})
endforeach()

# Include WRK directory
include_directories(${MORIS_PACKAGE_DIR}/WRK/src)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Library

# List source files
set(LIB_SOURCES
cl_XTK_Cell_CM.cpp
cl_XTK_Cell_No_CM.cpp
cl_MTK_Vertex_Interpolation_XTK_Impl.cpp
cl_MTK_Vertex_XTK_Impl.cpp
cl_XTK_Background_Mesh.cpp
cl_XTK_Child_Mesh.cpp
cl_XTK_Cut_Mesh.cpp
cl_XTK_Enrichment.cpp
cl_XTK_Field.cpp
cl_XTK_Ghost_Stabilization.cpp
cl_XTK_Hole_Seeder.cpp
cl_XTK_Mesh_Cleanup.cpp
cl_XTK_Model.cpp
cl_XTK_Multigrid.cpp
cl_XTK_Paramfile.cpp
cl_XTK_Interpolation_Vertex_Unzipped.cpp
cl_XTK_Interpolation_Cell_Unzipped.cpp
cl_XTK_Enriched_Interpolation_Mesh.cpp
cl_XTK_Vertex_Enrichment.cpp
cl_XTK_Enriched_Integration_Mesh.cpp
cl_XTK_Cell_Cluster.cpp
cl_XTK_Side_Cluster.cpp
cl_XTK_Interpolation_Cell.cpp
cl_XTK_Decomposition_Algorithm.cpp
cl_XTK_Integration_Mesh_Generator.cpp
cl_XTK_Regular_Subdivision_Interface.cpp
cl_XTK_Node_Hierarchy_Interface.cpp
cl_XTK_Node_Hierarchy_Template_Library.cpp
cl_XTK_Elevate_Order_Interface.cpp
cl_XTK_Cut_Integration_Mesh.cpp
cl_XTK_Decomposition_Algorithm_Factory.cpp
cl_XTK_Diagnostics.cpp
cl_XTK_Octree_Interface.cpp)

# List library dependencies
set(LIB_DEPENDENCIES 
	${GEN_MAIN}-lib
	${LINALG}-lib
	${MTK}-lib
	${ALG}-lib
	${CHR}-lib
	${INT}-lib
	${MDL}-lib
    ${MAP}-lib
    ${HMR}-lib
    ${TOL}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${XTK}-lib STATIC ${LIB_SOURCES})
target_include_directories(${XTK}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${XTK}_HEADER_INSTALL_DIR}> )
target_link_libraries(${XTK}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${XTK}-lib PROPERTIES OUTPUT_NAME ${XTK})

# Link third party libraries
foreach(TPL ${XTK_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${XTK}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${XTK}-lib 
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${XTK}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${XTK_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND XTK_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND XTK_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} QUIET HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND XTK_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} REQUIRED HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(XTK_COMPILE_DEFINITIONS ${XTK}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${XTK}Config.cmake.in
    ${${XTK}_CONFIG_BUILD_DIR}/${XTK}Config.cmake
    INSTALL_DESTINATION 
    ${${XTK}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${XTK}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${XTK}_CONFIG_BUILD_DIR}/${XTK}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )
 
 # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 # XTK Installation
 
# Install *Config and *ConfigVersion files
install(FILES ${${XTK}_CONFIG_BUILD_DIR}/${XTK}Config.cmake 
    ${${XTK}_CONFIG_BUILD_DIR}/${XTK}ConfigVersion.cmake
    DESTINATION ${${XTK}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${XTK}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${XTK}-lib
    DESTINATION ${${XTK}_LIB_INSTALL_DIR}
    EXPORT ${XTK}Targets )


# Export executable in to *Targets file
install(EXPORT ${XTK}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${XTK} )
