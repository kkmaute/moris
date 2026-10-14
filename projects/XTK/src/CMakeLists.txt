# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# XTK Source

# Package version.
set(${XTK}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_XTK_Matrix_Base_Utilities.hpp
    fn_verify_tet_topology.hpp
    cl_Interpolaton.hpp
    fn_approximate.hpp
    fn_bubble_sort.hpp
    fn_hex_8_volume.hpp
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
    cl_XTK_Side_Cluster.hpp
    #cl_XTK_Cluster_Group.hpp
    cl_XTK_Cell_Cluster_Group.hpp
    cl_XTK_Side_Cluster_Group.hpp
    cl_XTK_Child_Mesh_Modification_Template.hpp
    cl_XTK_Child_Mesh.hpp
    cl_XTK_Cut_Mesh.hpp
    cl_XTK_Decomposition_Data.hpp
    cl_XTK_Downward_Inheritance.hpp
    cl_XTK_Enrichment.hpp
    cl_XTK_Entity.hpp
    cl_XTK_Integration_Mesh_Generator.hpp
    cl_XTK_Decomposition_Algorithm.hpp
    cl_XTK_Regular_Subdivision_Interface.hpp
    cl_XTK_Node_Hierarchy_Interface.hpp
    cl_XTK_Elevate_Order_Interface.hpp
    cl_XTK_Octree_Interface.hpp
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
    #cl_XTK_Node.hpp
    cl_XTK_Output_Options.hpp
    cl_XTK_Paramfile.hpp
    cl_XTK_Pending_Node.hpp
    cl_XTK_Vertex_Enrichment.hpp
    cl_XTK_Field.hpp
    fn_compute_xtk_model_areas.hpp
    fn_compute_xtk_model_volumes.hpp
    fn_create_edges_from_element_to_node.hpp
    fn_create_faces_from_element_to_node.hpp
    fn_generate_element_to_element.hpp
    fn_local_child_mesh_flood_fill.hpp
    fn_mesh_flood_fill.hpp
    fn_prune_element_to_element.hpp
    fn_write_element_ownership_as_field.hpp
    fn_generate_element_to_element_2D.hpp
    fn_determine_cell_topology.hpp
    fn_sort_points_by_coordinates.hpp
    xtk_typedefs.hpp
    cl_XTK_Decomposition_Algorithm_Factory.hpp
    cl_XTK_Diagnostics.hpp
    cl_XTK_Subphase_Group.hpp
    cl_XTK_Integration_Mesh_Cleanup.hpp
    fn_XTK_convert_cell_to_map.hpp
    fn_XTK_convert_cell_to_multiset.hpp
    fn_XTK_find_most_frequent_int_in_cell.hpp
    fn_XTK_Multiset_Operations.hpp
    fn_XTK_punch_card_sum.hpp
    fn_XTK_find_in_cell.hpp
    fn_XTK_match_normal_to_side_ordinal.hpp
    cl_XTK_Basis_Processor.hpp
    cl_XTK_HMR_Helper.hpp
)

include_directories(${MORIS_PACKAGE_DIR}/include)

set(GEN_INCLUDES "ADV" "ADV/fields" "PDV")

# Include dependency directories.
foreach(GEN_INCLUDES IN LISTS GEN_INCLUDES)
    include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
endforeach()

# Include WRK directory
include_directories(${MORIS_PACKAGE_DIR}/WRK/src)

# XTK Library

# Sources.
set(LIB_SOURCES
    cl_XTK_Model.cpp
    cl_XTK_Cell_CM.cpp
    cl_XTK_Cell_No_CM.cpp
    cl_XTK_Interpolation_Cell.cpp
    cl_XTK_Cell_Cluster.cpp
    cl_XTK_Side_Cluster.cpp
    #cl_XTK_Cluster_Group.cpp
    cl_XTK_Cell_Cluster_Group.cpp
    cl_XTK_Side_Cluster_Group.cpp
    cl_MTK_Vertex_Interpolation_XTK_Impl.cpp
    cl_MTK_Vertex_XTK_Impl.cpp
    cl_XTK_Background_Mesh.cpp
    cl_XTK_Child_Mesh.cpp
    cl_XTK_Cut_Mesh.cpp
    cl_XTK_Cut_Integration_Mesh.cpp
    cl_XTK_Enriched_Interpolation_Mesh.cpp
    cl_XTK_Enriched_Integration_Mesh.cpp
    cl_XTK_Mesh_Cleanup.cpp
    cl_XTK_Enrichment.cpp
    cl_XTK_Vertex_Enrichment.cpp
    cl_XTK_Field.cpp
    cl_XTK_Ghost_Stabilization.cpp
    cl_XTK_Hole_Seeder.cpp
    cl_XTK_Multigrid.cpp
    cl_XTK_Paramfile.cpp
    cl_XTK_Interpolation_Vertex_Unzipped.cpp
    cl_XTK_Interpolation_Cell_Unzipped.cpp
    cl_XTK_Integration_Mesh_Generator.cpp
    cl_XTK_Decomposition_Algorithm.cpp
    cl_XTK_Decomposition_Algorithm_Factory.cpp
    cl_XTK_Regular_Subdivision_Interface.cpp
    cl_XTK_Node_Hierarchy_Interface.cpp
    cl_XTK_Node_Hierarchy_Template_Library.cpp
    cl_XTK_Elevate_Order_Interface.cpp
    cl_XTK_Octree_Interface.cpp
    cl_XTK_Diagnostics.cpp
    cl_XTK_Integration_Mesh_Cleanup.cpp
    cl_XTK_Basis_Processor.cpp
    cl_XTK_HMR_Helper.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES
    ${GEN}-lib
    ${LINALG}-lib
    ${MTK}-lib
    ${ALG}-lib
    ${CHR}-lib
    ${INT}-lib
    ${MDL}-lib
    ${MAP}-lib
    ${HMR}-lib
    ${TOL}-lib
    ${ENM}-lib
)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${XTK}-lib XTK_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${XTK_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${XTK}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${XTK_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${XTK_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${XTK})

# Third-party dependencies.
foreach(TPL IN LISTS XTK_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${XTK_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${XTK_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${XTK_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    XTK
    TARGET
    ${XTK_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${XTK_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${XTK}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${XTK}-lib DESTINATION ${${XTK}_LIB_INSTALL_DIR} EXPORT ${XTK}Targets)

# Export executable in to *Targets file
install(EXPORT ${XTK}Targets DESTINATION ${CONFIG_INSTALL_DIR}/${XTK})
