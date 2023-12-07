#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# MTK Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MTK}_VERSION ${MORIS_VERSION})

# List header dependencies                                                                                                                                                                                                                                                     
set(MTK_INCLUDES
        "cluster"
        "contact"
        "field"
        "geometry"
        "geometry/cell"
        "geometry/vertex"
        "integration"
        "interpolation"
        "intersection"
        "io"
        "mesh"
        "set"
)

set(MTK_INTERNAL_INCLUDES)

# Include depedency directories
foreach (MTK_INCLUDES ${MTK_INCLUDES})
    list(APPEND MTK_INTERNAL_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/${MTK_INCLUDES}")
    include_directories(${MORIS_PACKAGE_DIR}/MTK/src/${MTK_INCLUDES})
endforeach ()

# create global list for .../MTK/src/ includes
set_property(GLOBAL PROPERTY MTK_SRC_LIST "${MTK_INTERNAL_INCLUDES}")

# List header files
set(HEADERS
        MTK_Tools.hpp

        cl_MTK_Periodic_Boundary_Condition_Helper.hpp

        cluster/cl_MTK_Cell_Cluster.hpp
        cluster/cl_MTK_Cell_Cluster_DataBase.hpp
        cluster/cl_MTK_Cell_Cluster_Group.hpp
        cluster/cl_MTK_Cell_Cluster_Group_DataBase.hpp
        cluster/cl_MTK_Cell_Cluster_Input.hpp
        cluster/cl_MTK_Cell_Cluster_Proxy.hpp
        cluster/cl_MTK_Cluster.hpp
        cluster/cl_MTK_Cluster_Group.hpp
        cluster/cl_MTK_Double_Side_Cluster.hpp
        cluster/cl_MTK_Double_Side_Cluster_Input.hpp
        cluster/cl_MTK_Nonconformal_Side_Cluster.hpp
        cluster/cl_MTK_Side_Cluster.hpp
        cluster/cl_MTK_Side_Cluster_DataBase.hpp
        cluster/cl_MTK_Side_Cluster_Group.hpp
        cluster/cl_MTK_Side_Cluster_Group_DataBase.hpp
        cluster/cl_MTK_Side_Cluster_ISC_Impl.hpp
        cluster/cl_MTK_Side_Cluster_Input.hpp

        contact/cl_MTK_ContactMapper.hpp
        contact/cl_MTK_ContactMapper_Ray.hpp
        contact/cl_MTK_Contact_Mesh_Editor.hpp
        contact/cl_MTK_Spatial_Indexer.cpp
        contact/cl_MTK_Spatial_Indexer.h
        contact/cl_MTK_Spatial_Indexer_BruteForce.cpp
        contact/cl_MTK_Spatial_Indexer_BruteForce.h

        field/cl_MTK_Field.hpp
        field/cl_MTK_Field_Analytic.hpp
        field/cl_MTK_Field_Discrete.hpp
        field/cl_MTK_Fields_Info.hpp
        field/cl_MTK_Matrix_Field_Info.hpp
        field/cl_MTK_Scalar_Field_Info.hpp

        geometry/cell/cl_MTK_Cell.hpp
        geometry/cell/cl_MTK_Cell_DataBase.hpp
        geometry/cell/cl_MTK_Cell_ISC.hpp
        geometry/cell/cl_MTK_Cell_Info.hpp
        geometry/cell/cl_MTK_Cell_Info_Factory.hpp
        geometry/cell/cl_MTK_Cell_Info_Hex20.hpp
        geometry/cell/cl_MTK_Cell_Info_Hex27.hpp
        geometry/cell/cl_MTK_Cell_Info_Hex64.hpp
        geometry/cell/cl_MTK_Cell_Info_Hex8.hpp
        geometry/cell/cl_MTK_Cell_Info_Quad16.hpp
        geometry/cell/cl_MTK_Cell_Info_Quad8.hpp
        geometry/cell/cl_MTK_Cell_Info_Quad9.hpp
        geometry/cell/cl_MTK_Cell_Info_Tet10.hpp
        geometry/cell/cl_MTK_Cell_Info_Tet4.hpp
        geometry/cell/cl_MTK_Cell_Info_Tri3.hpp
        geometry/cell/cl_MTK_Cell_Info_Tri6.hpp
        geometry/cell/cl_MTK_Cell_Proxy.hpp
        geometry/cl_MTK_Edge.hpp
        geometry/cl_MTK_Facet.hpp
        geometry/cl_MTK_Topology.hpp
        geometry/vertex/cl_MTK_Vertex.hpp
        geometry/vertex/cl_MTK_Vertex_DataBase.hpp
        geometry/vertex/cl_MTK_Vertex_ISC_Impl.hpp
        geometry/vertex/cl_MTK_Vertex_Interpolation.hpp
        geometry/vertex/cl_MTK_Vertex_Interpolation_DataBase.hpp
        geometry/vertex/cl_MTK_Vertex_Proxy.hpp

        integration/cl_MTK_Integration_Coeffs.hpp
        integration/cl_MTK_Integration_Coeffs_Bar_1.hpp
        integration/cl_MTK_Integration_Coeffs_Bar_2.hpp
        integration/cl_MTK_Integration_Coeffs_Bar_3.hpp
        integration/cl_MTK_Integration_Coeffs_Bar_4.hpp
        integration/cl_MTK_Integration_Coeffs_Bar_5.hpp
        integration/cl_MTK_Integration_Coeffs_Bar_6.hpp
        integration/cl_MTK_Integration_Coeffs_Base.hpp
        integration/cl_MTK_Integration_Coeffs_Hex_2x2x2.hpp
        integration/cl_MTK_Integration_Coeffs_Hex_3x3x3.hpp
        integration/cl_MTK_Integration_Coeffs_Hex_4x4x4.hpp
        integration/cl_MTK_Integration_Coeffs_Hex_5x5x5.hpp
        integration/cl_MTK_Integration_Coeffs_Quad_2x2.hpp
        integration/cl_MTK_Integration_Coeffs_Quad_3x3.hpp
        integration/cl_MTK_Integration_Coeffs_Quad_4x4.hpp
        integration/cl_MTK_Integration_Coeffs_Quad_5x5.hpp
        integration/cl_MTK_Integration_Coeffs_Tet_1.hpp
        integration/cl_MTK_Integration_Coeffs_Tet_11.hpp
        integration/cl_MTK_Integration_Coeffs_Tet_15.hpp
        integration/cl_MTK_Integration_Coeffs_Tet_20.hpp
        integration/cl_MTK_Integration_Coeffs_Tet_35.hpp
        integration/cl_MTK_Integration_Coeffs_Tet_4.hpp
        integration/cl_MTK_Integration_Coeffs_Tet_5.hpp
        integration/cl_MTK_Integration_Coeffs_Tet_56.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_1.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_12.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_13.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_16.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_19.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_25.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_3.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_4.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_6.hpp
        integration/cl_MTK_Integration_Coeffs_Tri_7.hpp
        integration/cl_MTK_Integration_Rule.hpp
        integration/cl_MTK_Integrator.hpp
        integration/cl_MTK_Integrator_Test_Polynomial.hpp

        interpolation/cl_Interpolation.hpp
        interpolation/cl_Interpolation.hpp
        interpolation/cl_MTK_Interpolation_Function.hpp
        interpolation/cl_MTK_Interpolation_Function_Base.hpp
        interpolation/cl_MTK_Interpolation_Function_Constant_Bar2.hpp
        interpolation/cl_MTK_Interpolation_Function_Constant_Point.hpp
        interpolation/cl_MTK_Interpolation_Function_Factory.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Bar1.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Bar2.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Bar3.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Bar4.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Hex27.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Hex64.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Hex8.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Quad16.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Quad4.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Quad8.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Quad9.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Tet10.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Tet20.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Tet4.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Tri10.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Tri3.hpp
        interpolation/cl_MTK_Interpolation_Function_Lagrange_Tri6.hpp
        interpolation/cl_MTK_Interpolation_Rule.hpp
        interpolation/cl_MTK_Space_Interpolator.hpp
        interpolation/fn_MTK_Interpolation_Enum_Int_Conversion.hpp

        intersection/cl_MTK_Intersection_Detect.hpp
        intersection/cl_MTK_Intersection_Detect_2D.hpp

        io/cl_MTK_Exodus_IO_Helper.hpp
        io/cl_MTK_Reader_Exodus.hpp
        io/cl_MTK_Writer_Exodus.hpp
        io/cl_MTK_Json_Debug_Output.hpp
        io/json.hpp

        mesh/cl_MTK_Integration_Mesh.hpp
        mesh/cl_MTK_Integration_Mesh_Editor.hpp
        mesh/cl_MTK_Interpolation_Mesh.hpp
        mesh/cl_MTK_Interpolation_Mesh_Editor.hpp
        mesh/cl_MTK_Intersec_Mesh_Data.hpp
        mesh/cl_MTK_Intersection_Mesh.hpp
        mesh/cl_MTK_Mesh.hpp
        mesh/cl_MTK_Mesh_Checker.hpp
        mesh/cl_MTK_Mesh_Core.hpp
        mesh/cl_MTK_Mesh_DataBase_IG.hpp
        mesh/cl_MTK_Mesh_DataBase_IP.hpp
        mesh/cl_MTK_Mesh_Data_Input.hpp
        mesh/cl_MTK_Mesh_Factory.hpp
        mesh/cl_MTK_Mesh_Manager.hpp
        mesh/cl_MTK_Mesh_Pair.hpp
        mesh/cl_MTK_Mesh_Tools.hpp
        mesh/cl_MTK_Surface_Mesh.hpp

        set/cl_MTK_Block_Set.hpp
        set/cl_MTK_Block_Sets_Info.hpp
        set/cl_MTK_Double_Side_Set.hpp
        set/cl_MTK_Node_Sets_Info.hpp
        set/cl_MTK_Set.hpp
        set/cl_MTK_Set_Communicator.hpp
        set/cl_MTK_Sets_Info.hpp
        set/cl_MTK_Side_Set.hpp
        set/cl_MTK_Side_Sets_Info.hpp

        set/mapper/cl_MTK_PointMapper.hpp
        set/mapper/cl_MTK_PointMapper_Ray.hpp

        stk_impl/cl_MTK_Cell_Cluster_STK.hpp
        stk_impl/cl_MTK_Cell_STK.hpp
        stk_impl/cl_MTK_Integration_Mesh_STK.hpp
        stk_impl/cl_MTK_Interpolation_Mesh_STK.hpp
        stk_impl/cl_MTK_Mesh_Core_STK.hpp
        stk_impl/cl_MTK_Mesh_Data_STK.hpp
        stk_impl/cl_MTK_Side_Cluster_STK.hpp
        stk_impl/cl_MTK_Vertex_Interpolation_STK.hpp
        stk_impl/cl_MTK_Vertex_STK.hpp
)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Library

# List library source files
set(LIB_SOURCES
        cl_MTK_Periodic_Boundary_Condition_Helper.cpp

        cluster/cl_MTK_Cell_Cluster_DataBase.cpp
        cluster/cl_MTK_Cell_Cluster_Group.cpp
        cluster/cl_MTK_Cell_Cluster_Group_DataBase.cpp
        cluster/cl_MTK_Cluster_Group.cpp
        cluster/cl_MTK_Double_Side_Cluster.cpp
        cluster/cl_MTK_Nonconformal_Side_Cluster.cpp
        cluster/cl_MTK_Side_Cluster.cpp
        cluster/cl_MTK_Side_Cluster_DataBase.cpp
        cluster/cl_MTK_Side_Cluster_Group.cpp
        cluster/cl_MTK_Side_Cluster_Group_DataBase.cpp
        cluster/cl_MTK_Side_Cluster_ISC_Impl.cpp
        cluster/cl_MTK_Nonconformal_Side_Cluster.cpp

        field/cl_MTK_Field.cpp
        field/cl_MTK_Field_Analytic.cpp
        field/cl_MTK_Field_Discrete.cpp

        geometry/cell/cl_MTK_Cell.cpp
        geometry/cell/cl_MTK_Cell_DataBase.cpp
        geometry/cell/cl_MTK_Cell_ISC.cpp
        geometry/cell/cl_MTK_Cell_Info.cpp
        geometry/cell/cl_MTK_Cell_Info_Factory.cpp
        geometry/cell/cl_MTK_Cell_Info_Hex20.cpp
        geometry/cell/cl_MTK_Cell_Info_Hex27.cpp
        geometry/cell/cl_MTK_Cell_Info_Hex64.cpp
        geometry/cell/cl_MTK_Cell_Info_Hex8.cpp
        geometry/cell/cl_MTK_Cell_Info_Quad16.cpp
        geometry/cell/cl_MTK_Cell_Info_Quad4.cpp
        geometry/cell/cl_MTK_Cell_Info_Quad8.cpp
        geometry/cell/cl_MTK_Cell_Info_Quad9.cpp
        geometry/cell/cl_MTK_Cell_Info_Tet10.cpp
        geometry/cell/cl_MTK_Cell_Info_Tet4.cpp
        geometry/cell/cl_MTK_Cell_Info_Tri3.cpp
        geometry/cell/cl_MTK_Cell_Info_Tri6.cpp
        geometry/vertex/cl_MTK_Vertex_DataBase.cpp
        geometry/vertex/cl_MTK_Vertex_Interpolation_DataBase.cpp

        integration/cl_MTK_Integration_Rule.cpp
        integration/cl_MTK_Integrator.cpp
        integration/cl_MTK_Integrator_Test_Polynomial.cpp

        interpolation/cl_MTK_Interpolation_Function_Factory.cpp
        interpolation/cl_MTK_Interpolation_Rule.cpp
        interpolation/cl_MTK_Space_Interpolator.cpp

        intersection/cl_MTK_Intersection_Detect.cpp
        intersection/cl_MTK_Intersection_Detect_2D.cpp

        io/cl_MTK_Exodus_IO_Helper.cpp
        io/cl_MTK_Reader_Exodus.cpp
        io/cl_MTK_Writer_Exodus.cpp
        io/cl_MTK_Json_Debug_Output.cpp

        mesh/cl_MTK_Integration_Mesh.cpp
        mesh/cl_MTK_Integration_Mesh_Editor.cpp
        mesh/cl_MTK_Integration_Mesh_write_T_matrices.cpp
        mesh/cl_MTK_Interpolation_Mesh_Editor.cpp
        mesh/cl_MTK_Intersec_Mesh_Data.cpp
        mesh/cl_MTK_Intersection_Mesh.cpp
        mesh/cl_MTK_Mesh_Checker.cpp
        mesh/cl_MTK_Mesh_Core.cpp
        mesh/cl_MTK_Mesh_DataBase_IG.cpp
        mesh/cl_MTK_Mesh_DataBase_IP.cpp
        mesh/cl_MTK_Mesh_Factory.cpp
        mesh/cl_MTK_Mesh_Manager.cpp
        mesh/cl_MTK_Mesh_Pair.cpp
        mesh/cl_MTK_Surface_Mesh.cpp

        set/cl_MTK_Set_Communicator.cpp

        set/mapper/cl_MTK_PointMapper.cpp
        set/mapper/cl_MTK_PointMapper_Ray.cpp

        stk_impl/cl_MTK_Cell_Cluster_STK.cpp
        stk_impl/cl_MTK_Integration_Mesh_STK.cpp
        stk_impl/cl_MTK_Interpolation_Mesh_STK.cpp
        stk_impl/cl_MTK_Mesh_Core_STK.cpp
        stk_impl/cl_MTK_Side_Cluster_STK.cpp
        stk_impl/cl_MTK_Vertex_Interpolation_STK.cpp
)

# List library dependencies
set(LIB_DEPENDENCIES
        ${COM}-lib
        ${LINALG}-lib
        ${HMR}-lib
        ${TOL}-lib
        ${ENM}-lib
)

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MTK}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${MTK}-lib PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${MTK_INTERNAL_INCLUDES};>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mesh>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/set>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/geometry>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/geometry/cell>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/geometry/vertex>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/field>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cluster>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/io>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/intersection>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interpolation>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/integration>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/stk_impl>"
        "$<INSTALL_INTERFACE:${${MTK}_HEADER_INSTALL_DIR}>")
target_link_libraries(${MTK}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MTK}-lib PROPERTIES OUTPUT_NAME ${MTK})

foreach (TPL ${MTK_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${MTK}-lib
            PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MTK}-lib
            PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${MTK}-lib
            PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach ()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach (TPL ${MTK_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    string(APPEND MTK_CONFIG_TPL_DEPENDENCIES
            "include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n")
endforeach ()

# Makes sure all moris library dependencies are included in the config file
foreach (LIB_TARGET ${LIB_DEPENDENCIES})
    string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND MTK_CONFIG_MORIS_DEPENDENCIES
            "find_package(${LIB} QUIET HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n")
endforeach ()

foreach (LIB_TARGET ${LIB_DEPENDENCIES})
    string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND MTK_CONFIG_MORIS_DEPENDENCIES
            "find_package(${LIB} REQUIRED HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n")
endforeach ()

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
        COMPATIBILITY AnyNewerVersion)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# mtk Installation

# Install *Config and *ConfigVersion files
install(FILES ${${MTK}_CONFIG_BUILD_DIR}/${MTK}Config.cmake
        ${${MTK}_CONFIG_BUILD_DIR}/${MTK}ConfigVersion.cmake
        DESTINATION ${${MTK}_CONFIG_INSTALL_DIR})

# Install header files
install(FILES ${HEADERS}
        DESTINATION ${${MTK}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${MTK}-lib
        DESTINATION ${${MTK}_LIB_INSTALL_DIR}
        EXPORT ${MTK}Targets)

# Export library into *Targets file
install(EXPORT ${MTK}Targets
        DESTINATION ${${MTK}_CONFIG_INSTALL_DIR})

