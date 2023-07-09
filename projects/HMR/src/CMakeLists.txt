#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# HMR Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${HMR}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_HMR_Arguments.hpp
    cl_HMR_Background_Edge.hpp
    cl_HMR_Background_Element_Base.hpp
    cl_HMR_Background_Element.hpp
    cl_HMR_Background_Facet.hpp
    cl_HMR_Background_Mesh_1D.hpp
    cl_HMR_Background_Mesh_2D.hpp
    cl_HMR_Background_Mesh_3D.hpp
    cl_HMR_Background_Mesh_Base.hpp
    cl_HMR_Background_Mesh.hpp
    cl_HMR_Basis.hpp
    cl_HMR_BSpline_Element_Hex27.hpp
    cl_HMR_BSpline_Element_Hex64.hpp
    cl_HMR_BSpline_Element_Hex8.hpp
    cl_HMR_BSpline_Element_Quad16.hpp
    cl_HMR_BSpline_Element_Quad4.hpp
    cl_HMR_BSpline_Element_Quad9.hpp
    cl_HMR_BSpline_Element.hpp
    cl_HMR_BSpline_Mesh_Base.hpp
    cl_HMR_BSpline_Mesh.hpp
    cl_HMR_BSpline.hpp
    cl_HMR_Cell_Cluster.hpp
    cl_HMR_Database.hpp
    cl_HMR_Domain.hpp
    cl_HMR_Edge.hpp
    cl_HMR_Element.hpp
    cl_HMR_Facet.hpp
    cl_HMR_Factory.hpp
    cl_HMR_Field_Param.hpp
    cl_HMR_Field.hpp
    cl_HMR_File.hpp
    cl_HMR_Knot.hpp
    cl_HMR_Lagrange_Edge.hpp
    cl_HMR_Lagrange_Edge2.hpp
    cl_HMR_Lagrange_Edge3.hpp
    cl_HMR_Lagrange_Edge4.hpp
    cl_HMR_Lagrange_Element_Hex27.hpp
    cl_HMR_Lagrange_Element_Hex64.hpp
    cl_HMR_Lagrange_Element_Hex8.hpp
    cl_HMR_Lagrange_Element_Quad16.hpp
    cl_HMR_Lagrange_Element_Quad4.hpp
    cl_HMR_Lagrange_Element_Quad9.hpp
    cl_HMR_Lagrange_Element.hpp
    cl_HMR_Lagrange_Facet_Line2.hpp
    cl_HMR_Lagrange_Facet_Line3.hpp
    cl_HMR_Lagrange_Facet_Line4.hpp
    cl_HMR_Lagrange_Facet_Quad16.hpp
    cl_HMR_Lagrange_Facet_Quad4.hpp
    cl_HMR_Lagrange_Facet_Quad9.hpp
    cl_HMR_Lagrange_Facet.hpp
    cl_HMR_Lagrange_Mesh_Base.hpp
    cl_HMR_Lagrange_Mesh.hpp
    cl_HMR_Lagrange_Node_Interpolation.hpp
    cl_HMR_Lagrange_Node.hpp
    cl_HMR_Mesh_Base.hpp
    cl_HMR_Mesh_Integration.hpp
    cl_HMR_Mesh_Interpolation.hpp
    cl_HMR_Mesh.hpp
    cl_HMR_Parameters.hpp
    cl_HMR_Paramfile.hpp
    cl_HMR_Side_Cluster.hpp
    cl_HMR_Side_Set.hpp
    cl_HMR_State.hpp
    cl_HMR_STK.hpp
    cl_HMR_T_Matrix.hpp
    cl_HMR_T_Matrix_2.hpp
    cl_HMR.hpp
    fn_HMR_Background_Element_Edges_3D.hpp
    fn_HMR_Background_Element_Neighbors_2D.hpp
    fn_HMR_Background_Element_Neighbors_3D.hpp
    fn_HMR_Exec_dump_fields.hpp
    fn_HMR_Exec_dump_meshes.hpp
    fn_HMR_Exec_initialize_fields.hpp
    fn_HMR_Exec_load_parameters.hpp
    fn_HMR_Exec_load_user_library.hpp
    fn_HMR_Exec_perform_mapping.hpp
    fn_HMR_get_basis_neighbors_2d.hpp
    fn_HMR_get_basis_neighbors_3d.hpp
    fn_HMR_refinement_transition_locations.hpp
    HMR_Globals.hpp
    HMR_Tools.hpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Library

# List library source files
set(LIB_SOURCES
    cl_HMR_Arguments.cpp
    cl_HMR_Background_Edge.cpp
    cl_HMR_Background_Facet.cpp
    cl_HMR_Background_Mesh_Base.cpp
    cl_HMR_BSpline_Mesh_Base.cpp
    cl_HMR_Database.cpp
    cl_HMR.cpp
    cl_HMR_Edge.cpp
    cl_HMR_Element.cpp
    cl_HMR_Facet.cpp
    cl_HMR_Factory.cpp
    cl_HMR_Field.cpp
    cl_HMR_File.cpp
    cl_HMR_Lagrange_Mesh_Base.cpp
    cl_HMR_Lagrange_Node_Interpolation.cpp
    cl_HMR_Mesh_Base.cpp
    cl_HMR_Mesh.cpp
    cl_HMR_Parameters.cpp
    cl_HMR_Paramfile.cpp
    cl_HMR_STK.cpp
    cl_HMR_T_Matrix.cpp
    cl_HMR_T_Matrix_2.cpp
    cl_HMR_Cell_Cluster.cpp
    cl_HMR_Side_Cluster.cpp
   )

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${CNT}-lib
    ${IOS}-lib
    ${ALG}-lib
    ${MTK}-lib
    ${MAP}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${HMR}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${HMR}-lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${${HMR}_HEADER_INSTALL_DIR}> )
target_link_libraries(${HMR}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${HMR}-lib PROPERTIES OUTPUT_NAME ${HMR})

foreach(TPL ${HMR_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${HMR}-lib
        PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${HMR}-lib
        PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${HMR}-lib
        PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# project Executable

# List executable source files
set(EXE_SOURCES
    main.cpp )

# List executable dependencies
set(EXE_DEPENDENCIES
    ${LINALG}-lib
    ${CNT}-lib
    ${IOS}-lib
    ${HMR}-lib
    ${MTK}-lib
    ${DLA}-lib
    )
    
# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
if(BUILD_HMR)
    add_executable(${HMR}-exe ${EXE_SOURCES})
    target_link_libraries(${HMR}-exe PRIVATE ${EXE_DEPENDENCIES})
    set_target_properties(${HMR}-exe PROPERTIES OUTPUT_NAME hmr )
endif()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${HMR_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    string(APPEND HMR_CONFIG_TPL_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
    # skip if -lfoo format
    string(FIND ${LIB_TARGET} "-" pos)
    if(pos EQUAL "0")
        continue()
    endif()
    
    string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
    string(APPEND HMR_CONFIG_MORIS_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(HMR_COMPILE_DEFINITIONS ${HMR}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${HMR}Config.cmake.in
    ${${HMR}_CONFIG_BUILD_DIR}/${HMR}Config.cmake
    INSTALL_DESTINATION 
    ${${HMR}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${HMR}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${HMR}_CONFIG_BUILD_DIR}/${HMR}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# hmr Installation

# Install *Config and *ConfigVersion files
install(FILES ${${HMR}_CONFIG_BUILD_DIR}/${HMR}Config.cmake 
    ${${HMR}_CONFIG_BUILD_DIR}/${HMR}ConfigVersion.cmake
    DESTINATION ${${HMR}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${HMR}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${HMR}-lib
    DESTINATION ${${HMR}_LIB_INSTALL_DIR}
    EXPORT ${HMR}Targets )

# Export library into *Targets file
install(EXPORT ${HMR}Targets
    DESTINATION ${${HMR}_CONFIG_INSTALL_DIR} )

