# XTK Source --------------------------------------------------------------                                                                                                                                                                                                    
# -------------------------------------------------------------------------                                                                                                                                                                                                    

# Set version                                                                                                                                                                                                                                                                  
set(${XTK}_VERSION ${MORIS_VERSION})                                                                                                                                                                                                                                           

# List header files                                                                                                                                                                                                                                                            
set(HEADERS
    algorithms/fn_num_digits.hpp
    core/xtk_typedefs.hpp
    geomeng/cl_MGE_Enums.hpp
    geomeng/cl_MGE_Geometry_Engine.hpp
    geomeng/cl_MGE_Geometry_Object_Manager.hpp
    geomeng/cl_MGE_Geometry_Object.hpp
    geomeng/fn_Triangle_Geometry.hpp
    geometry/cl_Geometry.hpp
    geometry/cl_Analytic_Level_Set_Plane.hpp
    geometry/cl_Composite_Fiber_Straight_1.hpp
    geometry/cl_Composite_Fiber_Straight_2.hpp
    geometry/cl_Composite_Fiber_Straight_3.hpp
    geometry/cl_Composite_Fiber_Wave_1.hpp
    geometry/cl_Composite_Fiber_Wave_2.hpp
    geometry/cl_Composite_Fiber_Wave_3.hpp
    geometry/cl_Composite_Fiber.hpp
    geometry/cl_Discrete_Level_Set.hpp
    geometry/cl_Gyroid.hpp
    geometry/cl_Multi_Cylinder.hpp
    geometry/cl_Plane.hpp
    geometry/cl_Sphere.hpp
    geometry/cl_Spiral.hpp
    ios/cl_Logger.hpp
    linalg/cl_XTK_Matrix_Base_Utilities.hpp
    mesh/cl_Mesh_Bucket.hpp
    mesh/cl_Mesh_Data_MTK.hpp
    mesh/cl_Mesh_Entity.hpp
    mesh/cl_Mesh_External_Entity_Data_Stk.hpp
    mesh/cl_Mesh_Node_Set_Input.hpp
    mesh/cl_Mesh_Parallel_Data_Stk.hpp
    mesh/cl_Mesh_Side_Set_Input.hpp
    tools/fn_verify_tet_topology.hpp
    tools/cl_Interpolaton.hpp
    tools/cl_MPI_Tools.hpp
    tools/fn_approximate.hpp
    tools/fn_bubble_sort.hpp
    tools/fn_hex_8_volume.hpp
    tools/fn_matrix_to_offset.hpp
    tools/fn_Pairing.hpp
    tools/fn_tet_volume.hpp
    topology/cl_XTK_Basis_Function.hpp
    topology/cl_XTK_Topology.hpp
    topology/cl_XTK_Edge_Topology.hpp
    topology/cl_XTK_Hexahedron_8_Basis_Function.hpp
    topology/cl_XTK_Hexahedron_8_Topology.hpp
    topology/cl_XTK_Linear_Basis_Functions.hpp
    topology/cl_XTK_Quad_4_Basis_Function.hpp
    topology/cl_XTK_Quad_4_Topology.hpp
    xtk/cl_XTK_Active_Process_Manager.hpp
    xtk/cl_XTK_Background_Mesh.hpp
    xtk/cl_XTK_Child_Mesh_Modification_Template.hpp
    xtk/cl_XTK_Child_Mesh.hpp
    xtk/cl_XTK_Cut_Mesh.hpp
    xtk/cl_XTK_Downward_Inheritance.hpp
    xtk/cl_XTK_Enrichment.hpp
    xtk/cl_XTK_Entity_Tracker.hpp
    xtk/cl_XTK_Entity.hpp
    xtk/cl_XTK_Enums.hpp
    xtk/cl_XTK_External_Mesh_Data.hpp
    xtk/cl_XTK_Face_Registry.hpp
    xtk/cl_XTK_Input_Generator.hpp
    xtk/cl_XTK_Interface_Element.hpp
    xtk/cl_XTK_Model_Export.hpp
    xtk/cl_XTK_Model.hpp
    xtk/cl_XTK_Node.hpp
    xtk/cl_XTK_Output_Options.hpp
    xtk/cl_XTK_Parallel_Mesh_Data.hpp
    xtk/cl_XTK_Pending_Node.hpp
    xtk/cl_XTK_Phase_Table.hpp
    xtk/cl_XTK_Request_Handler.hpp
    xtk/cl_XTK_Sensitivity.hpp
    xtk/cl_XTK_Tetra4_Connectivity.hpp
    xtk/fn_assemble_boundary_subphase_constraint.hpp
    xtk/fn_compute_xtk_model_volumes.hpp
    xtk/fn_create_edges_from_element_to_node.hpp
    xtk/fn_create_faces_from_element_to_node.hpp
    xtk/fn_generate_cross_boundary_pairs.hpp
    xtk/fn_generate_element_to_element.hpp
    xtk/fn_generate_shared_face_element_graph.hpp
    xtk/fn_local_child_mesh_flood_fill.hpp
    xtk/fn_mesh_flood_fill.hpp
    xtk/fn_prune_element_to_element.hpp
    xtk/fn_write_element_ownership_as_field.hpp)                                                                                                                                                                                                                                                                   

# List header dependencies                                                                                                                                                                                                                                                     
set(XTK_INCLUDES                                                                                                                                                                                                                                                               
    "assert"                                                                                                                                                                                                                                                                   
    "containers"                                                                                                                                                                                                                                                               
    "core"                                                                                                                                                                                                                                                                     
    "geomeng"                                                                                                                                                                                                                                                                  
    "geometry"                                                                                                                                                                                                                                                                 
    "ios"                                                                                                                                                                                                                                                                      
    "linalg"                                                                                                                                                                                                                                                                   
    "mesh"                                                                                                                                                                                                                                                                     
    "tools"                                                                                                                                                                                                                                                                    
    "topology"                                                                                                                                                                                                                                                                 
    "xtk"                                                                                                                                                                                                                                                                      
    )                                                                                                                                                                                                                                                                          

# Include depedency directories                                                                                                                                                                                                                                                
foreach(XTK_INCLUDE ${XTK_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/XTK/src/${XTK_INCLUDE})                                                                                                                                                                                                           
endforeach()                                                                                                                                                                                                                                                                   

include_directories(${MORIS_PACKAGE_DIR}/XTK/include)                                                                                                                                                                                                                          



# List header dependencies
set(XTK_INCLUDES
    ${LINALG}
    ${MTK}
    ${MTK}/stk_impl
    ${COM}
    )


# Include depedency directories
foreach(XTK_INCLUDE ${XTK_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${XTK_INCLUDE}/src)
endforeach()
     
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                                                                                                    
# XTK Library                                                                                                                                                                                                                                                                  

# List source files                                                                                                                                                                                                                                                            
# set(LIB_SOURCES
#       linalg/cl_XTK_Linalg_Enums.${XTK_SOURCE_EXT)

# List library dependencies
# set(LIB_DEPENDENCIES "")

# Set the output path for library
# set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
#add_library(${XTK}-lib STATIC ${LIB_SOURCES})
#target_link_libraries(${XTK}-lib ${LIB_DEPENDENCIES})
#set_target_properties(${XTK}-lib PROPERTIES OUTPUT_NAME ${XTK})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Executable
# List source files
set(EXE_SOURCES
        main.cpp)

# List executable dependencies
set(EXE_DEPENDENCIES
    ${MTK}-lib
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    )    
 
# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${bin})

# Create executable
add_executable(${XTK} ${EXE_SOURCES})
target_link_libraries(${XTK} ${EXE_DEPENDENCIES})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# XTK Config Files

# Build *Config and *ConfigVersion files
 include(CMakePackageConfigHelpers)

 
 configure_package_config_file(
     ${XTK}Config.cmake.in
     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}Config.cmake
     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${XTK}
     PATH_VARS LIB_INSTALL_DIR )
 
 write_basic_package_version_file(
     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}ConfigVersion.cmake
     COMPATIBILITY AnyNewerVersion )
 
 # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 # XTK Installation
 
 # Install *Config and *ConfigVersion files
 install(FILES ${${XTK}_CONFIG_BUILD_DIR}/${XTK}Config.cmake 
     ${${XTK}_CONFIG_BUILD_DIR}/${XTK}ConfigVersion.cmake
     DESTINATION ${CONFIG_INSTALL_DIR}/${XTK} )
# Install header files
 install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${XTK} )

# Install executable
install(TARGETS ${XTK}
    DESTINATION ${EXE_INSTALL_DIR}/${XTK}
    EXPORT ${XTK}Targets )

# Export executable in to *Targets file
install(EXPORT ${XTK}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${XTK} )
