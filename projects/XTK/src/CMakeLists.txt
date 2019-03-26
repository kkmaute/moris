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
    xtk/cl_MTK_Mesh_XTK_Impl.hpp	
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
    ${ALG}
    ${MTK}
    ${MTK}/stk_impl
    ${MTK}/${MAP}
    ${COM}
    ${HMR}
    )


# Include depedency directories
foreach(XTK_INCLUDE ${XTK_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${XTK_INCLUDE}/src)
endforeach()
     
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                                                                                                                                                                                    
# XTK Library                                                                                                                                                                                                                                                                  

# List source files                                                                                                                                                                                                                                                            
set(LIB_SOURCES
       xtk/cl_XTK_Model.cpp
       xtk/cl_XTK_Background_Mesh.cpp
       xtk/cl_XTK_Child_Mesh.cpp
       xtk/cl_XTK_Enrichment.cpp
       xtk/cl_MTK_Cell_XTK_Impl.cpp
       xtk/cl_MTK_Mesh_XTK_Impl.cpp
       xtk/cl_MTK_Vertex_XTK_Impl.cpp
       xtk/cl_MTK_Vertex_Interpolation_XTK_Impl.cpp
       xtk/cl_XTK_Paramfile.cpp)

# List library dependencies
set(LIB_DEPENDENCIES 
	${COM}-lib
	${MTK}-lib
	${ALG}-lib
	${HMR}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS})

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${XTK}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${XTK}-lib ${LIB_DEPENDENCIES})
set_target_properties(${XTK}-lib PROPERTIES OUTPUT_NAME ${XTK})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#############################################################################
# Create executable for plato CDFEM
#############################################################################
# XTK Executable
# List source files
set(PLATO_EXE_SOURCES
        main_plato.cpp)

# List executable dependencies
set(PLATO_EXE_DEPENDENCIES
    ${XTK}-lib
    ${MTK}-lib
    ${COM}-lib
    ${ALG}-lib
	${MORIS_BASE_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}	
    )    
 
# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${bin})

add_executable(XTK_PLATO ${PLATO_EXE_SOURCES})
target_link_libraries(XTK_PLATO ${PLATO_EXE_DEPENDENCIES})
set_target_properties(XTK_PLATO PROPERTIES OUTPUT_NAME XTK_Plato.exe )


# Hole seeding executable
# XTK Executable
# List source files
set(HOLE_SEED_EXE_SOURCES
        main_hole_seed.cpp)

# List executable dependencies
set(HOLE_EXE_DEPENDENCIES
    ${XTK}-lib
    ${MTK}-lib
    ${COM}-lib
    ${ALG}-lib
	${MORIS_BASE_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}	
    )    
 
# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${bin})

#############################################################################
# Create executable for XFEM
#############################################################################
set(XTK_EXE_SOURCES
        main.cpp)

# List executable dependencies
set(XTK_EXE_DEPENDENCIES
    ${XTK}-lib
    ${MTK}-lib
    ${COM}-lib
    ${ALG}-lib
	${MORIS_BASE_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}	
    )    

add_executable(XTK_EXE ${XTK_EXE_SOURCES})
target_link_libraries(XTK_EXE ${XTK_EXE_DEPENDENCIES})
set_target_properties(XTK_EXE PROPERTIES OUTPUT_NAME XTK.exe )



#############################################################################
# Create executable for hole seeding
#############################################################################
set(HOLE_SEED_EXE_SOURCES
        main_hole_seed.cpp)

# List executable dependencies
set(HOLE_EXE_DEPENDENCIES
    ${XTK}-lib
    ${MTK}-lib
    ${COM}-lib
    ${ALG}-lib
	${MORIS_BASE_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}	
    )    
 
# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${bin})
add_executable(XTK_HOLE_SEED ${HOLE_SEED_EXE_SOURCES})
target_link_libraries(XTK_HOLE_SEED ${HOLE_EXE_DEPENDENCIES})
set_target_properties(XTK_HOLE_SEED PROPERTIES OUTPUT_NAME XTK_Hole_Seed.exe )



if(BUILD_HMR)
set(HMR_EXE_SOURCES
        main_hmr.cpp)

# List executable dependencies
set(HMR_EXE_DEPENDENCIES
	${XTK}-lib
    ${HMR}-lib
    ${MAP}-lib
    ${MTK}-lib
    ${MDL}-lib # for L2
    ${INT}-lib # for L2
    ${MSI}-lib # for L2
    ${NLA}-lib # for L2
#    ${ALG}-lib #MTK
#    ${TOL}-lib #MTK
    ${SDF}-lib # for SDF
#    ${XTK}-lib #MTK
#    ${STK}-lib #MTKz
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}  
    )    
 
# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${bin})

# Create executable
add_executable(XTK_HMR ${HMR_EXE_SOURCES})
target_link_libraries(XTK_HMR ${HMR_EXE_DEPENDENCIES})
set_target_properties(XTK_HMR PROPERTIES OUTPUT_NAME XTK_HMR.exe )

endif()


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
install(TARGETS XTK_PLATO
    DESTINATION ${EXE_INSTALL_DIR}/${XTK}
    EXPORT ${XTK}Targets )

# Install executable
install(TARGETS XTK_HOLE_SEED
    DESTINATION ${EXE_INSTALL_DIR}/${XTK}
    EXPORT ${XTK}Targets )
    
    
# Install executable
install(TARGETS XTK_EXE
    DESTINATION ${EXE_INSTALL_DIR}/${XTK}
    EXPORT ${XTK}Targets )
    
    
if(BUILD_HMR)
install(TARGETS XTK_HMR
    DESTINATION ${EXE_INSTALL_DIR}/${XTK}
    EXPORT ${XTK}Targets )
endif()    

# Export executable in to *Targets file
install(EXPORT ${XTK}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${XTK} )
