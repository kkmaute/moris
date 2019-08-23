# XTK Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
#set(TEST_SOURCES
#	main.cpp
#	assert/cl_Assert.cpp
#	geomeng/cl_Geometry_Engine.cpp
#	geomeng/cl_Geometry.cpp
#	ios/cl_Logger.cpp
#	tools/cl_MPI_Tools.cpp
#	tools/fn_bubble_sort.cpp
#	xtk/XTK_Tutorials.cpp
#	xtk/cl_XTK_Child_Mesh_NH_Permutations.cpp
#	xtk/cl_XTK_Child_Mesh.cpp
#	xtk/cl_XTK_Cut_Mesh_Modification.cpp
#	xtk/cl_XTK_Cut_Mesh_RegSub.cpp
#	xtk/cl_XTK_Cut_Mesh.cpp
#	xtk/cl_XTK_Downward_Inheritance.cpp
#	xtk/cl_XTK_Element_Tie_Constraints.cpp
#	xtk/cl_XTK_Enrichment.cpp
#	xtk/cl_XTK_Interface_Sides.cpp
#	xtk/cl_XTK_Mesh.cpp
#	xtk/cl_XTK_Model_Output_Options.cpp
#	xtk/cl_XTK_Model_Recursion.cpp
#	xtk/cl_XTK_Model_TET10.cpp
#	xtk/cl_XTK_Model.cpp
#	xtk/cl_XTK_Sensitivity.cpp
#	xtk/cl_XTK_Tet_Decomposition.cpp
#	xtk/fn_create_edges_from_element_to_node.cpp
#	xtk/fn_create_faces_from_element_to_node.cpp
#	xtk/fn_flood_fill.cpp
#	xtk/fn_generate_element_to_element.cpp)

set(TEST_SOURCES
	main.cpp
	xtk/UT_MTK_Mesh_XTK_Impl.cpp
	xtk/UT_XTK_Child_Mesh_NH_Permutations.cpp
	xtk/UT_XTK_Child_Mesh.cpp
	xtk/UT_XTK_Cut_Mesh_Modification.cpp
	xtk/UT_XTK_Cut_Mesh_RegSub.cpp
	xtk/UT_XTK_Cut_Mesh.cpp
	xtk/UT_XTK_Downward_Inheritance.cpp
	xtk/UT_XTK_Element_Tie_Constraints.cpp
	xtk/UT_XTK_Enrichment.cpp
	xtk/UT_XTK_Ghost_Stabilization.cpp
	#xtk/cl_XTK_Multigrid.cpp
	xtk/UT_XTK_Model.cpp
	xtk/UT_XTK_Model_Coincidence.cpp
	xtk/UT_XTK_Model_Recursion.cpp
	xtk/UT_XTK_Sensitivity.cpp
	#xtk/cl_XTK_Tetra4_Connectivity.cpp
	xtk/UT_FRF_Matching.cpp
	xtk/UT_MTK_Cell_Cluster_XTK_Impl.cpp
	xtk/UT_XTK_Model_Unzip.cpp
	xtk/UT_XTK_Model_Export.cpp
	xtk/UT_XTK_Model_Extract_Surf_Mesh.cpp
	xtk/UT_XTK_Model_Mesh_Clusters.cpp
	xtk/UT_AD_Scratch.cpp
	)

# List additional includes
include_directories(../src)
include_directories(../include)


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


# List test dependencies
# List test includes
set(TEST_INCLUDES
    #${LINALG}
    #${ALG}
    #${MTK}
    #${MTK}/stk_impl
    #${MTK}/${MAP}
    #${COM}
    #${HMR}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()


set(TEST_DEPENDENCIES
	test-libs
	${ALG}-lib
	${MTK}-lib
	${XTK}-lib
    #${HMR}-lib
    #${MAP}-lib
    #${MTK}-lib
    #${MDL}-lib # for L2
    #${INT}-lib # for L2
    #${MSI}-lib # for L2
    #${NLA}-lib # for L2
#    ${ALG}-lib #MTK
#    ${TOL}-lib #MTK
    #${SDF}-lib # for SDF
#    ${XTK}-lib #MTK
#    ${STK}-lib #MTKz
    #${COM}-lib
    ${MORIS_BASE_LIBS}
    #${MORIS_ACML_LAPACK_MKL_LIBS}
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_TRILINOS_LIBS}
    #${MORIS_ARPACK_LIBS}
    #${MORIS_BOOST_LIBS} 
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE TEST_DEPENDENCIES)
list(REMOVE_DUPLICATES TEST_DEPENDENCIES)
list(REVERSE TEST_DEPENDENCIES)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(XTK_TEST ${XTK}-test)

# Create executable
add_executable(${XTK_TEST} ${TEST_SOURCES})
target_link_libraries(${XTK_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${XTK_TEST} PROPERTIES OUTPUT_NAME ${XTK_TEST}${EXE_EXT})

# Mark as test
add_test(NAME ${XTK_TEST}
    COMMAND $<TARGET_FILE:${XTK_TEST}> )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(XTK_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${XTK_TEST_PROCS})
		set(XTK_PARALLEL_TEST ${XTK_TEST}-${PROCS}-procs)
		add_test(NAME ${XTK_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${XTK_TEST}>)
	endforeach()
endif()
