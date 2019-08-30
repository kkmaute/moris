# XTK Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files

set(TEST_SOURCES
	main.cpp
	xtk/UT_MTK_Mesh_XTK_Impl.cpp
	xtk/UT_XTK_Child_Mesh_NH_Permutations.cpp
	xtk/UT_XTK_Child_Mesh.cpp
	xtk/UT_XTK_Cut_Mesh_Modification.cpp
	xtk/UT_XTK_Cut_Mesh_RegSub.cpp
	xtk/UT_XTK_Cut_Mesh.cpp
	xtk/UT_XTK_Downward_Inheritance.cpp
	xtk/UT_XTK_Enrichment.cpp
	xtk/UT_XTK_Ghost_Stabilization.cpp
	#xtk/cl_XTK_Multigrid.cpp
	xtk/UT_XTK_Model.cpp
	#xtk/UT_XTK_Model_Coincidence.cpp
	#xtk/UT_XTK_Model_Recursion.cpp
	#xtk/cl_XTK_Tetra4_Connectivity.cpp
	xtk/UT_FRF_Matching.cpp
	xtk/UT_MTK_Cell_Cluster_XTK_Impl.cpp
	#xtk/UT_XTK_Model_Unzip.cpp
	xtk/UT_XTK_Model_Export.cpp
	#xtk/UT_XTK_Model_Extract_Surf_Mesh.cpp
	xtk/UT_XTK_Model_Mesh_Clusters.cpp
	xtk/UT_AD_Scratch.cpp
	)
	
if(BUILD_HMR)
set(XTK_HMR_TEST_SOURCES
	xtk/UT_XTK_HMR.cpp)
	
	set(HMR_TEST_DEPENDENCIES
	${HMR}-lib
    )
endif()

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
set(TEST_DEPENDENCIES
	test-libs
	${ALG}-lib
	${MTK}-lib
	${XTK}-lib
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(XTK_TEST ${XTK}-test)

# Create executable
add_executable(${XTK_TEST} ${TEST_SOURCES} ${XTK_HMR_TEST_SOURCES})
target_link_libraries(${XTK_TEST} PRIVATE ${TEST_DEPENDENCIES} ${HMR_TEST_DEPENDENCIES})
set_target_properties(${XTK_TEST} PROPERTIES OUTPUT_NAME ${XTK_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(XTK_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${XTK_TEST_PROCS})
		set(XTK_PARALLEL_TEST ${XTK_TEST}-${PROCS}-procs)
		add_test(NAME ${XTK_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} $<TARGET_FILE:${XTK_TEST}>)
	endforeach()
endif()
