# XTK Tests -------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
	main.cpp
	assert/cl_Assert.cpp
	geomeng/cl_Geometry_Engine.cpp
	geomeng/cl_Geometry.cpp
	ios/cl_Logger.cpp
	mesh/cl_Mesh_Data.cpp
	mesh/cl_Mesh_Side_Set_Input.cpp
	tools/cl_MPI_Tools.cpp
	tools/fn_bubble_sort.cpp
	xtk/XTK_Tutorials.cpp
	xtk/cl_XTK_Child_Mesh_NH_Permutations.cpp
	xtk/cl_XTK_Child_Mesh.cpp
	xtk/cl_XTK_Cut_Mesh_Modification.cpp
	xtk/cl_XTK_Cut_Mesh_RegSub.cpp
	xtk/cl_XTK_Cut_Mesh.cpp
	xtk/cl_XTK_Downward_Inheritance.cpp
	xtk/cl_XTK_Element_Tie_Constraints.cpp
	xtk/cl_XTK_Enrichment.cpp
	xtk/cl_XTK_FEM.cpp
	xtk/cl_XTK_Interface_Sides.cpp
	xtk/cl_XTK_Mesh.cpp
	xtk/cl_XTK_Model_Output_Options.cpp
	xtk/cl_XTK_Model_Recursion.cpp
	xtk/cl_XTK_Model_TET10.cpp
	xtk/cl_XTK_Model.cpp
	xtk/cl_XTK_Sensitivity.cpp
	xtk/cl_XTK_Tet_Decomposition.cpp
	xtk/fn_create_edges_from_element_to_node.cpp
	xtk/fn_create_faces_from_element_to_node.cpp
	xtk/fn_flood_fill.cpp
	xtk/fn_generate_element_to_element.cpp)

# List additional includes
include_directories(../src)
include_directories(../include)


# List test dependencies
# List test includes
set(TEST_INCLUDES
    ${LINALG} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()


set(TEST_DEPENDENCIES
	${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${XTK}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${XTK}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${XTK}-test
    COMMAND ${XTK}-test${EXE_EXT} )
