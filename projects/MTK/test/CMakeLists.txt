# MTK Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
	cl_MTK_Mesh.cpp
	cl_MTK_Blocks.cpp
	cl_MTK_Side_Set.cpp
	cl_MTK_Exodus_IO_Helper.cpp
	cl_MTK_Mesh_From_Data.cpp
	cl_MTK_Vertex.cpp
	cl_MTK_Cell.cpp
	cl_MTK_Side_Cluster.cpp
	cl_MTK_Mesh_New.cpp
	cl_MTK_Interpolation_Mesh.cpp
	cl_MTK_Cell_Cluster.cpp
	cl_MTK_Double_Side_Cluster.cpp
	cl_MTK_Writer_Exodus_Test.cpp
	UT_MTK_Mesh_Comm_Maps.cpp)

# List additional includes
include_directories(${MORIS_PACKAGE_DIR}/MTK/src/stk_impl)

# List test dependencies
set(TEST_DEPENDENCIES
    test-libs
    ${ALG}-lib
    ${MORIS_BASE_LIBS}
    ${MTK}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MTK_TEST ${MTK}-test)

# Create executable
add_executable(${MTK_TEST} ${TEST_SOURCES})
target_link_libraries(${MTK_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${MTK_TEST} PROPERTIES OUTPUT_NAME ${MTK_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(MTK_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${MTK_TEST_PROCS})
		set(MTK_PARALLEL_TEST ${MTK_TEST}-${PROCS}-procs)
		add_test(NAME ${MTK_PARALLEL_TEST} 
			COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} $<TARGET_FILE:${MTK_TEST}>)
	endforeach()
endif()