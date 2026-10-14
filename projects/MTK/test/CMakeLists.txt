# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# MTK Tests

include_directories(MTK_Test_Proxy)

# Sources.
set(TEST_SOURCES
    test_main.cpp
    cl_MTK_Mesh.cpp
    cl_MTK_Blocks.cpp
    cl_MTK_Side_Set.cpp
    #cl_MTK_Exodus_IO_Helper.cpp
    cl_MTK_Mesh_From_Data.cpp
    cl_MTK_Vertex.cpp
    cl_MTK_Cell.cpp
    cl_MTK_Side_Cluster.cpp
    cl_MTK_Mesh_New.cpp
    cl_MTK_Interpolation_Mesh.cpp
    cl_MTK_Cell_Cluster.cpp
    cl_MTK_Double_Side_Cluster.cpp
    cl_MTK_Writer_Exodus_Test.cpp
    UT_MTK_Cell_Info_Hex27.cpp
    UT_MTK_Cell_Info_Hex64.cpp
    UT_MTK_Cell_Info_Hex8.cpp
    UT_MTK_Cell_Info_Quad4.cpp
    UT_MTK_Cell_Info_Quad9.cpp
    UT_MTK_Cell_Info_Quad16.cpp
    UT_MTK_Cell_Info_Tet4.cpp
    UT_MTK_Cell_Info_Tet10.cpp
    UT_MTK_Cell_Info_Tri3.cpp
    UT_MTK_Cell_Info_Tri6.cpp
    UT_MTK_Connectivity.cpp
    UT_MTK_Mesh_Comm_Maps.cpp
    UT_MTK_Field.cpp
    UT_MTK_Cell_Area.cpp
    UT_MTK_Cell_Shape.cpp
    UT_MTK_Cell_Shape_Interpolation.cpp
    MTK_Test_Proxy/cl_MTK_Field_Proxy.cpp
    UT_MTK_Periodic_Boundary_Condition_Helper.cpp
    UT_MTK_Intersection_Detect.cpp
)

# List additional includes
include_directories(${MORIS_PACKAGE_DIR}/MTK/src/stk_impl)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${ALG}-lib ${MORIS_BASE_LIBS} ${MTK}-lib)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MTK_TEST ${MTK}-test)

# Define the executable and its dependencies.
add_executable(${MTK_TEST} ${TEST_SOURCES})
target_link_libraries(${MTK_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${MTK_TEST} PROPERTIES OUTPUT_NAME ${MTK_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${MTK_TEST} ${MORIS_TEST_PROCS})
