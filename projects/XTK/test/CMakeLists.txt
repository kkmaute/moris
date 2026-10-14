# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# XTK Tests

# Sources.

set(TEST_SOURCES
    main.cpp
    #xtk/UT_fn_flood_fill.cpp
    xtk/UT_MTK_Cell_Cluster_XTK_Impl.cpp
    xtk/UT_Snapping.cpp
    xtk/UT_XTK_Child_Mesh_Conform_2D.cpp
    xtk/UT_XTK_Child_Mesh_NH_Permutations.cpp
    xtk/UT_XTK_Child_Mesh_RegSub_2D.cpp
    xtk/UT_XTK_Cut_Mesh_Modification.cpp
    xtk/UT_XTK_Cut_Mesh_RegSub.cpp
    xtk/UT_XTK_Cut_Mesh.cpp
    xtk/UT_XTK_Downward_Inheritance.cpp
    xtk/UT_XTK_Enrichment_2D.cpp
    xtk/UT_XTK_Enrichment.cpp
    xtk/UT_XTK_Ghost_Stabilization.cpp
    xtk/UT_XTK_HMR_2D.cpp
    xtk/UT_XTK_HMR_Enrich_Inclusion.cpp
    xtk/UT_XTK_HMR_Enrich.cpp
    #xtk/UT_XTK_HMR_Multigrid.cpp
    xtk/UT_XTK_HMR.cpp
    #xtk/UT_XTK_Intersection_Detect_HMR_XTK.cpp
    xtk/UT_XTK_Mesh_Cleanup.cpp
    xtk/UT_XTK_Model_2D.cpp
    xtk/UT_XTK_Model_Mesh_Clusters.cpp
    xtk/UT_XTK_Model.cpp
    #xtk/UT_XTK_Multigrid.cpp
    xtk/UT_XTK_Octree_Refinement.cpp
    xtk/UT_XTK_Unenrichment.cpp
    xtk/UT_XTK_Basis_Extention.cpp
)

# List additional includes
include_directories(../src)
include_directories(../include)

set(XTK_INCLUDES
    "assert"
    "containers"
    "core"
    "ios"
    "linalg"
    "mesh"
    "tools"
    "topology"
    "xtk"
)
# Include dependency directories.
foreach(XTK_INCLUDE ${XTK_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/XTK/src/${XTK_INCLUDE})
endforeach()

# GEN directories
set(GEN_INCLUDES "ADV" "ADV/fields" "PDV")

foreach(GEN_INCLUDES ${GEN_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
endforeach()

# Include WRK directory
include_directories(${MORIS_PACKAGE_DIR}/WRK/src)

# Test dependencies.
set(TEST_DEPENDENCIES
    test-libs
    ${ALG}-lib
    ${MTK}-lib
    ${XTK}-lib
    ${MAP}-lib
    ${MDL}-lib
    ${INT}-lib
    ${HMR}-lib
    ${PRM}-lib
    ${GEN}-lib
    ${MORIS_BASE_LIBS}
)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(XTK_TEST ${XTK}-test)

# Define the executable and its dependencies.
add_executable(${XTK_TEST} ${TEST_SOURCES})
target_link_libraries(${XTK_TEST} PRIVATE ${TEST_DEPENDENCIES} ${HMR_TEST_DEPENDENCIES})
set_target_properties(${XTK_TEST} PROPERTIES OUTPUT_NAME ${XTK_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${XTK_TEST} ${MORIS_TEST_PROCS})
