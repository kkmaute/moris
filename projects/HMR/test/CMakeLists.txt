# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# HMR Tests

# Sources.
set(TEST_SOURCES
    main.cpp
    ut_HMR_Background_Mesh.cpp
    ut_HMR_Background_Mesh_Private.cpp
    ut_HMR_BSpline.cpp
    ut_HMR_BSpline_Mesh.cpp
    ut_HMR_BSpline_Mesh_Private.cpp
    ut_HMR_Lagrange_Mesh.cpp
    ut_HMR_T_Matrix.cpp
    ut_HMR_SideSets.cpp
    ut_HMR_L2.cpp
    ut_HMR_Integration_Mesh.cpp
    ut_HMR_Lagrange_Elements.cpp
    ut_HMR_IO.cpp
    ut_HMR_User_Defined_Refinement.cpp
)

set(TEST_DEPENDENCIES test-libs ${HMR}-lib ${ALG}-lib ${MORIS_BASE_LIBS})

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(HMR_TEST ${HMR}-test)

# Define the executable and its dependencies.
add_executable(${HMR_TEST} ${TEST_SOURCES})
target_link_libraries(${HMR_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${HMR_TEST} PROPERTIES OUTPUT_NAME ${HMR_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${HMR_TEST} ${MORIS_TEST_PROCS_HMR})
