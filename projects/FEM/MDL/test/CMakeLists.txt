# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Integration and Interpolation Tests
include_directories(DLA_Test_Proxys)

# Sources.
set(TEST_SOURCES
    cl_MDL_DiffusionElement_Mesh_Interaction_Test.cpp
    cl_MDL_DiffusionElement_Cut_Mesh.cpp
    #UT_MDL_Mesh_Model_Helper.cpp
    UT_MDL_XTK_DiffusionElement.cpp
    UT_MDL_XTK_HMR_DiffusionElement.cpp
    UT_MDL_XTK_HMR_Linear_Struc_2D.cpp
    UT_MDL_XTK_HMR_Linear_Struc_Contact_2D.cpp
    UT_MDL_XTK_HMR_Linear_Struc_Incompressible.cpp
    #UT_MDL_XTK_HMR_Linear_Struc_Thick_Wall_PV.cpp
    UT_MDL_XTK_HMR_ThermoElastic_2D.cpp
    UT_MDL_XTK_HMR_2_Material_Bar_Plane.cpp
    UT_MDL_XTK_HMR_4_Material_Bar_Plane_Hole.cpp
    UT_MDL_XTK_HMR_Material_Void_Bar_Plane.cpp
    UT_MDL_XTK_HMR_Symm_BCs.cpp
    #UT_MDL_Gyroid.cpp
    UT_MDL_FEM_Benchmark.cpp
    UT_MDL_FEM_Benchmark2.cpp
    UT_MDL_FEM_DQ_Dp.cpp
    UT_MDL_Fluid_Benchmark.cpp
    UT_XFEM_Measure.cpp
    #UT_MDL_Sensitivity_Test.cpp
    #UT_MDL_Transient.cpp
    test_main.cpp
    #MDL_Test_Proxys/cl_MDL_Element_Proxy.hpp
    #MDL_Test_Proxys/cl_MDL_Node_Proxy.hpp
)

set(GEN_INCLUDES "ADV" "ADV/fields" "PDV")

foreach(GEN_INCLUDES ${GEN_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
endforeach()

# Include WRK directory
include_directories(${MORIS_PACKAGE_DIR}/WRK/src)

# Test dependencies.
set(TEST_DEPENDENCIES
    test-libs
    ${MTK}-lib
    ${ALG}-lib
    ${INT}-lib
    ${HMR}-lib
    ${SDF}-lib
    ${MAP}-lib
    ${XTK}-lib
    ${GEN}-lib
    ${MORIS_BASE_LIBS}
    ${PRM}-lib
)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MDL_TEST ${MDL}-test)

# Define the executable and its dependencies.
add_executable(${MDL_TEST} ${TEST_SOURCES})
target_link_libraries(${MDL_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${MDL_TEST} PROPERTIES OUTPUT_NAME ${MDL_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${MDL_TEST} ${MORIS_TEST_PROCS})
