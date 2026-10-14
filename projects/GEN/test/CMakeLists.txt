# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Geometry Engine Tests

# Sources.
set(TEST_SOURCES
    test_main.cpp
    ut_GEN_PDV.cpp
    ut_GEN_Field.cpp
    ut_GEN_Property.cpp
    ut_GEN_Intersection_Level_Set.cpp
    ut_GEN_Intersection_Surface_Mesh.cpp
    cl_GEN_Geometry_Engine_Test.cpp
    fn_GEN_create_simple_mesh.cpp
    ut_GEN_STL_Geometry.cpp
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
    ${COR}-lib
    ${ALG}-lib
    ${GEN}-lib
    ${INT}-lib
    ${DLA}-lib
    ${HMR}-lib
    ${MORIS_BASE_LIBS}
    ${PRM}-lib
)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(GEN_TEST ${GEN}-test)

# Define the executable and its dependencies.
add_executable(${GEN_TEST} ${TEST_SOURCES})
target_link_libraries(${GEN_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${GEN_TEST} PROPERTIES OUTPUT_NAME ${GEN_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${GEN_TEST} ${MORIS_TEST_PROCS})

set(SO_TPLS "trilinos" ${ARMADILLO_EIGEN})

# get list of INT subfolders
get_property(INT_SRC_LIST GLOBAL PROPERTY INT_SRC_LIST)
get_property(MTK_SRC_LIST GLOBAL PROPERTY MTK_SRC_LIST)

set(SO_INCLUDES
    ${MORIS_PACKAGE_DIR}/ALG/src
    ${INT_SRC_LIST}
    ${MTK_SRC_LIST}
    ${MORIS_PACKAGE_DIR}/FEM/MSI/src
    ${MORIS_PACKAGE_DIR}/FEM/VIS/src
    ${MORIS_PACKAGE_DIR}/SOL/DLA/src
    ${MORIS_PACKAGE_DIR}/SOL/TSA/src
    ${MORIS_PACKAGE_DIR}/SOL/NLA/src
    ${MORIS_PACKAGE_DIR}/SOL/SOL_CORE/src
    ${MORIS_PACKAGE_DIR}/LINALG/src
    ${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
    ${MORIS_PACKAGE_DIR}/COM/src
    ${MORIS_PACKAGE_DIR}/MTK/src
    ${MORIS_PACKAGE_DIR}/PRM/src
    ${MORIS_PACKAGE_DIR}/MRS/COR/src
)

foreach(TPL ${SO_TPLS})
    string(TOUPPER ${TPL} TPL)
    list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
endforeach()
