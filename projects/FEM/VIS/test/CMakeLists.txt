# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Visualization Tests

# Sources.
set(TEST_SOURCES test_main.cpp UT_VIS_Visualization_Mesh.cpp UT_VIS_Output_Manager.cpp
                 UT_VIS_Output_Dof.cpp
)

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
    ${VIS}-lib
    ${MTK}-lib
    ${HMR}-lib
    ${XTK}-lib
    ${GEN}-lib
    ${MORIS_BASE_LIBS}
)

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(VIS_TEST ${VIS}-test)

# Define the executable and its dependencies.
add_executable(${VIS_TEST} ${TEST_SOURCES})
target_link_libraries(${VIS_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${VIS_TEST} PROPERTIES OUTPUT_NAME ${VIS_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${VIS_TEST} ${MORIS_TEST_PROCS})
