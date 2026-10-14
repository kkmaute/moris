# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# SDF Tests

# Sources.
set(TEST_SOURCES test_main.cpp ut_SDF_Line.cpp ut_SDF_Triangle.cpp ut_SDF_Raycast.cpp
                 ut_SDF_Core.cpp ut_SDF_Object.cpp
)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${SDF}-lib ${MTK}-lib ${MORIS_BASE_LIBS})

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(SDF_TEST ${SDF}-test)

# Define the executable and its dependencies.
add_executable(${SDF_TEST} ${TEST_SOURCES})
target_link_libraries(${SDF_TEST} ${TEST_DEPENDENCIES})
set_target_properties(${SDF_TEST} PROPERTIES OUTPUT_NAME ${SDF_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${SDF_TEST} ${MORIS_TEST_PROCS})
