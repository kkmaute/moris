# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Containers Tests

# Sources.
set(TEST_SOURCES
    test_main.cpp
    cl_Array.cpp
    cl_Bi_Map.cpp
    ut_Vector.cpp
    cl_Dist_Map.cpp
    cl_Map.cpp
    ut_Parameter_List.cpp
    cl_Tuple.cpp
)

# Test dependencies.
set(TEST_DEPENDENCIES test-libs ${CNT}-lib ${COM}-lib ${ALG}-lib ${MORIS_BASE_LIBS})

# Test output directory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(CNT_TEST ${CNT}-test)

# Define the executable and its dependencies.
add_executable(${CNT_TEST} ${TEST_SOURCES})
target_link_libraries(${CNT_TEST} PUBLIC ${TEST_DEPENDENCIES})
set_target_properties(${CNT_TEST} PROPERTIES OUTPUT_NAME ${CNT_TEST}${EXE_EXT})

# Register tests for the configured MPI process counts.
moris_register_unit_tests(${CNT_TEST} ${MORIS_TEST_PROCS})
