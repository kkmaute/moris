#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# HMR Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    main.cpp
    ut_HMR_Background_Mesh.cpp
    ut_HMR_Background_Mesh_Private.cpp
    ut_HMR_BSpline_Mesh.cpp
    ut_HMR_BSpline_Mesh_Private.cpp
    ut_HMR_Lagrange_Mesh.cpp
    ut_HMR_Node_Sharing.cpp
    ut_HMR_T_Matrix.cpp
    ut_HMR_SideSets.cpp
    ut_HMR_L2.cpp
    ut_HMR_Integration_Mesh.cpp
    ut_HMR_Lagrange_Elements.cpp
    ut_HMR_IO.cpp
    ut_HMR_User_Defined_Refinement.cpp)

set(TEST_DEPENDENCIES
    test-libs
    ${HMR}-lib
    ${ALG}-lib
    ${MORIS_BASE_LIBS}
    )
    
# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(HMR_TEST ${HMR}-test)

# Create executable
add_executable(${HMR_TEST} ${TEST_SOURCES})
target_link_libraries(${HMR_TEST} PRIVATE ${TEST_DEPENDENCIES})
set_target_properties(${HMR_TEST} PROPERTIES OUTPUT_NAME ${HMR_TEST}${EXE_EXT})

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(HMR_TEST_PROCS ${MORIS_TEST_PROCS_HMR})
    foreach(PROCS ${HMR_TEST_PROCS})
        set(HMR_PARALLEL_TEST ${HMR_TEST}-${PROCS}-procs)
        add_test(NAME ${HMR_PARALLEL_TEST} 
            COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS} $<TARGET_FILE:${HMR_TEST}>)
    endforeach()
endif()


