# STK Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Hierarchical_Mesh.cpp
    cl_Mesh_From_Data.cpp
    cl_Mesh_From_String.cpp
    cl_Mesh_Input_Buckets.cpp )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM}
    ${ALG}
    ${DLA}
    ${MOD}
    ${MTK}
    ${TOL}
    ${XTK}
    ${STK} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../${TEST_INCLUDE}/src)
endforeach()

# Additional includes
include_directories(../../${STK}/src/Hierarchical)

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(../../../snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${HIA}-lib
    ${MOD}-lib
    ${STK}-lib
    ${MORIS_LDLIBS}
    ${TOL}-lib )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${STK}-test ${TEST_SOURCES})
target_link_libraries(${STK}-test ${TEST_DEPENDENCIES})
