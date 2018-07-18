# Containers Tests --------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Array.cpp
    cl_Bi_Map.cpp
    cl_Cell.cpp
    cl_Dist_Map.cpp
    cl_Map.cpp
    cl_Param_List.cpp
    cl_Tuple.cpp )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM}
    ${ALG})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../../${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MORIS_LDLIBS}
    ${MORIS_BOOST_LIBS})

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${CON}-test ${TEST_SOURCES})
target_link_libraries(${CON}-test ${TEST_DEPENDENCIES})
