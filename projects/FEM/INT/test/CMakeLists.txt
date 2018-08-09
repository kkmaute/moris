# Integration and Interpolation Tests -------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_FEM_Integrator.cpp
    cl_FEM_Interpolation_Function_Lagrange_Hex20.cpp
    cl_FEM_Interpolation_Function_Lagrange_Hex27.cpp
    cl_FEM_Interpolation_Function_Lagrange_Hex64.cpp
    cl_FEM_Interpolation_Function_Lagrange_Hex8.cpp
    cl_FEM_Interpolation_Function_Lagrange_Quad16.cpp
    cl_FEM_Interpolation_Function_Lagrange_Quad4.cpp
    cl_FEM_Interpolation_Function_Lagrange_Quad8.cpp
    cl_FEM_Interpolation_Function_Lagrange_Quad9.cpp
    )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM}
    ${FEM}/${INT}
    ${MTK})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${INT}-lib
    ${MSI}-lib
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${INT}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${INT}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${INT}-test
    COMMAND ${INT}-test${EXE_EXT} )

