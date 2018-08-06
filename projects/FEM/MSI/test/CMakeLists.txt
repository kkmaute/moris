# MSI ----------------------------------------
# -------------------------------------------------------------------------
message("Hello!")
# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Dof_Manager_Test.cpp
    cl_Equation_Object_Test.cpp
    cl_Pdof_Host_Test.cpp
    cl_Model_Solver_Interface_Test.cpp)

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM}
    ${ALG}
    ${MSI}
    ${DLA})

include_directories(${MORIS_INCDIRS})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../../${TEST_INCLUDE}/src)
endforeach()

# List additional includes
include_directories(../src)
include_directories(../../../LNA/src)
include_directories(../../../IOS/src)
include_directories(../../../DLA/src)

# List test dependencies
set(TEST_DEPENDENCIES
    ${DLA}-lib
    ${IOS}-lib
    ${ALG}-lib
    ${COM}-lib
    ${MSI}-lib
    ${MORIS_LDLIBS}
    ${MORIS_MATH_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${MSI}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MSI}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MSI}-test
    COMMAND ${MSI}-test${EXE_EXT} )




