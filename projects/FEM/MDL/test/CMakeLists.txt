# Integration and Interpolation Tests -------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    cl_MDL_DiffusionElement_Mesh_Interaction_Test.cpp
    cl_MDL_DiffusionElement_Cut_Mesh.cpp
    test_main.cpp
    )

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${IOS}
    ${COM}
    ${MTK}
    ${HMR}
    ${ALG}
    ${SOL}/${DLA}
    ${SOL}/${NLA}
    ${SOL}/${TSA}
    ${SOL}/${SOL_CORE}
    ${FEM}/${INT}
    ${FEM}/${MSI}
    ${FEM}/${MDL}
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MDL}-lib
    ${MSI}-lib
    ${INT}-lib
    ${MTK}-lib
    ${TOL}-lib
    ${HMR}-lib
    ${MAP}-lib
    ${DLA}-lib
    ${NLA}-lib
    ${TSA}-lib
    ${SOL_CORE}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MDL_TEST ${MDL}-test)

# Create executable
add_executable(${MDL_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MDL_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MDL_TEST}
    COMMAND ${MDL_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
    set(MDL_TEST_PROCS ${MORIS_TEST_PROCS})
    foreach(PROCS ${MDL_TEST_PROCS})
        set(MDL_PARALLEL_TEST ${MDL_TEST}-${PROCS}-procs)
        add_test(NAME ${MDL_PARALLEL_TEST} 
            COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${MDL_TEST}${EXE_EXT}>)
    endforeach()
endif()

