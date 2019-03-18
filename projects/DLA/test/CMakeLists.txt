# Distributed Linear Algebra Tests ----------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Dist_Vector_Test.cpp #> VectorType error
    cl_Linear_Solver_Test.cpp
    cl_Map_Epetra.cpp
    cl_Sparse_Matrix_Test.cpp
    cl_DLA_Multigrid_Test.cpp
    map.mtx
    matrix.mtx
    vector.mtx )

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${DLA}
    ${ALG}
    ${HMR}
    ${MTK}
    ${NLA}
    ${GEN}/${SDF}
    ${MTK}/${MAP}
    ${FEM}/${MSI}
    ${FEM}/${INT} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${DLA}-lib
    ${HMR}-lib
    ${MAP}-lib
    ${MTK}-lib
    ${MDL}-lib 
    ${INT}-lib 
    ${MSI}-lib 
    ${NLA}-lib 
    ${SDF}-lib 
    ${MORIS_ARMADILLO_EIGEN_LIBS} )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(DLA_TEST ${DLA}-test)

# Create executable
add_executable(${DLA_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${DLA_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${DLA_TEST}
    COMMAND ${DLA_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(DLA_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${DLA_TEST_PROCS})
		set(DLA_PARALLEL_TEST ${DLA_TEST}-${PROCS}-procs)
		add_test(NAME ${DLA_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${DLA_TEST}${EXE_EXT}>)
	endforeach()
endif()

