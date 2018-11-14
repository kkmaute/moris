# MSI ----------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_MSI_Dof_Manager_Test.cpp
    cl_MSI_Equation_Object_Test.cpp
    cl_MSI_Pdof_Host_Test.cpp
    cl_MSI_Model_Solver_Interface_Test.cpp
    )

# List test includes
set(TEST_INCLUDES
    ${DLA}
    ${NLA}
    ${COM}
    ${LINALG}
    ${ALG}
    ${MTK}
    ${FEM}/${INT}
    ${FEM}/${MSI}
)

include_directories(${MORIS_INCDIRS})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
#    ${NLA}-lib
    ${INT}-lib
    ${MSI}-lib
    ${COM}-lib
    ${MORIS_BASE_LIBS}
#    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MSI_TEST ${MSI}-test)

# Create executable
add_executable(${MSI_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MSI_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MSI_TEST}
    COMMAND ${MSI_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(MSI_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${MSI_TEST_PROCS})
		set(MSI_PARALLEL_TEST ${MSI_TEST}-${PROCS}-procs)
		add_test(NAME ${MSI_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${MSI_TEST}${EXE_EXT}>)
	endforeach()
endif()




