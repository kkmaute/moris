# HMR Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    main.cpp
    cl_HMR_Background_Mesh.cpp
    cl_HMR_Background_Mesh_Private.cpp
    cl_HMR_BSpline_Mesh.cpp
    cl_HMR_BSpline_Mesh_Private.cpp
    cl_HMR_Lagrange_Mesh.cpp
    cl_HMR_SideSets.cpp
    cl_HMR_L2.cpp )
  
# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${HMR}
    ${MTK}
    ${ALG} #MTK
    ${TOL} #MTK
    ${XTK} #MTK
    )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
#set(TEST_DEPENDENCIES
#    ${MORIS_ARMADILLO_EIGEN_LIBS}
#    ${MORIS_SUPERLU_LIBS}
#    ${MORIS_TRILINOS_LIBS}
#    ${MORIS_BASE_LIBS}
#    ${MORIS_ARPACK_LIBS}
#    ${MORIS_BOOST_LIBS}
#    ${COM}-lib
#    ${HMR}-lib
#    ${MDL}-lib
#    ${MSI}-lib
#    ${INT}-lib
#    ${MAP}-lib
#    ${MTK}-lib
#    ${GEN}-lib
#    ${TOL}-lib
#    )

set(TEST_DEPENDENCIES
    ${HMR}-lib
    ${MAP}-lib
    ${MTK}-lib
    ${MDL}-lib # for L2
    ${INT}-lib # for L2
    ${MSI}-lib # for L2
    ${NLA}-lib # for L2
    ${SDF}-lib # for SDF
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    )
    
# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(HMR_TEST ${HMR}-test)

# Create executable
add_executable(${HMR_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${HMR_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${HMR_TEST}
    COMMAND ${HMR_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(HMR_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${HMR_TEST_PROCS})
		set(HMR_PARALLEL_TEST ${HMR_TEST}-${PROCS}-procs)
		add_test(NAME ${HMR_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${HMR_TEST}${EXE_EXT}>)
	endforeach()
endif()

