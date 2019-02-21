# Geometry Engine Tests ---------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_GeometryEngineFactory.cpp
    cl_GeometryObject.cpp
    cl_GeSDFCore.cpp
    cl_GeSDFGen.cpp
    cl_GeTriangle.cpp
    cl_GE_test.cpp
    )

# List test includes
set(TEST_INCLUDES
    ${MTK}
    ${MTK}/${MAP}
    ${XTK} #MTK
    ${COM}
    ${ALG}
    ${GEN}
    ${GEN}/${SDF}
    ${FEM}/${INT}
    ${FEM}/${MSI}
    ${FEM}/${MDL}
    ${HMR}
    ${MTK}
    ${LINALG}
)

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List snippet includes
set(SNIPPET_INCLUDES
    )

foreach(SNIPPET_INCLUDE ${SNIPPET_INCLUDES})
    include_directories(${MORIS_DIR}/snippets/${SNIPPET_INCLUDE})
endforeach()

# List test dependencies
set(TEST_DEPENDENCIES
    ${HMR}-lib
    ${MAP}-lib
    ${MTK}-lib
	${SDF}-lib
	${ALG}-lib
    ${MTK}-lib
    ${COM}-lib
    ${MDL}-lib # for L2
    ${INT}-lib # for L2
    ${MSI}-lib # for L2
    ${NLA}-lib # for L2
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_BASE_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(GEN_TEST ${GEN}-test)

# Create executable
add_executable(${GEN_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${GEN_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${GEN_TEST}
    COMMAND ${GEN_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(GEN_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${GEN_TEST_PROCS})
		set(GEN_PARALLEL_TEST ${GEN_TEST}-${PROCS}-procs)
		add_test(NAME ${GEN_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${GEN_TEST}${EXE_EXT}>)
	endforeach()
endif()

