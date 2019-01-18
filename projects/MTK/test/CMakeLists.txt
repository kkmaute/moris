# MTK Tests ---------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
	cl_MTK_Mesh.cpp
	cl_MTK_Exodus_IO_Helper.cpp
	cl_MTK_Mesh_From_Data.cpp
	cl_MTK_Vertex.cpp
	cl_MTK_Cell.cpp
	cl_MTK_Facet_Cluster.cpp)

# List test includes
set(TEST_INCLUDES
    ${LINALG}
    ${COM}
    ${ALG}
    ${SOL}/${DLA}
    ${MTK}
    ${TOL} )

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()

# List additional includes
include_directories(${MORIS_PACKAGE_DIR}/MTK/src/stk_impl)

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MTK}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    )

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(MTK_TEST ${MTK}-test)

# Create executable
add_executable(${MTK_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${MTK_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${MTK_TEST}
    COMMAND ${MTK_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(MTK_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${MTK_TEST_PROCS})
		set(MTK_PARALLEL_TEST ${MTK_TEST}-${PROCS}-procs)
		add_test(NAME ${MTK_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${MTK_TEST}${EXE_EXT}>)
	endforeach()
endif()