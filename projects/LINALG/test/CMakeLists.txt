# Linear Algebra Tests --------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
	cl_Matrix.cpp
	cl_Matrix3x1.cpp
	fn_all_true.cpp
	fn_chol_l.cpp
	fn_chol_u.cpp
	fn_comp_abs.cpp
	fn_cond.cpp
	fn_ctrans.cpp
	fn_det.cpp
	fn_diag.cpp
	fn_dot.cpp
	fn_eig_gen.cpp
	fn_eig_sym.cpp
	fn_eye.cpp
	fn_find_unique.cpp
	fn_find.cpp
	fn_histc.cpp
	fn_inv.cpp
	fn_iscol.cpp
	fn_isempty.cpp
	fn_isfinite.cpp
	fn_isrow.cpp
	fn_issquare.cpp
	fn_isvector.cpp
	fn_linsolve.cpp
	fn_linspace.cpp
	fn_lu.cpp
	fn_norm.cpp
	fn_qr.cpp
	fn_r2.cpp
	fn_reshape.cpp
	fn_save_and_load_matrix_binary_file.cpp
	fn_sort.cpp
	fn_sum.cpp
	fn_svd.cpp
	fn_trans.cpp
	fn_unique.cpp
	LINALG_Tutorial.cpp
	op_div.cpp
	op_elemwise_div.cpp
	op_elemwise_mult.cpp
	op_equal_equal.cpp
	op_greater_equal.cpp
	op_greater.cpp
	op_less_equal.cpp
	op_less.cpp
	op_minus.cpp
	op_move.cpp
	op_not_equal.cpp
	op_plus.cpp
	op_times.cpp
	test_main.cpp)

# List test includes
set(TEST_INCLUDES
    ${COM}
    ${LINALG}
    ${ALG})
    
    
foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../${TEST_INCLUDE}/src)
endforeach()

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TEST_INCLUDE}/src)
endforeach()


# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
	${COM}-lib
	
# 	${MORIS_TRILINOS_LIBS}
	
	${ALG}-lib
	${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_BASE_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_SUPERLU_LIBS}
	)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Set test name
set(LINALG_TEST ${LINALG}-test)

# Create executable
add_executable(${LINALG_TEST}${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${LINALG_TEST}${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${LINALG_TEST}
    COMMAND ${LINALG_TEST}${EXE_EXT} )

# Add parallel tests
if(MORIS_HAVE_PARALLEL_TESTS)
	set(LINALG_TEST_PROCS ${MORIS_TEST_PROCS})
	foreach(PROCS ${LINALG_TEST_PROCS})
		set(LINALG_PARALLEL_TEST ${LINALG_TEST}-${PROCS}-procs)
		add_test(NAME ${LINALG_PARALLEL_TEST} 
			COMMAND mpirun -n ${PROCS} $<TARGET_FILE:${LINALG_TEST}${EXE_EXT}>)
	endforeach()
endif()
