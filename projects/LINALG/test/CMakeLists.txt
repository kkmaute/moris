# Linear Algebra Tests --------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
	test_main.cpp
	cl_Matrix.cpp
	fn_chol_l.cpp
	fn_chol_u.cpp
	fn_comp_abs.cpp
	fn_ctrans.cpp 
	fn_det.cpp
	fn_dot.cpp
	fn_eye.cpp
	fn_inv.cpp
	fn_det.cpp
	fn_diag.cpp
	fn_eye.cpp
	fn_norm.cpp
	fn_trans.cpp
	fn_isempty.cpp
	fn_isrow.cpp
	fn_iscol.cpp
	fn_isfinite.cpp
	fn_issquare.cpp
	fn_isvector.cpp
	LINALG_Tutorial.cpp
	fn_norm.cpp
    fn_qr.cpp
    fn_r2.cpp
    fn_reshape.cpp
    fn_save_and_load_matrix_binary_file.cpp
    fn_sort.cpp
    fn_trans.cpp
    fn_unique.cpp
	op_elemwise_div.cpp
	op_elemwise_mult.cpp
	op_equal_equal.cpp
	op_greater_equal.cpp
	op_greater.cpp
	op_less_equal.cpp
	op_less.cpp
	op_minus.cpp
	op_not_equal.cpp 
	op_plus.cpp 
	op_times.cpp
	op_div.cpp)

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

# Create executable
add_executable(${LINALG}-test${EXE_EXT} ${TEST_SOURCES})
target_link_libraries(${LINALG}-test${EXE_EXT} ${TEST_DEPENDENCIES})

# Mark as test
add_test(NAME ${LINALG}-test
    COMMAND ${LINALG}-test${EXE_EXT} )
