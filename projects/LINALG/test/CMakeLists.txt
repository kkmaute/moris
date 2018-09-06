# Linear Algebra Tests --------------------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
	test_main.cpp
	cl_Matrix.cpp
	fn_comp_abs.cpp
	fn_ctrans.cpp 
	fn_det.cpp
	fn_dot.cpp
	fn_inv.cpp
	fn_det.cpp
	fn_eye.cpp
	fn_norm.cpp
	fn_trans.cpp
	fn_isempty.cpp
	fn_isrow.cpp
	fn_iscol.cpp
	fn_isfinite.cpp
	fn_issquare.cpp
	fn_isvector.cpp
	op_equal_equal.cpp
	op_greater_equal.cpp
	op_greater.cpp
	op_less.cpp
	op_minus.cpp 
	op_plus.cpp 
	op_times.cpp)

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
	${MORIS_LDLIBS}
	${ALG}-lib
	${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
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
