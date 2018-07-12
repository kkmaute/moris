# Linear Algebra Tests ----------------------------------------------------
# -------------------------------------------------------------------------

# List source files
set(TEST_SOURCES
    test_main.cpp
    cl_Col.cpp
    cl_Mat.cpp
    cl_Row.cpp
    cl_Sp_Mat.cpp
    cl_Tensor.cpp
    fn_chol_l.cpp
    fn_chol_u.cpp
    fn_cond.cpp
    fn_det.cpp
    fn_diag.cpp
    fn_eig_gen.cpp
    fn_eig_sym.cpp
    fn_eye.cpp
    fn_find.cpp
    fn_get_sparsity.cpp
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
    fn_mem_pointer.cpp
    fn_qr.cpp
    fn_reshape.cpp
    fn_save_to_AND_load_from_binary_file.cpp
    fn_sort.cpp
    fn_sum.cpp
    fn_svd.cpp
    fn_trans.cpp
    fn_unique.cpp
    op_elemwise_div.cpp
    op_elemwise_mult.cpp
    op_equal_equal.cpp
    op_greater.cpp
    op_greater_equal.cpp
    op_less.cpp
    op_less_equal.cpp
    op_minus.cpp
    op_not_equal.cpp
    op_ostream.cpp
    op_plus.cpp
    op_times.cpp
    perf_cl_Mat.cpp
    perf_cl_Sp_Mat.cpp )

# List test includes
set(TEST_INCLUDES
    ${LNA}
    ${IOS}
    ${COM}
    ${ALG})

foreach(TEST_INCLUDE ${TEST_INCLUDES})
    include_directories(../../${TEST_INCLUDE}/src)
endforeach()

# List additional includes

# List test dependencies
set(TEST_DEPENDENCIES
    ${COM}-lib
    ${MORIS_LDLIBS}
    ${ALG}-lib)

# Set the output path for test
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${LNA}-test ${TEST_SOURCES})
target_link_libraries(${LNA}-test ${TEST_DEPENDENCIES})
