# Linear Algebra Source ---------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    linalg.hpp
    cl_Arma_Mat.hpp
    cl_Arma_Mat.tpp
    cl_Arma_Sp_Mat.hpp
    cl_Arma_Sp_Mat.tpp
    cl_Base_Arma_Mat.hpp
    cl_Base_Eigen_Mat.hpp
    cl_Base_Mat.hpp
    cl_Eigen_Mat.hpp
    cl_Eigen_Mat.tpp
    cl_Eigen_Sp_Mat.hpp
    cl_Eigen_Sp_Mat.tpp
    cl_Mat.hpp
    cl_Sp_Mat.hpp
    cl_Tensor.hpp
    cl_TensorMapCreator.hpp
    fn_chol_l.hpp
    fn_chol_u.hpp
    fn_cond.hpp
    fn_ctrans.hpp
    fn_det.hpp
    fn_diag.hpp
    fn_dot.hpp
    fn_eig_gen.hpp
    fn_eig_sym.hpp
    fn_eye.hpp
    fn_find.hpp
    fn_find_unique.hpp
    fn_get_sparsity.hpp
    fn_histc.hpp
    fn_inv.hpp
    fn_iscol.hpp
    fn_isempty.hpp
    fn_isfinite.hpp
    fn_isrow.hpp
    fn_issquare.hpp
    fn_isvector.hpp
    fn_linsolve.hpp
    fn_linspace.hpp
    fn_load_vector_from_binary_file.hpp
    fn_lu.hpp
    fn_mem_pointer.hpp
    fn_norm.hpp
    fn_pinv.hpp
    fn_prod.hpp
    fn_qr.hpp
    fn_reshape.hpp
    fn_save_vector_to_binary_file.hpp
    fn_sort.hpp
    fn_sum.hpp
    fn_svd.hpp
    fn_trans.hpp
    fn_unique.hpp
    op_elemwise_div.hpp
    op_elemwise_mult.hpp
    op_equal_equal.hpp
    op_greater_equal.hpp
    op_greater.hpp
    op_less_equal.hpp
    op_less.hpp
    op_minus.hpp
    op_not_equal.hpp
    op_ostream.hpp
    op_plus.hpp
    op_times.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${LNA} )
