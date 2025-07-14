#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# LINALG --------------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${LINALG}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
cl_Matrix.hpp
fn_all_true.hpp
fn_chol_l.hpp
fn_chol_u.hpp
fn_comp_abs.hpp
fn_cond.hpp
fn_cross.hpp
fn_ctrans.hpp
fn_det.hpp
fn_diag_mat.hpp
fn_diag_vec.hpp
fn_dot.hpp
fn_eig_gen.hpp
fn_eig_sym.hpp
fn_eye.hpp
fn_find_unique.hpp
fn_find.hpp
fn_histc.hpp
fn_interp1.hpp
fn_intersect.hpp
fn_inv.hpp
fn_iscol.hpp
fn_isempty.hpp
fn_isfinite.hpp
fn_isrow.hpp
fn_issquare.hpp
fn_isvector.hpp
fn_linsolve.hpp
fn_linspace.hpp
fn_load_matrix_from_binary_file.hpp
fn_mat_to_vec.hpp
fn_lu.hpp
fn_max.hpp
fn_min.hpp
fn_norm.hpp
fn_print.hpp
fn_prod.hpp
fn_qr.hpp
fn_r2.hpp
fn_reshape.hpp
fn_save_matrix_to_binary_file.hpp
fn_sort.hpp
fn_sort_index.hpp
fn_sqrtmat.hpp
fn_sum.hpp
fn_svd.hpp
fn_sylvester.hpp
fn_trans.hpp
fn_unique.hpp
fn_vectorize.hpp
linalg_typedefs.hpp
op_div.hpp
op_elemwise_div.hpp
op_elemwise_mult.hpp
op_equal_equal.hpp
op_greater_equal.hpp
op_greater.hpp
op_less_equal.hpp
op_less.hpp
op_minus.hpp
op_move.hpp
op_not_equal.hpp
op_plus.hpp
op_times.hpp
fn_rank.hpp
fn_join_horiz.hpp
fn_join_vert.hpp
)
    
if(MORIS_USE_EIGEN)
list(APPEND HEADERS
Eigen_Impl/cl_Matrix_Eigen_3x1.hpp
Eigen_Impl/cl_Matrix_Eigen_3x3.hpp
Eigen_Impl/cl_Matrix_Eigen_Dynamic.hpp
Eigen_Impl/fn_chol_l_Eigen.hpp
Eigen_Impl/fn_chol_u_Eigen.hpp
Eigen_Impl/fn_comp_abs_Eigen.hpp
Eigen_Impl/fn_cond_Eigen.hpp
Eigen_Impl/fn_cross_Eigen.hpp
Eigen_Impl/fn_ctrans_Eigen.hpp
Eigen_Impl/fn_det_Eigen.hpp
Eigen_Impl/fn_diag_mat_Eigen.hpp
Eigen_Impl/fn_diag_vec_Eigen.hpp
Eigen_Impl/fn_dot_Eigen.hpp
Eigen_Impl/fn_eig_gen_Eigen.hpp
Eigen_Impl/fn_eig_sym_Eigen.hpp
Eigen_Impl/fn_eye_Eigen.hpp
Eigen_Impl/fn_find_Eigen.hpp
Eigen_Impl/fn_find_unique_Eigen.hpp
Eigen_Impl/fn_histc_Eigen.hpp
Eigen_Impl/fn_inv_Eigen.hpp
Eigen_Impl/fn_iscol_Eigen.hpp
Eigen_Impl/fn_isempty_Eigen.hpp
Eigen_Impl/fn_isfinite_Eigen.hpp
Eigen_Impl/fn_isrow_Eigen.hpp
Eigen_Impl/fn_issquare_Eigen.hpp
Eigen_Impl/fn_isvector_Eigen.hpp
Eigen_Impl/fn_linsolve_Eigen.hpp
Eigen_Impl/fn_linspace_Eigen.hpp
Eigen_Impl/fn_lu_Eigen.hpp
Eigen_Impl/fn_norm_Eigen.hpp
Eigen_Impl/fn_qr_Eigen.hpp
Eigen_Impl/fn_reshape_Eigen.hpp
Eigen_Impl/fn_sort_Eigen.hpp
Eigen_Impl/fn_sort_index_Eigen.hpp
Eigen_Impl/fn_sum_Eigen.hpp
Eigen_Impl/fn_svd_Eigen.hpp
Eigen_Impl/fn_trans_Eigen.hpp
Eigen_Impl/fn_unique_Eigen.hpp
Eigen_Impl/fn_vectorize_Eigen.hpp
Eigen_Impl/op_elemwise_div_Eigen.hpp
Eigen_Impl/op_elemwise_mult_Eigen.hpp
Eigen_Impl/op_equal_equal_Eigen.hpp
Eigen_Impl/op_greater_Eigen.hpp
Eigen_Impl/op_greater_equal_Eigen.hpp
Eigen_Impl/op_less_Eigen.hpp
Eigen_Impl/op_less_equal_Eigen.hpp
Eigen_Impl/op_minus_Eigen.hpp
Eigen_Impl/op_not_equal_Eigen.hpp
Eigen_Impl/op_plus_Eigen.hpp
Eigen_Impl/op_times_Eigen.hpp)
endif()

if(MORIS_USE_ARMA)
list(APPEND HEADERS
Arma_Impl/cl_Matrix_Arma_Dynamic.hpp
Arma_Impl/fn_chol_l_Arma.hpp
Arma_Impl/fn_chol_u_Arma.hpp
Arma_Impl/fn_comp_abs_Arma.hpp
Arma_Impl/fn_cond_Arma.hpp
Arma_Impl/fn_cross_Arma.hpp
Arma_Impl/fn_ctrans_Arma.hpp
Arma_Impl/fn_det_Arma.hpp
Arma_Impl/fn_diag_mat_Arma.hpp
Arma_Impl/fn_diag_vec_Arma.hpp
Arma_Impl/fn_dot_Arma.hpp
Arma_Impl/fn_eig_gen_Arma.hpp
Arma_Impl/fn_eig_sym_Arma.hpp
Arma_Impl/fn_eye_Arma.hpp
Arma_Impl/fn_find_Arma.hpp
Arma_Impl/fn_find_unique_Arma.hpp
Arma_Impl/fn_histc_Arma.hpp
Arma_Impl/fn_interp1_Arma.hpp
Arma_Impl/fn_intersect_Arma.hpp
Arma_Impl/fn_inv_Arma.hpp
Arma_Impl/fn_iscol_Arma.hpp
Arma_Impl/fn_isempty_Arma.hpp
Arma_Impl/fn_isfinite_Arma.hpp
Arma_Impl/fn_isrow_Arma.hpp
Arma_Impl/fn_issquare_Arma.hpp
Arma_Impl/fn_isvector_Arma.hpp
Arma_Impl/fn_linsolve_Arma.hpp
Arma_Impl/fn_linspace_Arma.hpp
Arma_Impl/fn_lu_Arma.hpp
Arma_Impl/fn_norm_Arma.hpp
Arma_Impl/fn_qr_Arma.hpp
Arma_Impl/fn_reshape_Arma.hpp
Arma_Impl/fn_sort_Arma.hpp
Arma_Impl/fn_sort_index_Arma.hpp
Arma_Impl/fn_sqrtmat_Arma.hpp
Arma_Impl/fn_sum_Arma.hpp
Arma_Impl/fn_svd_Arma.hpp
Arma_Impl/fn_sylvester_Arma.hpp
Arma_Impl/fn_trans_Arma.hpp
Arma_Impl/fn_unique_Arma.hpp
Arma_Impl/fn_vectorize_Arma.hpp
Arma_Impl/op_elemwise_div_Arma.hpp
Arma_Impl/op_elemwise_mult_Arma.hpp
Arma_Impl/op_equal_equal_Arma.hpp
Arma_Impl/op_greater_Arma.hpp
Arma_Impl/op_greater_equal_Arma.hpp
Arma_Impl/op_less_Arma.hpp
Arma_Impl/op_less_equal_Arma.hpp
Arma_Impl/op_minus_Arma.hpp
Arma_Impl/op_not_equal_Arma.hpp
Arma_Impl/op_plus_Arma.hpp
Arma_Impl/op_times_Arma.hpp
Arma_Impl/fn_rank_Arma.hpp
Arma_Impl/fn_join_horiz_Arma.hpp
Arma_Impl/fn_join_vert_Arma.hpp
Arma_Impl/fn_max_Arma.hpp
Arma_Impl/fn_min_Arma.hpp
  )
endif()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Library
set(LIB_DEPENDENCIES
	${ASR}-lib
	${CNT}-lib
	)

add_library(${LINALG}-lib INTERFACE)
target_include_directories(${LINALG}-lib INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${LINALG_IMPLEMENTATION_INCLUDES}>
	$<INSTALL_INTERFACE:${${LINALG}_HEADER_INSTALL_DIR}>
	)
target_link_libraries(${LINALG}-lib INTERFACE ${LIB_DEPENDENCIES})

foreach(TPL ${LINALG_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${LINALG}-lib
    	INTERFACE ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${LINALG}-lib
    	INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${LINALG}-lib
    	INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${LINALG_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND LINALG_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	# skip if -lfoo format
	string(FIND ${LIB_TARGET} "-" pos)
	if(pos EQUAL "0")
		continue()
	endif()
	
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	#string(APPEND LINALG_CONFIG_MORIS_DEPENDENCIES
	#	"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )

        string(APPEND LINALG_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} REQUIRED HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${LINALG}Config.cmake.in
    ${${LINALG}_CONFIG_BUILD_DIR}/${LINALG}Config.cmake
    INSTALL_DESTINATION 
    ${${LINALG}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${LINALG}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${LINALG}_CONFIG_BUILD_DIR}/${LINALG}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Installation

# Install *Config and *ConfigVersion files
install(FILES ${${LINALG}_CONFIG_BUILD_DIR}/${LINALG}Config.cmake 
    ${${LINALG}_CONFIG_BUILD_DIR}/${LINALG}ConfigVersion.cmake
    DESTINATION ${${LINALG}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${LINALG}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${LINALG}-lib
	EXPORT ${LINALG}Targets
    DESTINATION ${${LINALG}_LIB_INSTALL_DIR} )

# Export library into *Targets file
install(EXPORT ${LINALG}Targets
    DESTINATION ${${LINALG}_CONFIG_INSTALL_DIR} )
