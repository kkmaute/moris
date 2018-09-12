# LINALG --------------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${LINALG}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
	cl_Matrix.hpp
	fn_comp_abs.hpp
	fn_det.hpp
	fn_eye.hpp
	fn_norm.hpp
	fn_print.hpp
	fn_trans.hpp
	fn_isempty.hpp
	fn_iscol.hpp
	fn_isrow.hpp
	fn_isvector.hpp
	fn_issquare.hpp
	fn_isfinite.hpp
	linalg_typedefs.hpp
	op_equal_equal.hpp
	op_greater_equal.hpp
	op_greater.hpp
	op_less_equal.hpp
	op_less.hpp
	op_minus.hpp
	op_not_equal.hpp
	op_plus.hpp
	op_times.hpp)
	
if(MORIS_USE_EIGEN)
list(APPEND HEADERS
	Eigen_Impl/cl_Matrix_Eigen_3x3.hpp
	Eigen_Impl/cl_Matrix_Eigen_Dynamic.hpp
	Eigen_Impl/fn_comp_abs_Eigen.hpp
 	Eigen_Impl/fn_det_Eigen.hpp
	Eigen_Impl/fn_eye_Eigen.hpp
	Eigen_Impl/fn_norm_Eigen.hpp
	Eigen_Impl/fn_trans_Eigen.hpp
	Eigen_Impl/fn_isempty_Eigen.hpp
	Eigen_Impl/fn_iscol_Eigen.hpp
	Eigen_Impl/fn_isrow_Eigen.hpp
	Eigen_Impl/fn_isfinite_Eigen.hpp
	Eigen_Impl/fn_issquare_Eigen.hpp
	Eigen_Impl/fn_isvector_Eigen.hpp
	Eigen_Impl/op_equal_equal_Eigen.hpp
	Eigen_Impl/op_greater_equal_Eigen.hpp
	Eigen_Impl/op_greater_Eigen.hpp
	Eigen_Impl/op_less_equal_Eigen.hpp
	Eigen_Impl/op_less_Eigen.hpp
	Eigen_Impl/op_minus_Eigen.hpp
	Eigen_Impl/op_not_equal_Eigen.hpp
	Eigen_Impl/op_plus_Eigen.hpp
	Eigen_Impl/op_times_Eigen.hpp)
endif()

if(MORIS_USE_ARMA)
list(APPEND HEADERS
	Arma_Impl/cl_Matrix_Arma_Dynamic.hpp
	Arma_Impl/fn_comp_abs_Arma.hpp
	Arma_Impl/fn_det_Arma.hpp
	Arma_Impl/fn_eye_Arma.hpp
	Arma_Impl/fn_norm_Arma.hpp
	Arma_Impl/fn_trans_Arma.hpp
	Arma_Impl/fn_isempty_Arma.hpp
	Arma_Impl/fn_isrow_Arma.hpp
	Arma_Impl/fn_iscol_Arma.hpp
	Arma_Impl/fn_isvector_Arma.hpp
	Arma_Impl/fn_issquare_Arma.hpp
	Arma_Impl/fn_isfinite_Arma.hpp
	Arma_Impl/op_equal_equal_Arma.hpp
	Arma_Impl/op_greater_equal_Arma.hpp
	Arma_Impl/op_greater_Arma.hpp
	Arma_Impl/op_less_equal_Arma.hpp
	Arma_Impl/op_less_Arma.hpp
	Arma_Impl/op_minus_Arma.hpp
	Arma_Impl/op_not_equal.hpp
	Arma_Impl/op_plus_Arma.hpp
	Arma_Impl/op_times_Arma.hpp)
endif()


# List header dependencies
set(LINALG_INCLUDES
    ${COM}
    ${LNA})

foreach(LINALG_INCLUDE ${LINALG_INCLUDES})
    include_directories(../../${LINALG_INCLUDE}/src)
endforeach()

# List dependencies

include_directories(${LINALG_IMPLEMENTATION_INCLUDES})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Executable

# LINALG Executable
# List source files
set(EXE_SOURCES
        main.cpp)

# List library dependencies
set(EXE_DEPENDENCIES 
	${MORIS_ARMADILLO_EIGEN_LIBS} 
	${MORIS_LDLIBS} 
	${COM}-lib
	${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_SUPERLU_LIBS})
	
# Set the output path for library
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${bin})

# Create library
add_executable(${LINALG} ${EXE_SOURCES})
target_link_libraries(${LINALG} ${EXE_DEPENDENCIES})




include_directories(${MORIS_PACKAGE_DIR}/XTK/include)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${LINALG} )
