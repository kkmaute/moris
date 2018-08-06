# Install script for directory: /home/schmidt/codes/moris/projects/LNA/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/LNA" TYPE FILE FILES
    "/home/schmidt/codes/moris/projects/LNA/src/linalg.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Arma_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Arma_Mat.tpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Arma_Sp_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Arma_Sp_Mat.tpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Base_Arma_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Base_Eigen_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Base_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Eigen_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Eigen_Mat.tpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Eigen_Sp_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Eigen_Sp_Mat.tpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Sp_Mat.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_Tensor.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/cl_TensorMapCreator.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_chol_l.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_chol_u.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_cond.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_ctrans.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_det.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_diag.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_dot.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_eig_gen.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_eig_sym.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_eye.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_find.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_find_unique.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_get_sparsity.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_histc.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_inv.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_iscol.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_isempty.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_isfinite.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_isrow.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_issquare.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_isvector.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_linsolve.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_linspace.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_load_vector_from_binary_file.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_lu.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_mem_pointer.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_pinv.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_prod.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_qr.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_reshape.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_save_vector_to_binary_file.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_sort.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_sum.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_svd.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_trans.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/fn_unique.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_elemwise_div.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_elemwise_mult.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_equal_equal.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_greater_equal.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_greater.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_less_equal.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_less.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_minus.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_not_equal.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_ostream.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_plus.hpp"
    "/home/schmidt/codes/moris/projects/LNA/src/op_times.hpp"
    )
endif()

