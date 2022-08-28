/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_diag_vec_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DIAG_VEC_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DIAG_VEC_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "fn_isvector.hpp"
#include "Eigen/Dense"

namespace moris
{
template< typename ET>
auto
diag_vec( Eigen::MatrixBase<ET> & aA,
          size_t     const & ak = 0 )
->decltype(aA.diagonal())
{
    return aA.diagonal();
}

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DIAG_VEC_EIGEN_HPP_ */

