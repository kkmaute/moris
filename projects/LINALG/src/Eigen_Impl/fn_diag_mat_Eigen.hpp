/*
 * fn_diag_Eigen.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DIAG_MAT_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DIAG_MAT_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "fn_isvector.hpp"
#include "Eigen/Dense"


namespace moris
{

template< typename ET>
auto
diag_mat( Eigen::MatrixBase<ET> & aA,
          size_t     const & ak = 0 )
->decltype(aA.asDiagonal())
{
    return aA.asDiagonal();
}
}



#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_DIAG_MAT_EIGEN_HPP_ */
