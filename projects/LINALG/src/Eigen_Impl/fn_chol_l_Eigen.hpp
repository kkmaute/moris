/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_chol_l_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CHOL_L_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CHOL_L_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
    template<typename ET>
    Eigen::MatrixXd
    chol_l( const Eigen::MatrixBase<ET> & aA )
    {
        Eigen::LLT<ET> lltOfA( aA );
        return lltOfA.matrixL();
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CHOL_L_EIGEN_HPP_ */

