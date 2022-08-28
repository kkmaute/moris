/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_chol_u_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CHOL_U_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CHOL_U_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
    template<typename ET>
    Eigen::MatrixXd
    chol_u( const Eigen::MatrixBase<ET> & aA )
    {
        Eigen::LLT<ET> lltOfA( aA );
        return lltOfA.matrixU();
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_CHOL_U_EIGEN_HPP_ */

