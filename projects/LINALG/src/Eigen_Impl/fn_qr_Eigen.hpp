/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_qr_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_QR_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_QR_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include <Eigen/Dense>
#include "fn_trans_Eigen.hpp"
namespace moris
{
    template< typename ET >
    void
    qr(       Eigen::MatrixBase<ET> & aQ,
              Eigen::MatrixBase<ET> & aR,
        const Eigen::MatrixBase<ET> & aA )
    {
        Eigen::HouseholderQR< ET > tQR( aA );
        aQ = tQR.householderQ();
        aR = trans( aQ )*aA;
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_QR_EIGEN_HPP_ */

