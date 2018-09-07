/*
 * fn_qr_Eigen.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: messe
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
