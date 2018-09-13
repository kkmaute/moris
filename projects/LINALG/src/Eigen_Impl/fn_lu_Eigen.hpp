/*
 * fn_lu_Eigen.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LU_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LU_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include <Eigen/Dense>
#include "fn_trans_Eigen.hpp"
namespace moris
{

    void
    lu(       DDRMat & aL,
              DDRMat & aU,
              DDRMat & aP,
        const DDRMat & aA )
    {
        Eigen::PartialPivLU< DDRMat > tMyLU( aA );

        aU = tMyLU.matrixLU().triangularView< Eigen::Upper >();
        aL = tMyLU.matrixLU().triangularView< Eigen::UnitLower >();

        //aL = matrix_t::Identity( tMyLU.rows(), tMyLU.cols() );

       //aL.triangularView< Eigen::StrictlyLower >() = tMyLU.matrixLU();

        aP = tMyLU.permutationP();
    }
}



#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LU_EIGEN_HPP_ */
