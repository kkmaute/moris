/*
 * fn_svd_Eigen.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SVD_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SVD_EIGEN_HPP_

#include <Eigen/Dense>

namespace moris {

    template< typename T, typename S >
    void
    svd(      Eigen::MatrixBase< T >  & aU,
              Eigen::MatrixBase< S >  & aS,
              Eigen::MatrixBase< T >  & aV,
        const Eigen::MatrixBase< T >  & aM)
    {
        Eigen::JacobiSVD< T > tMySVD(
                aM,
                  Eigen::DecompositionOptions::ComputeThinU
                | Eigen::DecompositionOptions::ComputeThinV );

        aS = tMySVD.singularValues();
        aU = tMySVD.matrixU();
        aV = tMySVD.matrixV();

    }

    template< typename T, typename S >
    void
    svd(      Eigen::MatrixBase< S >      & aS,
            const Eigen::MatrixBase< T >  & aM)
    {
        Eigen::JacobiSVD< T > tMySVD(
                aM,
                Eigen::DecompositionOptions::ComputeThinU
                | Eigen::DecompositionOptions::ComputeThinV );

        aS = tMySVD.singularValues();
    }

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_SVD_EIGEN_HPP_ */
