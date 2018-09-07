/*
 * op_plus_Eigen.hpp
 *
 *  Created on: Aug 27, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_PLUS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_PLUS_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
    template< typename T1, typename T2, typename ET >
    auto
    operator+( const Eigen::MatrixBase<ET> &  aA,
               Matrix< T1, T2 > & aB )
    ->decltype( aA + aB.matrix_data() )
    {
        return  aA + aB.matrix_data();
    }

    template< typename T1, typename T2, typename ET >
    auto
    operator+( Matrix< T1, T2 > & aA,
               const Eigen::MatrixBase<ET> &  aB)
    ->decltype( aA.matrix_data() + aB )
    {
        return  aA.matrix_data() + aB;
    }


    template< typename T1, typename T2 >
    auto
    operator+( const Matrix< T1, T2 > & aMatrix,
               const T1 &  aScalar )
    ->decltype( aMatrix.matrix_data() + aScalar*Eigen::MatrixXd::Ones( aMatrix.n_rows(), aMatrix.n_cols() ) )
    {
        return  aMatrix.matrix_data() + aScalar*Eigen::MatrixXd::Ones( aMatrix.n_rows(), aMatrix.n_cols() );
    }

    template< typename T1, typename T2 >
    auto
    operator+( const T1 &  aScalar,
               const Matrix< T1, T2 > & aMatrix )
    ->decltype( aScalar*Eigen::MatrixXd::Ones( aMatrix.n_rows(), aMatrix.n_cols() ) + aMatrix.matrix_data() )
    {
        return  aScalar*Eigen::MatrixXd::Ones( aMatrix.n_rows(), aMatrix.n_cols() ) + aMatrix.matrix_data();
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_PLUS_EIGEN_HPP_ */
