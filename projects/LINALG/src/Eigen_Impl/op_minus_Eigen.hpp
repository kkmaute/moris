/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_minus_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_MINUS_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_MINUS_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
    template< typename Matrix_Type, typename ET >
    auto
    operator-( const Eigen::MatrixBase<ET> &  aA,
               Matrix< Matrix_Type > & aB )
    ->decltype( aA - aB.matrix_data() )
    {
        return  aA - aB.matrix_data();
    }

    template< typename Matrix_Type, typename ET >
    auto
    operator-( Matrix< Matrix_Type > & aA,
               const Eigen::MatrixBase<ET> &  aB)
    ->decltype( aA.matrix_data() - aB )
    {
        return  aA.matrix_data() - aB;
    }

    template< typename Matrix_Type >
    auto
    operator-( const Matrix< Matrix_Type > & aMatrix,
               const typename Matrix< Matrix_Type >::Data_Type &  aScalar )
    ->decltype( aMatrix.matrix_data() - aScalar*Eigen::MatrixXd::Ones( aMatrix.n_rows(), aMatrix.n_cols() ) )
    {
        return  aMatrix.matrix_data() - aScalar*Eigen::MatrixXd::Ones( aMatrix.n_rows(), aMatrix.n_cols() );
    }

    template< typename Matrix_Type >
    auto
    operator-( const typename Matrix< Matrix_Type >::Data_Type &  aScalar,
               const Matrix< Matrix_Type > & aMatrix )
    ->decltype( aScalar*Eigen::MatrixXd::Ones( aMatrix.n_rows(), aMatrix.n_cols() ) - aMatrix.matrix_data() )
    {
        return  aScalar*Eigen::MatrixXd::Ones( aMatrix.n_rows(), aMatrix.n_cols() ) - aMatrix.matrix_data();
    }

}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_MINUS_EIGEN_HPP_ */

