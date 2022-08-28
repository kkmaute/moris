/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_elemwise_div_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_DIV_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_DIV_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{

//--------------------------------------------------------------------------------

    template< typename Matrix_Type, typename ET >
    auto
    operator/( const Matrix< Matrix_Type > & aA,
               const Eigen::MatrixBase<ET> & aB )
        ->decltype( aA.matrix_data().cwiseQuotient( aB ) )
    {
        return aA.matrix_data().cwiseQuotient( aB );
    }

//--------------------------------------------------------------------------------

    template< typename Matrix_Type, typename ET >
    auto
    operator/( const Eigen::MatrixBase<ET> & aA,
               const Matrix< Matrix_Type > & aB )
        ->decltype( aA.cwiseQuotient( aB.matrix_data() ) )
    {
        return aA.cwiseQuotient( aB.matrix_data() );
    }

//--------------------------------------------------------------------------------

    template< typename ET1, typename ET2 >
    auto
    operator/( const Eigen::MatrixBase<ET1> & aA,
               const Eigen::MatrixBase<ET2> & aB )
        ->decltype( aA.cwiseQuotient( aB ) )
    {
        return aA.cwiseQuotient( aB );
    }

//--------------------------------------------------------------------------------
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_DIV_EIGEN_HPP_ */

