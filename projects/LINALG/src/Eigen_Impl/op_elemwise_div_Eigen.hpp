/*
 * op_elemwise_div.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_DIV_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_DIV_EIGEN_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
//--------------------------------------------------------------------------------

    template< typename Type, typename Matrix_Type >
    auto
    operator/(
            const  Matrix< Type, Matrix_Type > & aA,
            const  Matrix< Type, Matrix_Type > & aB )
        -> decltype ( aA.matrix_data().cwiseQuotient( aB.matrix_data() ) )
    {
        return aA.matrix_data().cwiseQuotient( aB.matrix_data() );
    }

//--------------------------------------------------------------------------------

    template< typename Type, typename Matrix_Type, typename ET >
    auto
    operator/( const Matrix< Type, Matrix_Type > & aA, const Eigen::MatrixBase<ET> & aB )
        ->decltype( aA.matrix_data().cwiseQuotient( aB ) )
    {
        return aA.matrix_data().cwiseQuotient( aB );
    }

//--------------------------------------------------------------------------------

    template< typename Type, typename Matrix_Type, typename ET >
    auto
    operator/( const Eigen::MatrixBase<ET> & aA, const Matrix< Type, Matrix_Type > & aB )
        ->decltype( aA.cwiseQuotient( aB.matrix_data() ) )
    {
        return aA.cwiseQuotient( aB.matrix_data() );
    }

//--------------------------------------------------------------------------------

    template< typename ET >
    auto
    operator/( const Eigen::MatrixBase<ET> & aA, const Eigen::MatrixBase<ET> & aB )
        ->decltype( aA.cwiseQuotient( aB ) )
    {
        return aA.cwiseQuotient( aB );
    }

//--------------------------------------------------------------------------------
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_OP_ELEMWISE_DIV_EIGEN_HPP_ */
