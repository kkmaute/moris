/*
 * op_minus.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_OP_MINUS_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_OP_MINUS_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename Matrix_Type, typename ET >
    auto
    operator-( const ET &  aA,
               const Matrix< Matrix_Type > & aB )
    ->decltype( aA - aB.matrix_data() )
    {
        return  aA - aB.matrix_data();
    }

    template< typename Matrix_Type, typename ET >
    auto
    operator-( const Matrix< Matrix_Type > & aA,
               const ET &  aB)
    ->decltype( aA.matrix_data() - aB )
    {
        return  aA.matrix_data() - aB;
    }

    template< typename ET, typename Scalar >
    auto
    scalar_subtraction( const ET & aA,
                        const Scalar &  aB)
    ->decltype( aA.matrix_data() - aB )
    {
        return  aA.matrix_data() - aB;
    }

    template< typename ET, typename Scalar >
    auto
    scalar_subtraction( const Scalar &  aA,
                        const ET & aB)
    ->decltype( aA - aB.matrix_data() )
    {
        return  aA - aB.matrix_data();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_OP_MINUS_ARMA_HPP_ */
