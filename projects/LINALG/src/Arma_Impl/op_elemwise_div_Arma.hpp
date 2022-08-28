/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_elemwise_div_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_OP_ELEMWISE_DIV_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_OP_ELEMWISE_DIV_ARMA_HPP_

namespace moris
{

//--------------------------------------------------------------------------------

    template< typename Type >
    auto
    operator/( const arma::Mat<Type> & aA, const arma::Mat<Type> & aB )
        ->decltype( aA.data() / aB.data() )
    {
        return  aA.data() / aB.data();
    }

//--------------------------------------------------------------------------------

    template< typename Matrix_Type, typename ET >
    auto
    operator/( const Matrix< Matrix_Type > & aA, const ET & aB )
        ->decltype( aA.matrix_data() / aB )
    {
        return  aA.matrix_data() / aB;
    }

//--------------------------------------------------------------------------------

    template< typename Matrix_Type, typename ET >
    auto
    operator/( const ET & aA,
               const Matrix< Matrix_Type > & aB )
        ->decltype( aA / aB.matrix_data() )
    {
        return aA / aB.matrix_data();
    }

//--------------------------------------------------------------------------------
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_OP_ELEMWISE_DIV_ARMA_HPP_ */

