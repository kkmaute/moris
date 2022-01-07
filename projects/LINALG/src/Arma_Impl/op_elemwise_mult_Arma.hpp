/*
 * op_elemwise_mult.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_OP_ELEMWISE_MULT_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_OP_ELEMWISE_MULT_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
//--------------------------------------------------------------------------------

    template< typename Matrix_Type >
    auto
    operator%(
            const  Matrix< Matrix_Type > & aA,
            const  Matrix< Matrix_Type > & aB )
    -> decltype ( aA.matrix_data() % aB.matrix_data() )
    {
        return aA.matrix_data() % aB.matrix_data();
    }

//--------------------------------------------------------------------------------

    template< typename Type >
    auto
    operator%( const arma::Mat<Type> & aA, const arma::Mat<Type> & aB )
        ->decltype( aA.data() % aB.data() )
    {
        return  aA.data() % aB.data();
    }

//--------------------------------------------------------------------------------

    template<typename Matrix_Type, typename ET >
    auto
    operator%( const Matrix< Matrix_Type > & aA, const ET & aB )
        ->decltype( aA.matrix_data() % aB )
    {
        return  aA.matrix_data() % aB;
    }

//--------------------------------------------------------------------------------

    template<typename Matrix_Type, typename ET >
    auto
    operator%( const ET & aA, const Matrix< Matrix_Type > & aB )
        ->decltype( aA % aB.matrix_data() )
    {
        return aA % aB.matrix_data();
    }

//--------------------------------------------------------------------------------
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_OP_ELEMWISE_MULT_ARMA_HPP_ */
