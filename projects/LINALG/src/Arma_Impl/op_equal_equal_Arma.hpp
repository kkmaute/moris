/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_equal_equal_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_OP_EQUAL_EQUAL_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_OP_EQUAL_EQUAL_ARMA_HPP_

namespace moris
{
template< typename Matrix_Type,
          typename T1,
          typename ET >
auto
operator==( const arma::eOp<T1, ET> &  aA,
           Matrix< Matrix_Type > & aB )
->decltype( aA == aB.matrix_data() )
{
    return  aA == aB.matrix_data();
}

template< typename Matrix_Type,
          typename T1,
          typename ET >
auto
operator==( Matrix< Matrix_Type > & aA,
            const arma::eOp<T1, ET> &  aB)
->decltype( aA.matrix_data()  == aB )
{
    return  aA.matrix_data() == aB;
}

template< typename Matrix_Type, typename ET >
auto
operator==( typename Matrix< Matrix_Type >::Data_Type & aA,
            const ET &  aB)
->decltype( aA  == aB )
{
    return  aA == aB;
}

template< typename Matrix_Type, typename ET >
auto
operator==( const ET &  aA,
            typename Matrix< Matrix_Type >::Data_Type & aB)
->decltype( aA  == aB )
{
    return  aA == aB;
}

}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_OP_EQUAL_EQUAL_ARMA_HPP_ */

