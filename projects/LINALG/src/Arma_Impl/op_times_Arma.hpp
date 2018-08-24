/*
 * op_times_Arma.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_OP_TIMES_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_OP_TIMES_ARMA_HPP_


#include "cl_Matrix.hpp"
#include <armadillo>


namespace moris
{
template< typename T1, typename T2, typename ET >
auto
operator*( const ET &  aA,
           Mat_New< T1, T2 > & aB )
->decltype( aA * aB.matrix_data() )
{
    return  aA * aB.matrix_data();
}

template< typename T1, typename T2, typename ET >
auto
operator*( Mat_New< T1, T2 > & aA,
           const ET &  aB)
->decltype( aA.matrix_data() * aB )
{
    return  aA.matrix_data() * aB;
}

}




#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_OP_TIMES_ARMA_HPP_ */
