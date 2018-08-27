/*
 * op_times.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_OP_TIMES_HPP_
#define PROJECTS_LINALG_OP_TIMES_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_times_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_times_Arma.hpp"
#endif


namespace moris
{
/**
 * @brief Matrix multiplication
 *
 * @param[in] aA Input matrix
 * @param[in] aB Input matrix
 *
 * This function returns a matrix such that @f$ \mathbf{C}_{ij} = \mathbf{A}_{ik} \mathbf{B}_{kj} @f$
 * where C is the product of A and B. If A is an m-by-p and B is a p-by-n matrix, then C is an m-by-n matrix.
 *
 * Example:
 * @include LNA/src/op_times.inc
 *
 */
template< typename T1, typename T2 >
auto
operator*( Mat_New< T1, T2 > & aA,
           Mat_New< T1, T2 > & aB )
->decltype( aA.matrix_data() * aB.matrix_data() )
{
    return  aA.matrix_data() * aB.matrix_data();
}

}



#endif /* PROJECTS_LINALG_OP_TIMES_HPP_ */
