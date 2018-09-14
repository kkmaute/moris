/*
 * op_elemwise_mult.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: messe
 */

#ifndef PROJECTS_LINALG_SRC_OP_ELEMWISE_DIV_HPP_
#define PROJECTS_LINALG_SRC_OP_ELEMWISE_DIV_HPP_

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_elemwise_div_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_elemwise_div_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief Element wise division operator.
 *
 * @param[in] A Elements of A are the dividend.
 * @param[in] B Elements of B are the divisor.
 *
 * @return Creates a matrix corresponding to element-wise
 * division of the two input matrices.
 *USE_EIGEN
 * Example:
 * @include LNA/src/op_elemwise_div.inc
 *
 */
template< typename Matrix_Data >
auto
operator/(
        moris::Matrix< Matrix_Data > const & A,
        moris::Matrix< Matrix_Data > const & B )
-> decltype( operator/( A.matrix_data(), B.matrix_data() ) )
{
    return operator/( A.matrix_data(), B.matrix_data() );
}
}

#endif /* PROJECTS_LINALG_SRC_OP_ELEMWISE_DIV_HPP_ */
