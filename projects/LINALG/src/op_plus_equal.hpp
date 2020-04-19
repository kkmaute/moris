/*
 * op_plus_equal.hpp
 *
 *  Created on: Feb 21, 2020
 *      Author: wunsch
 */

#ifndef PROJECTS_LINALG_SRC_OP_PLUS_EQUAL_HPP_
#define PROJECTS_LINALG_SRC_OP_PLUS_EQUAL_HPP_

#ifdef MORIS_USE_EIGEN
//#include "Eigen_Impl/op_plus_equal_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
//#include "op_plus_equal_Arma.hpp"
#endif

namespace moris
{

/**
 * @brief Plus equal operator
 *
 *
 * aA += aB
 *
 */


template< typename Matrix_Type >
auto
operator+=(
        Matrix< Matrix_Type > & aA,
        Matrix< Matrix_Type > const & aB )
        ->decltype( aA.matrix_data() += aB.matrix_data() )
        {
    return aA.matrix_data() += aB.matrix_data();
        }

}
#endif /* PROJECTS_LINALG_SRC_OP_PLUS_EQUAL_HPP_ */
