/*
 * fn_cross_Arma.hpp
 *
 *  Created on: Jan 16, 2019
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CROSS_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CROSS_ARMA_HPP_



#include "cl_Matrix.hpp"
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>

namespace moris
{
    template< typename ET1, typename ET2>
    auto
    cross(
            Matrix<ET1> const & aA,
            ET2         const & aB)
    -> decltype( arma::cross( aA.matrix_data(), aB ) )
    {
        return arma::cross( aA.matrix_data(), aB );
    }

    template< typename ET1, typename ET2>
    auto
    cross(
            ET1         const & aA,
            Matrix<ET2> const & aB)
    -> decltype( arma::cross( aA, aB.matrix_data() ) )
    {
        return arma::cross( aA, aB.matrix_data() );
    }

    template< typename ET1, typename ET2>
    auto
    cross(
            ET1 const & aA,
            ET2 const & aB)
    -> decltype( arma::cross( aA, aB ) )
    {
        return arma::cross( aA, aB );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CROSS_ARMA_HPP_ */
