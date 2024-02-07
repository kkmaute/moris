/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_cross_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CROSS_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CROSS_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET1 >
    auto
    cross(
            arma::Mat< ET1 > const &aA,
            arma::Mat< ET1 > const &aB )
            -> decltype( arma::cross( aA, aB ) )
    {
        return arma::cross( aA, aB );
    }

    template< typename ET1, typename ET2 >
    auto
    cross(
            Matrix< ET1 > const &aA,
            ET2 const           &aB )
            -> decltype( arma::cross( aA.matrix_data(), aB ) )
    {
        return arma::cross( aA.matrix_data(), aB );
    }

    template< typename ET1, typename ET2 >
    auto
    cross(
            ET1 const           &aA,
            Matrix< ET2 > const &aB )
            -> decltype( arma::cross( aA, aB.matrix_data() ) )
    {
        return arma::cross( aA, aB.matrix_data() );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_CROSS_ARMA_HPP_ */
