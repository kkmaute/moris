/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_prod_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_PROD_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_PROD_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template<typename ET >
    auto
    prod( const ET &  aA)
    ->decltype( arma::prod( aA ))
    {
        return arma::prod( aA );
    }

    template<typename ET >
    auto
    prod( Matrix<ET> const &  aA)
    ->decltype( arma::prod( aA.matrix_data() ))
    {
        return arma::prod( aA.matrix_data() );
    }

    template<typename ET1, typename ET2 >
    auto
    prod(   ET1 const &  aA,
            ET2 const &  aB)
    ->decltype( arma::prod( aA, aB ))
    {
        return arma::prod( aA,  aB );
    }

    template<typename ET1, typename ET2 >
    auto
    prod( Matrix<ET1> const &  aA,
                 ET2  const &  aB)
    ->decltype( arma::prod( aA.matrix_data(), aB ))
    {
        return arma::prod( aA.matrix_data(),  aB );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_PROD_ARMA_HPP_ */

