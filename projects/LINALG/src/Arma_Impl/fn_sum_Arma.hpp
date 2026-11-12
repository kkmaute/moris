/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_sum_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{

    // MORIS sum without a dimension reduces every entry to a scalar, unlike
    // arma::sum's default column reduction. accu also accepts lazy expressions.
    template< typename Type >
    Type
    sum( const arma::Mat< Type >& aA )
    {
        return arma::accu( aA );
    }

    template< typename Type, typename Expression >
    Type
    sum( const arma::Base< Type, Expression >& aA )
    {
        return arma::accu( aA.get_ref() );
    }

    template< typename Type >
    auto
    sum( const arma::Mat< Type >& aA, const uint aDim )
            -> decltype( arma::sum( aA, aDim ) )
    {
        return arma::sum( aA, aDim );
    }

    template< typename Type, typename Expression >
    auto
    sum( const arma::Base< Type, Expression >& aA, const uint aDim )
            -> decltype( arma::sum( aA.get_ref(), aDim ) )
    {
        return arma::sum( aA.get_ref(), aDim );
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_ */
