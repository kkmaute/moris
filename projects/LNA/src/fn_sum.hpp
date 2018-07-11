/*
 * fn_sum.hpp
 *
 *  Created on: Apr 18, 2017
 *      Author: gleim
 */

#ifndef SRC_LINALG_FN_SUM_HPP_
#define SRC_LINALG_FN_SUM_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    auto
    sum(
            moris::Mat< T > const & aA )
    -> decltype( arma::accu( aA.data() ) )
    {
        return arma::accu( aA.data() );
    }

    template< typename T >
    auto
    sum(
            moris::Sp_Mat< T > const & aA )
    -> decltype( arma::accu( aA.data() ) )
    {
        return arma::accu( aA.data() );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    auto
    sum(
            moris::Mat< T > const & aA )
    -> decltype( aA.data().sum() )
    {
        return aA.data().sum();
    }

    template< typename T >
    auto
    sum(
            moris::Sp_Mat< T > const & aA )
    -> decltype( aA.data().sum() )
    {
        return aA.data().sum();
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Calculate the sum of a matrix.
     *
     *@param[in] aA A given matrix
     *
     * Example:
     * @include LNA/src/fn_sum.inc
     *
     */
    template< typename T>
    auto
    sum(
            moris::Mat< T > const & aA )
    -> decltype( moris::Math::sum( aA ) )
    {
        return moris::Math::sum( aA );
    }

    template< typename T>
    auto
    sum(
            moris::Sp_Mat< T > const & aA )
    -> decltype( moris::Math::sum( aA ) )
    {
        return moris::Math::sum( aA );
    }
}


#endif /* SRC_LINALG_FN_SUM_HPP_ */
