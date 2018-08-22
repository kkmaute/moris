/*
 * fn_norm.hpp
 *
 *  Created on: Aug 22, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LNA_SRC_FN_NORM_HPP_
#define PROJECTS_LNA_SRC_FN_NORM_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    auto
    norm(
            moris::Mat< T > const & aA )
    -> decltype( arma::norm( aA.data(), 2 ) )
    {
        return arma::norm( aA.data(), 2 );
    }

    template< typename T >
    auto
    norm( T  const & aA )
    -> decltype( arma::norm( aA, 2 ) )
    {
        return arma::norm( aA, 2 );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    auto
    norm( moris::Mat< T > const & aA )
    -> decltype( aA.data().norm() )
    {
        return aA.data().norm() ;
    }

    template< typename T >
    auto
    norm(const Eigen::MatrixBase<T> &  aA )
    -> decltype(  aA.norm() )
    {
        return aA.norm();
    }
}

namespace moris
{
template< typename T >
auto
norm( T const & aA )
-> decltype( moris::Math::norm( aA ) )
{
    return moris::Math::norm( aA );
}
}

#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Calculate the norm of a matrix.
     *
     *@param[in] aA A given matrix
     *
     *
     */
    template< typename T >
    auto
    norm( moris::Mat< T > const & aA )
    -> decltype( moris::Math::norm( aA ) )
    {
        return moris::Math::norm( aA );
    }
}



#endif /* PROJECTS_LNA_SRC_FN_NORM_HPP_ */
