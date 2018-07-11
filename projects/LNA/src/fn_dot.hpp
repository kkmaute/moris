#ifndef MORIS_LINALG_FN_DOT_HPP_
#define MORIS_LINALG_FN_DOT_HPP_

// C++ header files.
#include <utility>

// MORIS library header files.
#include "cl_Mat.hpp"


// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    auto
    dot(
            moris::Mat< T > const & aA,
            moris::Mat< T > const & aB)
    -> decltype( arma::accu( aA.data() % aB.data() ) )
    {
        return arma::accu( aA.data() % aB.data() );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    auto
    dot(
            moris::Mat< T > const & aA,
            moris::Mat< T > const & aB)
        -> decltype( aA.data().cwiseProduct(aB.data()).sum() )
        {
    		return aA.data().cwiseProduct(aB.data()).sum();
        }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Calculates the dot (scalar) product between two matrices or vectors of the same size
     * @f$ s=\mathbf{A}_{i,j}*\mathbf{B}_{i,j}@f$ \n
     *
     *@param[in] aA A given matrix
     *@param[in] aB A given matrix
     *
     * Example:
     * @include LNA/src/fn_dot.inc
     *
     */
    template< typename T, bool M = ( std::is_same< T, moris::real >::value), typename std::enable_if< M >::type* = nullptr >
    auto
    dot(
            moris::Mat< T > const & aA,
            moris::Mat< T > const & aB)
    -> decltype( moris::Math::dot( aA, aB ) )
    {
        return moris::Math::dot( aA, aB );
    }
}

#endif  /* MORIS_LINALG_FN_DOT_HPP_ */
