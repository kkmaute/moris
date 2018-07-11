#ifndef MORIS_LINALG_FN_DET_HPP_
#define MORIS_LINALG_FN_DET_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    auto
    det(
            moris::Mat< T > const & aA )
    -> decltype( arma::det( aA.data() ) )
    {
        return arma::det( aA.data() );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    auto
    det(
            moris::Mat< T > const & aA )
    -> decltype( aA.data().determinant() )
    {
        return aA.data().determinant();
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Calculate the determinant of square matrix.
     *
     *@param[in] aA A given square matrix
     *
     * Example:
     * @include LNA/src/fn_det.inc
     *
     */
    template< typename T, bool M = ( std::is_same< T, moris::real >::value), typename std::enable_if< M >::type* = nullptr >
    auto
    det(
            moris::Mat< T > const & aA )
    -> decltype( moris::Math::det( aA ) )
    {
        return moris::Math::det( aA );
    }
}

#endif  /* MORIS_LINALG_FN_DET_HPP_ */
