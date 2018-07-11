#ifndef MORIS_LINALG_FN_ISEMPTY_HPP_
#define MORIS_LINALG_FN_ISEMPTY_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    bool
    isempty(
            moris::Mat< T > const & aA )
    {
        return aA.data().is_empty( );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    bool
    isempty(
            moris::Mat< T > const & aA )
    {
        // If either of the dimensions is zero, then the matrix is empty
        // Note: eigen returns a matrix [0x1] or [1x0] if is empty.
        if ( aA.data().rows() == 0 || aA.data().cols() == 0 )
        {
            return true;
        }
        return false;
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Checks if a given object is empty
     *
     * @param[in] aA Input object to be checked
     *
     * This function operates on given object (i.e. matrix, vector) and returns
     * true if the object has no element. It also return false if the object has
     * one or more elements.
     *
     * Example:
     * @include LNA/src/fn_isempty.inc
     *
     */
    template< typename T >
    bool
    isempty(
            moris::Mat< T > const & aA )
    {
        return moris::Math::isempty( aA );
    }
}

#endif /* MORIS_LINALG_FN_ISEMPTY_HPP_ */
