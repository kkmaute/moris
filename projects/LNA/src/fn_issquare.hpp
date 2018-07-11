#ifndef MORIS_LINALG_FN_ISSQUARE_HPP_
#define MORIS_LINALG_FN_ISSQUARE_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    bool
    issquare(
            moris::Mat< T > const & aA )
    {
        return aA.data().is_square();
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    bool
    issquare(
            moris::Mat< T > const & aA )
    {
        if ( aA.data().rows() == aA.data().cols() )
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
     * @brief Checks if a given object is square
     *
     * @param[in] aA Input object to be checked
     *
     * This function operates on given object (i.e. matrix, vector) and returns
     * true if the object can be interpreted as a square matrix (i.e. number of
     * columns equal to number of rows). It also return false if the object
     * isn't square matrix.
     *
     * Example:
     * @include LNA/src/fn_issquare.inc
     *
     */
    template< typename T >
    bool
    issquare(
            moris::Mat< T > const & aA )
    {
        return moris::Math::issquare( aA );
    }
}

#endif  /* MORIS_LINALG_FN_ISSQUARE_HPP_ */
