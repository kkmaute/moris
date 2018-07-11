#ifndef MORIS_LINALG_FN_ISCOL_HPP_
#define MORIS_LINALG_FN_ISCOL_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T>
    bool
    iscol(
            moris::Mat< T > const & aA )
    {
        return aA.data().is_colvec();
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T>
    bool
    iscol(
            moris::Mat< T > const & aA )
    {
        if ( aA.data().rows() >= 1 && aA.data().cols() == 1 )
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
     * @brief Checks if a given object is column vector
     *
     * @param[in] aA Input object to be checked
     *
     * This function operates on given object (i.e. matrix, vector) and returns
     * true if the object can be interpreted as a column vector. It also return
     * false if the object isn't column vector.
     *
     * Example:
     * @include LNA/src/fn_iscol.inc
     *
     */
    template< typename T>
    bool
    iscol(
            moris::Mat< T > const & aA )
    {
        return moris::Math::iscol( aA );
    }
}

#endif  /* MORIS_LINALG_FN_ISCOL_HPP_ */
