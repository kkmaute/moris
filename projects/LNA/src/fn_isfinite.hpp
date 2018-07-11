#ifndef MORIS_LINALG_FN_ISFINITE_HPP_
#define MORIS_LINALG_FN_ISFINITE_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    bool
    isfinite(
        moris::Mat< T > const & aA)
    {
        return aA.data().is_finite();
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    bool
    isfinite(
        moris::Mat< T > const & aA)
    {
        return aA.data().allFinite ();
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Checks if an object is finite
     *
     * @param[in] aA object to be checked.
     *
     * @return Checks if an object only contains finite entries or not,
     *         @f$ \mathbf{A}_{ii} = finite@f$\n
     *         Returns true if all elements of the object are finite.\n
     *         Returns false if at least one of the elements of the object
     *         is non-finite (±infinity or NaN).
     */
    template< typename T >
    bool
    isfinite(
        moris::Mat< T > const & aA)
    {
        return moris::Math::isfinite( aA );
    }
}

#endif /* MORIS_FN_ISFINITE_HPP_ */
