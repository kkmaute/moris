#ifndef MORIS_LINALG_FN_INV_HPP_
#define MORIS_LINALG_FN_INV_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{

    template< typename T>
    auto
    inv(
            const moris::Mat< T > & aA) -> decltype( arma::inv( aA.data() ) )
            {
        return arma::inv( aA.data() );
            }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{

    template< typename T >
    auto
    inv(
            const moris::Mat< T > & aA ) -> decltype( (aA.data()).inverse() )
            {
        return (aA.data()).inverse();
            }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{

    /**
     * @brief Compute the inverse of Matrix A.
     *
     * @param[in] aA Matrix.
     *
     * @return The inverse of the Matrix aA.
     */
    template< typename T>
    auto
    inv(
            const moris::Mat< T > & aA ) -> decltype( moris::Math::inv( aA ) )
            {
        return moris::Math::inv( aA );
            }

}

#endif /* MORIS_FN_INV_HPP_ */
