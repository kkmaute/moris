#ifndef MORIS_LINALG_FN_CTRANS_HPP_
#define MORIS_LINALG_FN_CTRANS_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    auto
    ctrans(
            const moris::Mat< T > & A )
    -> decltype( arma::trans( A.data() ) )
    {
        return arma::trans( A.data() ); // trans is an armadillo functionality
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    auto
    ctrans(
            const moris::Mat< T > & A )
    -> decltype( A.data().adjoint() )
    {
        return A.data().adjoint(); // adjoint is an Eigen functionality
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /*
     * @brief Computes the conjugate transpose of Matrix A. This functionality
     * is more costly than the 'trans' function, but may be necessary for complex arithmetic
     *
     * @param[in] A Matrix.
     *
     * @return The conjugate transpose of the Matrix A.
     *
     * Examples:
     * @include LNA/src/fn_ctrans/ctrans_real.inc
     * @include LNA/src/fn_ctrans/ctrans_complex.inc
     */
    template< typename T >
    auto
    ctrans(
            const moris::Mat< T > & A )
    -> decltype( moris::Math::ctrans( A ) )
    {
        return moris::Math::ctrans( A );
    }
}

#endif  /* MORIS_LINALG_FN_CTRANS_HPP_ */
