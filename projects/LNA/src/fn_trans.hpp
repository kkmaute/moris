#ifndef MORIS_LINALG_FN_TRANS_HPP_
#define MORIS_LINALG_FN_TRANS_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"
#include "cl_Sp_Mat.hpp"

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    auto
    trans(
            moris::Base_Mat< T > const & A )
    -> decltype( arma::strans( A.data() ) )
    {
        return arma::strans( A.data() );
    }

    template< typename T >
    auto
    trans(
            T const & A )
    -> decltype( arma::strans( A ) )
    {
        return arma::strans( A );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    auto
    trans(
            moris::Base_Mat< T > const & A )
    -> decltype( A.data().transpose() )
    {
        return A.data().transpose();
    }

    template< typename T >
    auto
    trans(
            T const & A )
    -> decltype( A.transpose() )
    {
        return A.transpose();
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Computes transpose of Matrix A.
     *
     * @param[in] A Matrix.
     *
     * @return The transpose of the Matrix A.
     *
     * @note In case of complex numbers, no conjugate is taken.
     *
     * Example:
     * @include LNA/src/fn_trans/trans_real.inc
     * @include LNA/src/fn_trans/trans_complex.inc
     */
    template< typename T >
    auto
    trans(
            const moris::Base_Mat< T > & A )
    -> decltype( moris::Math::trans( A ) )
    {
        return moris::Math::trans( A );
    }

    /**
     * @brief Computes transpose of Matrix A.
     *
     * @param[in] A Matrix.
     *
     * @return The transpose of the Matrix A.
     *
     * This version of moris::trans can used in combination with other template
     * functions.
     */
    template< typename T >
    auto
    trans(
            T const & A )
    -> decltype( moris::Math::trans( A ) )
    {
        return moris::Math::trans( A );
    }

}

#endif  /* MORIS_LINALG_FN_TRANS_HPP_ */
