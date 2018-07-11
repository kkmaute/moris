#ifndef MORIS_LINALG_FN_QR_HPP_
#define MORIS_LINALG_FN_QR_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"
#include "fn_trans.hpp"
#include "cl_Tuple.hpp" // CON/src

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    void
    qr(
            moris::Mat< T >       & aQ,
            moris::Mat< T >       & aR,
            moris::Mat< T > const & aA )
    {
        bool tSuccess = arma::qr( aQ.data(), aR.data(), aA.data() );

        if ( !tSuccess )
        {
            throw std::runtime_error( "QR decomposition failed!" );
        }
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    void
    qr(
            moris::Mat< T >       & aQ,
            moris::Mat< T >       & aR,
            moris::Mat< T > const & aA )
    {
        Eigen::HouseholderQR<
            Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >
                tMyQR( aA.data() );

        aQ.data() = tMyQR.householderQ();
        aR        = moris::trans(aQ)*aA.data();
    }
}
#endif

namespace moris
{
    /**
     * @brief QR decomposition of dense matrix
     *
     * @param[out] aQ Q-matrix
     * @param[out] aR R-matrix
     * @param[in]  aA Input matrix
     *
     * The QR decomposition of a matrix aA is such that \f$ A = Q*R \f$ where
     * \f$ Q^T*Q = I \f$ and R is an upper triangular non-square matrix. If A
     * is m-by-n, then Q is m-by-m and R is m-by-n.
     *
     * Example:
     * @include LNA/src/fn_qr.inc
     *
     */
    template< typename T >
    void
    qr(
            moris::Mat< T >       & aQ,
            moris::Mat< T >       & aR,
            moris::Mat< T > const & aA )
    {
        moris::Math::qr( aQ, aR, aA );
        return;
    }

    /**
     * @brief QR decomposition of dense matrix
     *
     * @param[in] aA Input matrix
     *
     * Performs the QR decomposition of the given matrix and returns the Q and R
     * matrices, in this order, as a Tuple.
     *
     * @usesTuple
     *
     * Example:
     * @include LNA/src/fn_qr.inc
     *
     */
    template< typename T >
    moris::Tuple< moris::Mat<T>, moris::Mat<T> >
    qr(
            moris::Mat< T > const & aA )
    {
        moris::Mat< T> tQ;
        moris::Mat< T> tR;
        moris::qr(tQ, tR, aA);

        return moris::make_tuple(tQ, tR);
    }

}

#endif /* MORIS_FN_QR_HPP_ */

