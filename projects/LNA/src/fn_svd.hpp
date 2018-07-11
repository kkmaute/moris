#ifndef MORIS_LINALG_FN_SVD_HPP_
#define MORIS_LINALG_FN_SVD_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T1, typename T2 >
    void
    svd(
            moris::Mat< T2 >       & aS,
            moris::Mat< T1 > const & aA )
    {
        arma::Col< T2 > tTmpS; ///< arma column type required for arma::svd.

        bool tSuccess = arma::svd( tTmpS, aA.data() );

        if ( !tSuccess )
        {
            throw std::runtime_error ("SVD failed. ");
        }
        aS = tTmpS;

        return;
    }

    template< typename T1, typename T2 >
    void
    svd(
            moris::Mat< T1 >       & aU,
            moris::Mat< T2 >       & aS,
            moris::Mat< T1 >       & aV,
            moris::Mat< T1 > const & aA )
    {
        arma::Col< T2 > tTmpS; ///< arma column type required for arma::svd.

        bool tSuccess = arma::svd( aU.data(), tTmpS, aV.data(), aA.data() );

        if ( !tSuccess )
        {
            throw std::runtime_error ("SVD failed. ");
        }

        aS = tTmpS;
        return;
    }
}
#endif
// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T1, typename T2 >
    void
    svd(
            moris::Mat< T2 >       &  aS,
            moris::Mat< T1 > const & aA )
    {
        Eigen::JacobiSVD< Eigen::Matrix< T1, Eigen::Dynamic, Eigen::Dynamic > > tMySDV( aA.data() );
        aS = tMySDV.singularValues();

        return;
    }

    template< typename T1, typename T2 >
    void
    svd(
            moris::Mat< T1 >       & aU,
            moris::Mat< T2 >       & aS,
            moris::Mat< T1 >       & aV,
            moris::Mat< T1 > const & aA )
    {
        Eigen::JacobiSVD< Eigen::Matrix< T1, Eigen::Dynamic, Eigen::Dynamic > > tMySVD(
                aA.data(), Eigen::DecompositionOptions::ComputeThinU | Eigen::DecompositionOptions::ComputeThinV );

        aS = tMySVD.singularValues();
        aU = tMySVD.matrixU();
        aV = tMySVD.matrixV();

        return;
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Singular Value decomposition of a dense matrix. This method
     * requests just the singular value decomposition vector (aS), which are
     * arranged in descending order.
     *
     * @param[in] aA An input Matrix.
     *
     * @param[out] aS A diagonal matrix of the same dimensions as aA, with
     *                nonnegative diagonal elements in decreasing order.
     *
     * Example:
     * @include LNA/src/fn_svd/fn_svd_short.inc
    */
    template< typename T1, typename T2 >
    void
    svd(
            moris::Mat< T1 >       & aS,
            moris::Mat< T2 > const & aA )
    {
        moris::Math::svd( aS, aA);
        return;
    }

    /**
     * @brief Singular Value decomposition of dense matrix aA. This method
     * requests the singular value decomposition vector (aS), along with the
     * left and right unitary matrices aU and aV
     *
     * @param[in] aA Matrix.
     *
     * @param[out] aU, aS, aV A diagonal matrix (aS) of the same dimensions as
     * aA, with nonnegative diagonal elements in decreasing order and unitary
     * matrices (aU,aV) such that @f$\mathbf{A}_{ij}
     *  = \mathbf{U}_{ik}  \mathbf{s}_{ks} \mathbf{V}_{js} @f$
     *
     *  Example:
     *  @include LNA/src/fn_svd/fn_svd_long.inc
    */
    template< typename T1, typename T2 >
    void
    svd(
            moris::Mat< T2 >       & aU,
            moris::Mat< T1 >       & aS,
            moris::Mat< T2 >       & aV,
            moris::Mat< T2 > const & aA )
    {
        moris::Math::svd( aU, aS, aV, aA );
        return;
    }
}

#endif  /* MORIS_LINALG_FN_SVD_HPP_ */
