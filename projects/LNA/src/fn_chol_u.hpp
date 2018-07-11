#ifndef MORIS_LINALG_FN_CHOL_U_HPP_
#define MORIS_LINALG_FN_CHOL_U_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T>
    void
    cholu(
            moris::Mat< T >       & aU,
            moris::Mat< T > const & aA,
            const std::string &     aStr = "upper" )
    {
        aU = arma::chol( aA.data(), aStr.c_str() );

        return;
    }
}
#endif
// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    void
    cholu(
            moris::Mat< T >       & aU,
            moris::Mat< T > const & aA,
            const std::string &     aStr = "upper" )
    {
        Eigen::LLT<Eigen::MatrixXd> lltOfA( aA.data() );
        aU = lltOfA.matrixU();

        return;
    }
}
#endif
// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Cholesky decomposition of a Hermitian positive-definite matrix
     *
     *@param[out] aU Upper triangle matrix
     *@param[in]  aA Input Matrix
     *@param[in]  aStr String flag for upper triangle decomposition
     *
     * This Cholesky decomposition uses only the diagonal and the upper triangle
     * of A to produce an upper triangular U so that
     * @f$ \mathbf{A}_{ij} = \mathbf{U}_{ki} \mathbf{U}_{kj} @f$
     * If A is not positive definite, an error message is printed.
     *
     * Example:
     * @include LNA/src/fn_chol_u.inc
     *
     */
    template< typename T >
    void
    cholu(
            moris::Mat< T >       & aU,
            moris::Mat< T > const & aA,
            const std::string &     aStr = "upper" )
    {
        moris::Math::cholu( aU, aA );

        return;
    }
}

#endif  /* MORIS_LINALG_FN_CHOL_U_HPP_ */
