#ifndef MORIS_LINALG_FN_CHOL_L_HPP_
#define MORIS_LINALG_FN_CHOL_L_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA

namespace arma_Math
{
    template< typename T>
    void
    choll(
            moris::Mat< T >       & aL,
            moris::Mat< T > const & aA,
            const std::string &     aStr = "lower" )
    {
        aL = arma::chol( aA.data(), aStr.c_str() );

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
    choll(
            moris::Mat< T >       & aL,
            moris::Mat< T > const & aA,
            const std::string &     aStr = "lower" )
    {
        Eigen::LLT<Eigen::MatrixXd> lltOfA( aA.data() );
        aL = lltOfA.matrixL();

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
     *@param[out] aL Lower triangle matrix
     *@param[in]  aA Input matrix
     *@param[in]  aStr String flag for lower triangle decomposition
     *
     * This Cholesky decomposition uses only the diagonal and the lower triangle
     * of A to produce a lower triangular L so that
     * @f$ \mathbf{A}_{ij} = \mathbf{L}_{ik} \mathbf{L}_{jk} @f$
     * If A is not positive definite, an error message is printed.
     *
     * Example:
     * @include LNA/src/fn_chol_l.inc
     *
     */
    template< typename T >
    void
    choll(
            moris::Mat< T >       & aL,
            moris::Mat< T > const & aA,
            const std::string &     aStr = "lower" )
    {
        moris::Math::choll( aL, aA );

        return;
    }
}



#endif  /* MORIS_LINALG_FN_CHOL_L_HPP_ */
