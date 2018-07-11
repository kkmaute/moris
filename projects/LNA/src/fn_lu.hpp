#ifndef MORIS_LINALG_FN_LU_HPP_
#define MORIS_LINALG_FN_LU_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"
#include "core.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    inline
    void
    lu(
            moris::Mat< moris::real >       & aL,
            moris::Mat< moris::real >       & aU,
            moris::Mat< moris::real >       & aP,
            moris::Mat< moris::real > const & aA )
    {
        bool tSuccess = arma::lu( aL.data(), aU.data(), aP.data(), aA.data() );

        if ( ! tSuccess )
        {
            throw std::runtime_error( "LU decomposition with partial pivoting failed." );
        }

        return;
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef  MORIS_USE_EIGEN
namespace eigen_Math
{
    inline
    void
    lu(
            moris::Mat< moris::real >       & aL,
            moris::Mat< moris::real >       & aU,
            moris::Mat< moris::real >       & aP,
            moris::Mat< moris::real > const & aA )
    {
        ///< This operator works for double case. Data type with this operator for Eigen should be double.
        typedef Eigen::Matrix< moris::real, Eigen::Dynamic, Eigen::Dynamic > matrix_t;

        Eigen::PartialPivLU< matrix_t > tMyLU( aA.data() );

        matrix_t tU = tMyLU.matrixLU().triangularView< Eigen::Upper >();

        matrix_t tL = matrix_t::Identity( tMyLU.rows(), tMyLU.cols() );

        tL.triangularView< Eigen::StrictlyLower >() = tMyLU.matrixLU();

        matrix_t tP = tMyLU.permutationP();

        aL = tL;
        aU = tU;
        aP = tP;

        return;
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief LU decomposition with partial pivoting.
     *
     *@param[out] aL Lower-triangular output matrix
     *@param[out] aU Upper-triangular output matrix
     *@param[out] aP Permutation matrix
     *@param[in]  aA Input matrix
     *
     * The LU decomposition of a matrix aA is such that \f$ P*A = L*U \f$ where L is a lower
     * triangular matrix, U is an upper triangular matrix, and P is a permutation matrix.
     *
     * Example:
     * @include LNA/src/fn_lu.inc
     *
     */
    inline
    void
    lu(
            moris::Mat< moris::real >       & aL,
            moris::Mat< moris::real >       & aU,
            moris::Mat< moris::real >       & aP,
            moris::Mat< moris::real > const & aA )
    {
        moris::Math::lu( aL, aU, aP, aA );

        return;
    }
}

#endif /* MORIS_FN_LU_HPP_ */
