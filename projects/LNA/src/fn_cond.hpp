#ifndef MORIS_LINALG_FN_COND_HPP_
#define MORIS_LINALG_FN_COND_HPP_

// MORIS library header files.
#include "cl_Mat.hpp" // LNA/src
#include "fn_svd.hpp" // LNA/src

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T>
    auto
    cond(moris::Mat< T > const & aA)
    -> decltype( arma::cond( aA.data() ) )
    {
        return arma::cond( aA.data() );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T>
    moris::real
    cond(moris::Mat< T > const & aA)
    {
        moris::Mat< T > aS;
        moris::svd( aS, aA );

        moris::size_t lowerEigvalInd = 0;
        moris::size_t upperEigvalInd = aS.size( 0 )-1;

        return aS( lowerEigvalInd, 0 ) / aS( upperEigvalInd, 0 );
    }
}
#endif
// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Calculates the two-norm condition number of a matrix
     *
     * @return Calculates the two-norm condition number of a matrix (ratio
     * of the largest singular value to the smallest one).\n
     * Large condition numbers suggest that matrix A is nearly singular.
     *
     * @note Neither Armadillo nor Eigen compute the p-norm condition number.\n
     * Moreover, a one-norm condition number estimate function is not provided
     * either. If there is a need for a different condition number than a two-
     * norm, or a condition number estimate, it needs to be implemented.
     */
    template< typename T>
    auto
    cond(moris::Mat< T > const & aA)
    -> decltype( moris::Math::cond( aA) )
    {
        return moris::Math::cond( aA );
    }
}

#endif /* MORIS_LINALG_FN_COND_HPP_ */
