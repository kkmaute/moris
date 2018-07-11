#ifndef MORIS_LINALG_FN_EIG_SYM_HPP_
#define MORIS_LINALG_FN_EIG_SYM_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

#ifdef MORIS_USE_ARMA
namespace arma_Math
{

    template< typename T>
    void
    eig_sym(
            moris::Mat< T > &       eigenvalues,
            moris::Mat< T > &       eigenvectors,
            const moris::Mat< T > & aA )
    {
        arma::Col< T > tTmpS; ///< arma column type required for arma::eig_sym.

        bool tSuccess = arma::eig_sym( tTmpS, eigenvectors.data(), aA.data() );
        if ( !tSuccess )
        {
             throw std::runtime_error ("eig_sym failed. ");
        }

        eigenvalues = tTmpS;
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{

    template< typename T>
    void
    eig_sym(
            moris::Mat< T > &       eigenvalues,
            moris::Mat< T > &       eigenvectors,
            const moris::Mat< T > & aA )
    {
        Eigen::SelfAdjointEigenSolver< Eigen::Matrix < T, Eigen::Dynamic, Eigen::Dynamic > > eigensolver( aA.data() );

        eigenvalues = eigensolver.eigenvalues();
        eigenvectors = eigensolver.eigenvectors();
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{

    /**
     * @brief Eigenvalue decomposition of dense symmetric/Hermitian Matrix A.
     *
     * @param[in] aA Matrix.
     *
     * @param[out] eigenvalues The eigenvalues of Matrix aA.
     *
     * @param[out] eigenvectors The eigenvectors of Matrix aA.
     */

    template< typename T>
    void
    eig_sym(
            moris::Mat< T > &       eigenvalues,
            moris::Mat< T > &       eigenvectors,
            const moris::Mat< T > & aA)
    {
        return moris::Math::eig_sym( eigenvalues, eigenvectors, aA );
    }
}

#endif /* MORIS_LINALG_FN_EIG_SYM_HPP_ */
