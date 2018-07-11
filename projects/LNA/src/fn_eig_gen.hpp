#ifndef MORIS_LINALG_FN_EIG_GEN_HPP_
#define MORIS_LINALG_FN_EIG_GEN_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

#ifdef MORIS_USE_ARMA

namespace arma_Math
{
    template< typename T >
    void
    eig_gen(
            moris::Mat< moris::cplx > & eigenvalues,
            moris::Mat< moris::cplx > & eigenvectors,
            const moris::Mat< T > &     aA )
    {
        arma::Col< moris::cplx > tTmpS; ///< arma column type required for arma::eig_gen.

        bool tSuccess = arma::eig_gen( tTmpS, eigenvectors.data(), aA.data() );
        if ( !tSuccess )
        {
             throw std::runtime_error ("eig_gen failed. ");
        }

        eigenvalues = tTmpS;
    }
}

#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    void
    eig_gen(
            moris::Mat< moris::cplx > & eigenvalues,
            moris::Mat< moris::cplx > & eigenvectors,
            const moris::Mat< T > &     aA )
    {
        Eigen::EigenSolver< Eigen::Matrix <T, Eigen::Dynamic, Eigen::Dynamic > > eigensolver( aA.data() );

        eigenvalues = eigensolver.eigenvalues();
        eigenvectors = eigensolver.eigenvectors();
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{

    /**
     * @brief Eigenvalue decomposition of dense general Matrix A.
     *
     * @param[in] aA Matrix.
     *
     * @param[out] eigenvalues The eigenvalues of Matrix aA.
     *
     * @param[out] eigenvectors The eigenvectors of Matrix aA.
     */
    template< typename T >
    void
    eig_gen(
            moris::Mat< moris::cplx > & eigenvalues,
            moris::Mat< moris::cplx > & eigenvectors,
            const moris::Mat< T > &     aA)
    {
        moris::Math::eig_gen( eigenvalues, eigenvectors, aA );
    }
}

#endif /* MORIS_LINALG_FN_EIG_GEN_HPP_ */
