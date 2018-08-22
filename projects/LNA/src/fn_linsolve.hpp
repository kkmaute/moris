#ifndef MORIS_LINALG_FN_LINSOLVE_HPP_
#define MORIS_LINALG_FN_LINSOLVE_HPP_

// Third-party header files.
// Currently conflict with Armadillo using superlu, thus superlu support for Eigen is disabled
#ifdef MORIS_USE_EIGEN_XXX
#include "Eigen/SuperLUSupport"
#endif
#ifdef MORIS_USE_EIGEN
#include <Eigen/UmfPackSupport>
#endif

#ifdef MORIS_USE_ARMA
#define VIENNACL_WITH_ARMADILLO 1
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/bicgstab.hpp"
#include "viennacl/linalg/gmres.hpp"
#include "viennacl/io/matrix_market.hpp"
#endif

// MORIS library header files.
#include "cl_Mat.hpp"
#include "cl_Sp_Mat.hpp"
#include "cl_Logger.hpp" // IOS/src

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    auto
    solve(
            moris::Mat< T > const & A,
            moris::Mat< T > const & B,
            std::string     const & aSolver = "default" )
    -> decltype( arma::solve( A.data(), B.data() ) )
    {
        return arma::solve( A.data(), B.data() );
    }

    template< typename T >
    auto
    solve(
            moris::Sp_Mat< T > const & A,
            moris::Mat< T >    const & B,
            std::string        const & aSolver = "superlu" )
    -> decltype( arma::spsolve( A.data(), B.data() ) )
    {
        if ( aSolver == "superlu" )
        {
            // set options for superlu
            arma::superlu_opts superlu_settings;
            superlu_settings.permutation  = arma::superlu_opts::COLAMD;
            superlu_settings.refine       = arma::superlu_opts::REF_NONE;
            superlu_settings.equilibrate  = true;
            superlu_settings.symmetric    = false;
            superlu_settings.pivot_thresh = 1.0;

            return arma::spsolve( A.data(), B.data(), aSolver.c_str(), superlu_settings );
        }
        else if ( aSolver == "cg" )
        {
            typedef double ScalarType;
            const arma::Col<ScalarType>   arma_rhs    = arma::conv_to< arma::Col<ScalarType>   >::from( B.data() );
            const arma::Col<ScalarType>   arma_result = viennacl::linalg::solve( A.data(), arma_rhs,
                    viennacl::linalg::cg_tag(1e-8,500));
            const arma::Col<ScalarType>   residual    =  A.data() * arma_result - arma_rhs;
            std::fprintf(stdout,"cg: relative residual: %e\n", viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(arma_rhs) );
            return arma::conv_to<arma::Mat< T >>::from(arma_result);
        }
        else if ( aSolver == "bicgstab" )
        {
            typedef double ScalarType;
            const arma::Col<ScalarType>   arma_rhs    = arma::conv_to< arma::Col<ScalarType>   >::from( B.data() );
            const arma::Col<ScalarType>   arma_result = viennacl::linalg::solve( A.data(), arma_rhs,
                    viennacl::linalg::bicgstab_tag(1e-8,500,250));
            const arma::Col<ScalarType>   residual    =  A.data() * arma_result - arma_rhs;
            std::fprintf(stdout,"bicgstab: relative residual: %e\n", viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(arma_rhs) );
            return arma::conv_to<arma::Mat< T >>::from(arma_result);
        }
        else if ( aSolver == "gmres" )
        {
            typedef double ScalarType;
            const arma::Col<ScalarType>   arma_rhs    = arma::conv_to< arma::Col<ScalarType>   >::from( B.data() );
            const arma::Col<ScalarType>   arma_result = viennacl::linalg::solve( A.data(), arma_rhs,
                    viennacl::linalg::gmres_tag(1e-8,500,50));
            const arma::Col<ScalarType>   residual    =  A.data() * arma_result - arma_rhs;
            std::fprintf(stdout,"gmres: relative residual: %e\n", viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(arma_rhs) );
            return arma::conv_to<arma::Mat< T >>::from(arma_result);
        }
        else
        {
            MORIS_LOG_WARNING << "Armadillo doesn't support this solver, switching to SuperLU solver.";
            return arma::spsolve( A.data(), B.data(), "superlu" );
        }
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    auto
    solve(
            moris::Mat< T > const & A,
            moris::Mat< T > const & B,
            std::string     const & aSolver = "default" )
    -> decltype( Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >( A.data().colPivHouseholderQr().solve( B.data() ) ) )
    {
        Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > tempX;
        tempX = A.data().colPivHouseholderQr().solve( B.data() ); // uses rank-revealing QR decomposition of a matrix with column-pivoting
        return tempX;
    }

    template< typename T >
    auto
    solve(
            moris::Sp_Mat< T > const & A,
            moris::Mat< T >    const & B,
            std::string        const & aSolver = "superlu" )
    -> decltype( Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >() )
    {
//        Eigen::ComputationInfo tSuccess = Eigen::Success;
        Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > X;

        if ( aSolver == "superlu_XXX" )
        {
#ifdef MORIS_USE_EIGEN_XXX
            // Set the solver type
            Eigen::SuperLU<Eigen::SparseMatrix< T > > solver;

            // Compute the decomposition of sparse matrix
            solver.compute( A.data() ) ;

            if( solver.info() != tSuccess ) {
                MORIS_LOG_ERROR << "Decomposition failed in Eigen sparse matrix solver";}

            // Solve the linear system
            X = solver.solve( B.data() );

//            if( solver.info() != tSuccess ) {
//                MORIS_LOG_ERROR << "Solving failed in Eigen sparse matrix solver";}
#endif
        }
        else if ( aSolver == "sparselu" )
        {
            typedef Eigen::SparseMatrix< T, 0, long int > tSpMat;

            // Set the solver type
            Eigen::SparseLU< tSpMat, Eigen::COLAMDOrdering< long int > >  solver;

            // Compute the ordering permutation vector from the structural pattern of A
            solver.analyzePattern( A.data() );

            // Compute the numerical factorization
            solver.factorize( A.data() );

            // Use the factors to solve the linear system
            X = solver.solve( B.data() );
        }
        else if ( aSolver == "umfpack" )
        {
//            // Set the solver type
//            Eigen::UmfPackLU<Eigen::SparseMatrix< T > > solver;
//
//            // Compute the decomposition of sparse matrix
//            solver.compute( A.data() ) ;
//
//            if( solver.info() != tSuccess ) {
//                MORIS_LOG_ERROR << "Decomposition failed in Eigen sparse matrix solver";}
//
//            // Solve the linear system
//            X = solver.solve( B.data() );
//
//            if( solver.info() != tSuccess ) {
//                MORIS_LOG_ERROR << "Solving failed in Eigen sparse matrix solver";}
        }
        else
        {
//            MORIS_LOG_WARNING << "Eigen doesn't support this solver, switching to UmfPackLU solver.";

            // Set the solver type
//            Eigen::UmfPackLU<Eigen::SparseMatrix< T > > solver;
//
//            // Compute the decomposition of sparse matrix
//            solver.compute( A.data() ) ;
//
//            if( solver.info() != tSuccess ) {
//                MORIS_LOG_ERROR << "Decomposition failed in Eigen sparse matrix solver";}
//
//            // Solve the linear system
//            X = solver.solve( B.data() );
//
//            if( solver.info() != tSuccess ) {
//                MORIS_LOG_ERROR << "Solving failed in Eigen sparse matrix solver";}
        }

        return X;
    }
}

#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Solve for a linear set of equations Ax = B.
     *
     * @param[in] A The LHS Matrix
     * @param[in] B The RHS Vector
     * @param[in] aSolver Optional solver type
     *
     * @return The vector of solutions, x. Similar to B/A in Matlab
     *
     * @note If A is square, solve() is faster and more accurate than using X = inv(A)*B .
     * If A is non-square, solve() will try to provide approximate solutions to under-determined
     * as well as over-determined systems.
     * Eigen provides various options for decomposition of matrices to facilitate a linear solve.
     * We use the default option of QR decomposition with column pivoting.
     *
     * Example:
     * @include LNA/src/fn_linsolve/linsolve_real.inc
     * @include LNA/src/fn_linsolve/linsolve_complex.inc
     */
    template< typename T >
    auto
    solve(
            moris::Mat< T > const & A,
            moris::Mat< T > const & B,
            std::string     const & aSolver = "default" )
    -> decltype( moris::Math::solve( A, B ) )
    {
        return moris::Math::solve( A, B );
    }

    template< typename T >
    auto
    solve(
            moris::Sp_Mat< T > const & A,
            moris::Mat< T >    const & B,
            std::string        const & aSolver = "superlu" )
    -> decltype( moris::Math::solve( A, B ) )
    {
        return moris::Math::solve( A, B, aSolver.c_str() );
    }
}

#endif  /* MORIS_LINALG_FN_LINSOLVE_HPP_ */
