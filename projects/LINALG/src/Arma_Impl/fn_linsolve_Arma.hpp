/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_linsolve_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSOLVE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSOLVE_ARMA_HPP_

#include "cl_Matrix.hpp"

//#include "viennacl/linalg/cg.hpp"
//#include "viennacl/linalg/bicgstab.hpp"
//#include "viennacl/linalg/gmres.hpp"
//#include "viennacl/io/matrix_market.hpp"

namespace moris
{
    template< typename ET >
    auto
    solve( const ET          & aA,
           const ET          & aB,
           const std::string & aSolver = "default" )
    -> decltype( arma::solve( aA, aB ) )
    {
#ifdef MORIS_HAVE_DEBUG
        return arma::solve( aA, aB );
#else
        return arma::solve( aA, aB, arma::solve_opts::fast );
#endif
    }

//    template< typename ET >
//    auto
//    solve( const sparse_Mat          & aA,
//           const ET          & aB,
//           const std::string & aSolver = "default" )
//    -> decltype( arma::spsolve( A.data(), B.data() ) )
//    {
//        if ( aSolver == "superlu" )
//        {
//            // set options for superlu
//            arma::superlu_opts superlu_settings;
//            superlu_settings.permutation  = arma::superlu_opts::COLAMD;
//            superlu_settings.refine       = arma::superlu_opts::REF_NONE;
//            superlu_settings.equilibrate  = true;
//            superlu_settings.symmetric    = false;
//            superlu_settings.pivot_thresh = 1.0;
//
//            return arma::spsolve( aA, aB, aSolver.c_str(), superlu_settings );
//        }
//        else if ( aSolver == "cg" )
//        {
//            typedef double ScalarType;
//            const arma::Col<ScalarType>   arma_rhs    = arma::conv_to< arma::Col<ScalarType>   >::from( aB );
//            const arma::Col<ScalarType>   arma_result = viennacl::linalg::solve( aA, arma_rhs, viennacl::linalg::cg_tag(1e-8,500) );
//            const arma::Col<ScalarType>   residual    =  aA * arma_result - arma_rhs;
//            std::fprintf(stdout,"cg: relative residual: %e\n", viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(arma_rhs) );
//            return arma::conv_to<arma::Mat< T >>::from(arma_result);
//        }
//        else if ( aSolver == "bicgstab" )
//        {
//            typedef double ScalarType;
//            const arma::Col<ScalarType>   arma_rhs    = arma::conv_to< arma::Col<ScalarType>   >::from( aB );
//            const arma::Col<ScalarType>   arma_result = viennacl::linalg::solve( aA, arma_rhs, viennacl::linalg::bicgstab_tag(1e-8,500,250));
//            const arma::Col<ScalarType>   residual    =  aA * arma_result - arma_rhs;
//            std::fprintf(stdout,"bicgstab: relative residual: %e\n", viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(arma_rhs) );
//            return arma::conv_to<arma::Mat< T >>::from(arma_result);
//        }
//        else if ( aSolver == "gmres" )
//        {
//            typedef double ScalarType;
//            const arma::Col<ScalarType>   arma_rhs    = arma::conv_to< arma::Col<ScalarType>   >::from( aB );
//            const arma::Col<ScalarType>   arma_result = viennacl::linalg::solve( aA, arma_rhs, viennacl::linalg::gmres_tag(1e-8,500,50));
//            const arma::Col<ScalarType>   residual    =  aA * arma_result - arma_rhs;
//            std::fprintf(stdout,"gmres: relative residual: %e\n", viennacl::linalg::norm_2(residual) / viennacl::linalg::norm_2(arma_rhs) );
//            return arma::conv_to<arma::Mat< T >>::from(arma_result);
//        }
//        else
//        {
//            MORIS_LOG_WARNING << "Armadillo doesn't support this solver, switching to SuperLU solver.";
//            return arma::spsolve( aA, aB, "superlu" );
//        }
//    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_LINSOLVE_ARMA_HPP_ */

