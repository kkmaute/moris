/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_linsolve_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LINSOLVE_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LINSOLVE_EIGEN_HPP_
#include <Eigen/Dense>

#ifdef MORIS_USE_EIGEN_XXX
#include "Eigen/SuperLUSupport"
#endif
#ifdef MORIS_USE_EIGEN
#include <Eigen/UmfPackSupport>
#endif

namespace moris
{
    template< typename ET >
    auto
    solve( const Eigen::MatrixBase<ET> & aA,
           const Eigen::MatrixBase<ET> & aB,
           const std::string           & aSolver = "default" )
    -> decltype( Eigen::Matrix< real, Eigen::Dynamic, Eigen::Dynamic >( aA.colPivHouseholderQr().solve( aB ) ) )
    {
        Eigen::Matrix< real, Eigen::Dynamic, Eigen::Dynamic > tempX;
        tempX = aA.colPivHouseholderQr().solve( aB ); // uses rank-revealing QR decomposition of a matrix with column-pivoting
        return tempX;
    }

//    template< typename ET >
//    auto
//    solve( const sparse_Mat & aA,
//           const Eigen::MatrixBase<ET> & aB,
//           const std::string           & aSolver = "superlu" )
//    -> decltype( Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic >() )
//	{
//
//	}
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LINSOLVE_EIGEN_HPP_ */

