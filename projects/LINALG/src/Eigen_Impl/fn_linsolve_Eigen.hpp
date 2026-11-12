/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include <cmath>
#include "cl_Matrix.hpp"
#include "fn_Eigen_Utilities.hpp"

namespace moris
{
    template< class A, class B >
    auto solve( const Eigen::MatrixBase< A >& a, const Eigen::MatrixBase< B >& b, const std::string& solver = "default" )
    {
        MORIS_ERROR( a.rows() == b.rows(), "solve: incompatible right-hand side." );
        using Plain = linalg_internal::eigen_plain< A >;

        // Borrow compatible native storage. Ref evaluates other expressions once
        // and keeps that storage alive through factorization and residual checks.
        const Eigen::Ref< const Plain >                             tA( a.derived() );
        const Eigen::Ref< const linalg_internal::eigen_plain< B > > tB( b.derived() );
        Plain                                                       tResult;
        if ( tA.rows() == 0 || tA.cols() == 0 )
        {
            tResult.setZero( tA.cols(), tB.cols() );
            return tResult;
        }

        // Retain minimum-norm least-squares behavior for square, rectangular,
        // and rank-deficient systems.
        tResult = tA.completeOrthogonalDecomposition().solve( tB );
        return tResult;
    }
}    // namespace moris
