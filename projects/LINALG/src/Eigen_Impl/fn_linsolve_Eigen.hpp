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
        using Real  = typename Eigen::NumTraits< typename A::Scalar >::Real;

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

        if ( tA.rows() == tA.cols() && tA.allFinite() && tB.allFinite() )
        {
            Eigen::PartialPivLU< Plain > tLU( tA );
            const Real                   tEpsilon = std::numeric_limits< Real >::epsilon();

            // PartialPivLU assumes invertibility. Reject zero pivots before
            // estimating conditioning, and retain COD's rank policy for nearly
            // singular systems even when LU could produce a small residual.
            if ( tLU.matrixLU().allFinite() && ( tLU.matrixLU().diagonal().array().abs() > Real( 0 ) ).all() )
            {
                const Real tRcond = tLU.rcond();
                if ( std::isfinite( tRcond ) && tRcond > std::sqrt( tEpsilon ) )
                {
                    tResult = tLU.solve( tB );
                    if ( tResult.allFinite() )
                    {
                        Plain tResidual( tA.rows(), tB.cols() );
                        tResidual.noalias() = tA * tResult;
                        tResidual -= tB;

                        // The 1-norm avoids squaring extreme values and shares
                        // the norm used by LU's reciprocal-condition estimate.
                        const Real tANorm     = tA.cwiseAbs().colwise().sum().maxCoeff();
                        const Real tTolerance = Real( 32 ) * tEpsilon * Real( tA.rows() );
                        bool       tAccept    = true;
                        for ( Eigen::Index j = 0; j < tB.cols(); ++j )
                        {
                            // Check each RHS independently so a large column
                            // cannot hide a poor solution in another column.
                            const Real tScale = tANorm * tResult.col( j ).cwiseAbs().sum() + tB.col( j ).cwiseAbs().sum();
                            const Real tError = tResidual.col( j ).cwiseAbs().sum();
                            if ( !std::isfinite( tScale ) || !std::isfinite( tError )
                                    || ( tScale == Real( 0 ) ? tError != Real( 0 ) : tError / tScale > tTolerance ) )
                            {
                                tAccept = false;
                                break;
                            }
                        }
                        if ( tAccept )
                        {
                            return tResult;
                        }
                    }
                }
            }
        }

        // Preserve the existing minimum-norm least-squares behavior for
        // rectangular, rank-deficient, ill-conditioned, or rejected LU solves.
        tResult = tA.completeOrthogonalDecomposition().solve( tB );
        return tResult;
    }
}    // namespace moris
