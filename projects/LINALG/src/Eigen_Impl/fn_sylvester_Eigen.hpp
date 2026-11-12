/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"
#include "fn_Eigen_Utilities.hpp"

namespace moris
{
    template< class A, class B, class C >
    linalg_internal::eigen_plain< A > sylvester( const Eigen::MatrixBase< A >& aA, const Eigen::MatrixBase< B >& aB, const Eigen::MatrixBase< C >& aC )
    {
        using Real    = typename Eigen::NumTraits< typename A::Scalar >::Real;
        using Complex = std::complex< Real >;
        using Plain   = Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic >;

        MORIS_ERROR( aA.rows() == aA.cols() && aB.rows() == aB.cols() && aC.rows() == aA.rows() && aC.cols() == aB.rows(), "sylvester: incompatible dimensions." );

        Eigen::ComplexSchur< Plain > tSchurA( aA.template cast< Complex >() ), tSchurB( aB.template cast< Complex >() );
        if ( tSchurA.info() != Eigen::Success || tSchurB.info() != Eigen::Success )
        {
            throw std::runtime_error( "sylvester: Schur decomposition failed" );
        }

        Plain tRhs = -tSchurA.matrixU().adjoint() * aC.template cast< Complex >() * tSchurB.matrixU();
        Plain tY   = Plain::Zero( aC.rows(), aC.cols() );
        // Test singularity relative to the operator, independent of its units.
        const Real tScale = tSchurA.matrixT().norm() + tSchurB.matrixT().norm();
        for ( Eigen::Index j = 0; j < aC.cols(); ++j )
        {
            for ( Eigen::Index i = aC.rows(); i-- > 0; )
            {
                Complex tValue = tRhs( i, j );
                for ( Eigen::Index k = i + 1; k < aC.rows(); ++k )
                {
                    tValue -= tSchurA.matrixT()( i, k ) * tY( k, j );
                }
                for ( Eigen::Index k = 0; k < j; ++k )
                {
                    tValue -= tY( i, k ) * tSchurB.matrixT()( k, j );
                }
                const Complex tDivisor = tSchurA.matrixT()( i, i ) + tSchurB.matrixT()( j, j );
                MORIS_ERROR( std::abs( tDivisor ) > std::numeric_limits< Real >::epsilon() * tScale, "sylvester: singular equation." );
                tY( i, j ) = tValue / tDivisor;
            }
        }

        Plain tResult = tSchurA.matrixU() * tY * tSchurB.matrixU().adjoint();
        if constexpr ( Eigen::NumTraits< typename A::Scalar >::IsComplex )
        {
            return tResult;
        }
        else
        {
            return tResult.real();
        }
    }

    template< class A, class B, class C, class M >
    void sylvester( const Eigen::MatrixBase< A >& aA, const Eigen::MatrixBase< B >& aB, const Eigen::MatrixBase< C >& aC, Matrix< M >& aX )
    {
        aX = sylvester( aA, aB, aC );
    }
}    // namespace moris
