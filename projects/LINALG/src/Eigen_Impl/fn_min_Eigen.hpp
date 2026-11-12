/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "fn_Eigen_Utilities.hpp"

namespace moris::linalg_internal
{
    template< typename Scalar >
    struct eigen_minimum
    {
        Scalar operator()( const Scalar& aA, const Scalar& aB ) const
        {
            return eigen_less( aA, aB ) ? aA : aB;
        }
    };
}    // namespace moris::linalg_internal
namespace moris
{
    template< typename A, typename B >
    auto min( const Eigen::MatrixBase< A >& aA, const Eigen::MatrixBase< B >& aB )
    {
        if constexpr ( Eigen::NumTraits< typename A::Scalar >::IsComplex )
        {
            return aA.binaryExpr( aB, linalg_internal::eigen_minimum< typename A::Scalar >{} );
        }
        else
        {
            return aA.cwiseMin( aB );
        }
    }

    template< typename A >
    auto min( const Eigen::MatrixBase< A >& aA )
    {
        if constexpr ( A::IsVectorAtCompileTime )
        {
            if constexpr ( Eigen::NumTraits< typename A::Scalar >::IsComplex )
            {
                return aA.redux( linalg_internal::eigen_minimum< typename A::Scalar >{} );
            }
            else
            {
                return aA.minCoeff();
            }
        }
        else if constexpr ( Eigen::NumTraits< typename A::Scalar >::IsComplex )
        {
            return aA.colwise().redux( linalg_internal::eigen_minimum< typename A::Scalar >{} );
        }
        else
        {
            return aA.colwise().minCoeff();
        }
    }
}    // namespace moris
