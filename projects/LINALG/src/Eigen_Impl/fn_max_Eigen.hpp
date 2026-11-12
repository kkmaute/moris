/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "fn_Eigen_Utilities.hpp"

namespace moris::linalg_internal
{
    template< typename Scalar >
    struct eigen_maximum
    {
        Scalar operator()( const Scalar& aA, const Scalar& aB ) const
        {
            return eigen_less( aB, aA ) ? aA : aB;
        }
    };
}    // namespace moris::linalg_internal
namespace moris
{
    template< typename A, typename B >
    auto max( const Eigen::MatrixBase< A >& aA, const Eigen::MatrixBase< B >& aB )
    {
        if constexpr ( Eigen::NumTraits< typename A::Scalar >::IsComplex )
        {
            return aA.binaryExpr( aB, linalg_internal::eigen_maximum< typename A::Scalar >{} );
        }
        else
        {
            return aA.cwiseMax( aB );
        }
    }

    template< typename A >
    auto max( const Eigen::MatrixBase< A >& aA )
    {
        if constexpr ( A::IsVectorAtCompileTime )
        {
            if constexpr ( Eigen::NumTraits< typename A::Scalar >::IsComplex )
            {
                return aA.redux( linalg_internal::eigen_maximum< typename A::Scalar >{} );
            }
            else
            {
                return aA.maxCoeff();
            }
        }
        else if constexpr ( Eigen::NumTraits< typename A::Scalar >::IsComplex )
        {
            return aA.colwise().redux( linalg_internal::eigen_maximum< typename A::Scalar >{} );
        }
        else
        {
            return aA.colwise().maxCoeff();
        }
    }
}    // namespace moris
