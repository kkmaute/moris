/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace moris::linalg_internal
{
    template< typename E, typename = void >
    struct has_eigen_data : std::false_type
    {
    };

    template< typename E >
    struct has_eigen_data< E, std::void_t< decltype( std::declval< const E& >().data() ), decltype( std::declval< const E& >().innerStride() ), decltype( std::declval< const E& >().outerStride() ) > > : std::true_type
    {
    };

    template< typename E, typename = void >
    struct has_eigen_operands : std::false_type
    {
    };

    template< typename E >
    struct has_eigen_operands< E, std::void_t< decltype( std::declval< const E& >().lhs() ), decltype( std::declval< const E& >().rhs() ) > > : std::true_type
    {
    };

    template< typename E, typename = void >
    struct has_eigen_nested : std::false_type
    {
    };

    template< typename E >
    struct has_eigen_nested< E, std::void_t< decltype( std::declval< const E& >().nestedExpression() ) > > : std::true_type
    {
    };

    template< typename E, typename = void >
    struct has_eigen_functor_values : std::false_type
    {
    };

    template< typename E >
    struct has_eigen_functor_values< E, std::void_t< decltype( std::declval< const E& >().functor().values ) > > : std::true_type
    {
    };

    // Inspect native expression operands without evaluating them. A conservative
    // bounding interval includes strided blocks as well as contiguous matrices.
    // The destination is writable storage whose coefficients may be uninitialized;
    // only its address is inspected. Keeping that contract avoids GCC treating
    // the pointer as an input buffer under -Wmaybe-uninitialized in release builds.
    template< typename E, typename T >
    bool eigen_references( const E& aExpression, T* aData, size_t aCount )
    {
        if constexpr ( ( Eigen::internal::traits< E >::Flags & Eigen::DirectAccessBit ) && has_eigen_data< E >::value )
        {
            if ( !aCount || !aExpression.size() )
            {
                return false;
            }
            const auto tInner        = aExpression.innerStride();
            const auto tOuter        = aExpression.outerStride();
            const auto tInnerEnd     = ( aExpression.innerSize() - 1 ) * tInner;
            const auto tOuterEnd     = ( aExpression.outerSize() - 1 ) * tOuter;
            const auto tFirst        = std::min< Eigen::Index >( 0, tInnerEnd ) + std::min< Eigen::Index >( 0, tOuterEnd );
            const auto tLast         = std::max< Eigen::Index >( 0, tInnerEnd ) + std::max< Eigen::Index >( 0, tOuterEnd );
            const auto tFirstAddress = reinterpret_cast< std::uintptr_t >( aExpression.data() + tFirst );
            const auto tLastAddress  = reinterpret_cast< std::uintptr_t >( aExpression.data() + tLast );
            const auto tBegin        = std::min( tFirstAddress, tLastAddress );
            const auto tEnd          = std::max( tFirstAddress, tLastAddress ) + sizeof( typename E::Scalar );
            const auto tOwnBegin     = reinterpret_cast< std::uintptr_t >( aData );
            return tBegin < tOwnBegin + aCount * sizeof( T ) && tOwnBegin < tEnd;
        }
        else if constexpr ( has_eigen_operands< E >::value )
        {
            return eigen_references( aExpression.lhs(), aData, aCount ) || eigen_references( aExpression.rhs(), aData, aCount );
        }
        else if constexpr ( has_eigen_nested< E >::value )
        {
            return eigen_references( aExpression.nestedExpression(), aData, aCount );
        }
        else if constexpr ( has_eigen_functor_values< E >::value )
        {
            return eigen_references( aExpression.functor().values, aData, aCount );
        }
        return false;
    }

    template< typename Op, typename A, typename B, typename T >
    bool eigen_assignment_hazard( const Eigen::CwiseBinaryOp< Op, A, B >&, T*, size_t, size_t );
    template< typename Op, typename A, typename T >
    bool eigen_assignment_hazard( const Eigen::CwiseUnaryOp< Op, A >&, T*, size_t, size_t );
    template< typename A, typename T >
    bool eigen_assignment_hazard( const Eigen::Transpose< A >&, T*, size_t, size_t );

    template< typename E, typename T >
    bool eigen_assignment_hazard( const E& aExpression, T* aData, size_t aRows, size_t aCols )
    {
        if constexpr ( ( Eigen::internal::traits< E >::Flags & Eigen::DirectAccessBit ) && has_eigen_data< E >::value )
        {
            return eigen_references( aExpression, aData, aRows * aCols )
                && ( static_cast< const void* >( aExpression.data() ) != static_cast< const void* >( aData )
                        || aExpression.rows() != static_cast< Eigen::Index >( aRows )
                        || aExpression.cols() != static_cast< Eigen::Index >( aCols )
                        || ( aRows > 1 && aCols > 1 && E::IsRowMajor )
                        || aExpression.innerStride() != 1
                        || aExpression.outerStride() != static_cast< Eigen::Index >( E::IsRowMajor ? aCols : aRows ) );
        }
        // Products, transposes, inverses and other reordered reads need a
        // temporary only when they actually reference the destination.
        return eigen_references( aExpression, aData, aRows * aCols );
    }

    template< typename A, typename T >
    bool eigen_assignment_hazard( const Eigen::Transpose< A >& aExpression,
            T*                                                 aData,
            size_t                                             aRows,
            size_t                                             aCols )
    {
        return eigen_references( aExpression, aData, aRows * aCols );
    }

    template< typename Op, typename A, typename B, typename T >
    bool eigen_assignment_hazard( const Eigen::CwiseBinaryOp< Op, A, B >& aExpression,
            T*                                                            aData,
            size_t                                                        aRows,
            size_t                                                        aCols )
    {
        return eigen_assignment_hazard( aExpression.lhs(), aData, aRows, aCols )
            || eigen_assignment_hazard( aExpression.rhs(), aData, aRows, aCols );
    }

    template< typename Op, typename A, typename T >
    bool eigen_assignment_hazard( const Eigen::CwiseUnaryOp< Op, A >& aExpression,
            T*                                                        aData,
            size_t                                                    aRows,
            size_t                                                    aCols )
    {
        return eigen_assignment_hazard( aExpression.nestedExpression(), aData, aRows, aCols );
    }
}    // namespace moris::linalg_internal
