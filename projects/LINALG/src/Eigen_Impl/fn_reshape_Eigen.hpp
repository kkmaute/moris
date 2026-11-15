/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include "fn_Eigen_Utilities.hpp"
#include <type_traits>

namespace moris::linalg_internal
{
    template< typename E >
    struct eigen_reshape_map : std::false_type
    {
    };

    template< typename Plain, int Options, typename Stride >
    struct eigen_reshape_map< Eigen::Map< Plain, Options, Stride > > : std::true_type
    {
    };

    template< typename E >
    struct eigen_reshape : eigen_reduction_operand< E >
    {
        using eigen_reduction_operand< E >::values;
        using Nested = std::decay_t< decltype( values ) >;
        using Scalar = typename E::Scalar;

        Eigen::Index  rows;
        bool          rowOrder;
        Eigen::Index  sourceRows;
        Eigen::Index  sourceCols;
        Eigen::Index  sourceSize;
        const Scalar* data = nullptr;

        eigen_reshape( const E& aValues, Eigen::Index aRows, bool aRowOrder )
                : eigen_reduction_operand< E >( aValues )
                , rows( aRows )
                , rowOrder( aRowOrder )
                , sourceRows( values.rows() )
                , sourceCols( values.cols() )
                , sourceSize( values.size() )
        {
            using Nesting = Eigen::internal::nested_eval< E, 1 >;
            // A plain matrix or transpose reference can reallocate or change
            // dimensions before deferred evaluation. Only owned snapshots and
            // value-owned Maps provide a stable address and shape to cache.
            if constexpr ( ( Nesting::Evaluate
                                   || ( eigen_reshape_map< Nested >::value && !std::is_reference< typename Nesting::type >::value ) )
                           && ( Eigen::internal::traits< Nested >::Flags & Eigen::DirectAccessBit )
                           && has_eigen_data< Nested >::value )
            {
                if ( values.innerStride() == 1
                        && ( values.outerSize() == 1 || values.outerStride() == values.innerSize() ) )
                {
                    // Evaluated operands share ownership in the base. Copying
                    // this functor therefore keeps its cached data pointer valid.
                    data = values.data();
                }
            }
        }

        Scalar operator()( Eigen::Index i, Eigen::Index j ) const
        {
            const auto tIndex      = i + j * rows;
            const auto tSourceSize = data != nullptr ? sourceSize : values.size();
            if ( tIndex >= tSourceSize )
            {
                return Scalar( 0 );
            }
            if ( data != nullptr && rowOrder == static_cast< bool >( Nested::IsRowMajor ) )
            {
                return data[ tIndex ];
            }

            const auto tSourceRows = data != nullptr ? sourceRows : values.rows();
            const auto tSourceCols = data != nullptr ? sourceCols : values.cols();
            const auto tRow        = rowOrder ? tIndex / tSourceCols : tIndex % tSourceRows;
            const auto tCol        = rowOrder ? tIndex % tSourceCols : tIndex / tSourceRows;
            if ( data != nullptr )
            {
                return data[ Nested::IsRowMajor ? tRow * sourceCols + tCol : tCol * sourceRows + tRow ];
            }
            return values.coeff( tRow, tCol );
        }
    };
}    // namespace moris::linalg_internal

namespace moris
{
    template< typename E >
    auto reshape( const Eigen::MatrixBase< E >& aA, size_t aRows, size_t aCols )
    {
        // Retain the operand expression; padding/truncation follows Armadillo.
        return linalg_internal::eigen_plain< E >::NullaryExpr( aRows, aCols, linalg_internal::eigen_reshape< E >{ aA.derived(), static_cast< Eigen::Index >( aRows ), false } );
    }
}    // namespace moris
