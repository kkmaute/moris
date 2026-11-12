/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"
#include <algorithm>
#include <complex>
#include <numeric>
#include <vector>
#include <stdexcept>
#include <limits>
#include <memory>

namespace moris::linalg_internal
{
    template< class E >
    using eigen_plain = Eigen::Matrix< typename E::Scalar, Eigen::Dynamic, Eigen::Dynamic >;
    template< class T >
    bool eigen_less( const T& a, const T& b )
    {
        if constexpr ( Eigen::NumTraits< T >::IsComplex )
        {
            return std::abs( a ) < std::abs( b );
        }
        else
        {
            return a < b;
        }
    }

    template< class E >
    std::vector< Eigen::Index > eigen_sorted_indices( const E& a, bool descending = false )
    {
        std::vector< Eigen::Index > indices( a.size() );
        std::iota( indices.begin(), indices.end(), 0 );
        std::stable_sort( indices.begin(), indices.end(), [ & ]( auto i, auto j ) { return descending ? eigen_less( a.coeff( j ), a.coeff( i ) ) : eigen_less( a.coeff( i ), a.coeff( j ) ); } );
        return indices;
    }
    // Follow Eigen's nesting policy for unevaluated operands. Products that
    // require evaluation instead share an immutable snapshot: Eigen copies
    // functors as expressions compose, which must not copy a dense product.
    template< class E, bool Evaluate = Eigen::internal::nested_eval< E, 1 >::Evaluate >
    struct eigen_reduction_operand
    {
        typename Eigen::internal::nested_eval< E, 1 >::type values;

        explicit eigen_reduction_operand( const E& aValues )
                : values( aValues )
        {
        }
    };

    template< class E >
    struct eigen_reduction_operand< E, true >
    {
      private:
        std::shared_ptr< const typename E::PlainObject > mStorage;

      public:
        const typename E::PlainObject& values;

        explicit eigen_reduction_operand( const E& aValues )
                : mStorage( std::make_shared< const typename E::PlainObject >( aValues ) )
                , values( *mStorage )
        {
        }
    };

    template< class E, bool Product >
    struct eigen_reduce : eigen_reduction_operand< E >
    {
        using eigen_reduction_operand< E >::values;
        unsigned dim;

        eigen_reduce( const E& aValues, unsigned aDim )
                : eigen_reduction_operand< E >( aValues )
                , dim( aDim )
        {
        }

        typename E::Scalar operator()( Eigen::Index i, Eigen::Index j ) const
        {
            typename E::Scalar value = Product ? 1 : 0;
            const auto         n     = dim == 0 ? values.rows() : values.cols();
            for ( Eigen::Index k = 0; k < n; ++k )
            {
                if constexpr ( Product )
                {
                    value *= dim == 0 ? values.coeff( k, j ) : values.coeff( i, k );
                }
                else
                {
                    value += dim == 0 ? values.coeff( k, j ) : values.coeff( i, k );
                }
            }
            return value;
        }
    };

    template< class E, bool Product >
    auto eigen_reduction( const E& a, unsigned dim )
    {
        MORIS_ERROR( dim < 2, "Reduction dimension must be zero or one." );
        return eigen_plain< E >::NullaryExpr( dim == 0 ? 1 : a.rows(), dim == 0 ? a.cols() : 1, eigen_reduce< E, Product >{ a, dim } );
    }
}    // namespace moris::linalg_internal
