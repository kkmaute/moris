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
    template< class A, class B, class M, class I, class J >
    void intersect( const Eigen::MatrixBase< A >& a, const Eigen::MatrixBase< B >& b, Matrix< M >& result, Matrix< I >& ai, Matrix< J >& bi )
    {
        auto                              first  = linalg_internal::eigen_sorted_indices( a.derived() );
        auto                              second = linalg_internal::eigen_sorted_indices( b.derived() );
        std::vector< typename A::Scalar > values;
        std::vector< uint >               ia, ib;
        size_t                            i = 0, j = 0;
        while ( i < first.size() && j < second.size() )
        {
            const auto x = a.derived().coeff( first[ i ] ), y = b.derived().coeff( second[ j ] );
            if ( x == y )
            {
                if ( values.empty() || values.back() != x )
                {
                    values.push_back( x );
                    ia.push_back( first[ i ] );
                    ib.push_back( second[ j ] );
                }
                ++i;
                ++j;
            }
            else if ( linalg_internal::eigen_less( x, y ) )
            {
                ++i;
            }
            else
            {
                ++j;
            }
        }
        const bool tRowVector = a.rows() == 1 && b.rows() == 1;
        result.set_size( tRowVector ? 1 : values.size(), tRowVector ? values.size() : 1 );
        ai.set_size( values.size(), 1 );
        bi.set_size( values.size(), 1 );
        for ( size_t k = 0; k < values.size(); ++k )
        {
            result( k ) = values[ k ];
            ai( k )     = ia[ k ];
            bi( k )     = ib[ k ];
        }
    }

    template< class A, class B >
    auto intersect( const Eigen::MatrixBase< A >& a, const Eigen::MatrixBase< B >& b )
    {
        using Plain = linalg_internal::eigen_plain< A >;
        Matrix< Plain >  result;
        Matrix< DDUMat > ai, bi;
        intersect( a, b, result, ai, bi );
        return Plain( result.matrix_data() );
    }
}    // namespace moris
