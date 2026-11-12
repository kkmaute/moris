/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"
#include "fn_Eigen_Utilities.hpp"
#include <cstring>

namespace moris
{
    template< class E, class M, class I = uint >
    void sort( const Eigen::MatrixBase< E >& a, Matrix< M >& result, const char* direction, I dim = 0 )
    {
        MORIS_ERROR( ( std::strcmp( direction, "ascend" ) == 0 || std::strcmp( direction, "descend" ) == 0 ) && dim < 2, "sort: invalid direction or dimension." );
        linalg_internal::eigen_plain< E > copy       = a;
        const bool                        descending = std::strcmp( direction, "descend" ) == 0;
        for ( Eigen::Index outer = 0; outer < ( dim == 0 ? a.cols() : a.rows() ); ++outer )
        {
            std::vector< typename E::Scalar > values( dim == 0 ? a.rows() : a.cols() );
            for ( size_t k = 0; k < values.size(); ++k )
            {
                values[ k ] = dim == 0 ? copy( k, outer ) : copy( outer, k );
            }
            std::sort( values.begin(), values.end(), [ & ]( const auto& x, const auto& y ) { return descending ? linalg_internal::eigen_less( y, x ) : linalg_internal::eigen_less( x, y ); } );
            for ( size_t k = 0; k < values.size(); ++k )
            {
                if ( dim == 0 )
                {
                    copy( k, outer ) = values[ k ];
                }
                else
                {
                    copy( outer, k ) = values[ k ];
                }
            }
        }
        result = copy;
    }

    template< class E, class M >
    void sort( const Eigen::MatrixBase< E >& a, Matrix< M >& result )
    {
        sort( a, result, "ascend", 0 );
    }
}    // namespace moris
