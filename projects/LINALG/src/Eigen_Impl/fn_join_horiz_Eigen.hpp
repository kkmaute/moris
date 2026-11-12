/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include "cl_Matrix.hpp"
#include "fn_Eigen_Utilities.hpp"

namespace moris
{
    template< class A, class B >
    auto join_horiz( const Eigen::MatrixBase< A >& a, const Eigen::MatrixBase< B >& b )
    {
        // A default-constructed 0x0 matrix is an empty accumulator, as in Armadillo.
        MORIS_ERROR( a.rows() == b.rows() || ( a.rows() == 0 && a.cols() == 0 )
                             || ( b.rows() == 0 && b.cols() == 0 ),
                "join_horiz: row counts differ." );
        linalg_internal::eigen_plain< A > result( std::max( a.rows(), b.rows() ), a.cols() + b.cols() );
        if ( a.size() != 0 )
        {
            result.leftCols( a.cols() ) = a;
        }
        if ( b.size() != 0 )
        {
            result.rightCols( b.cols() ) = b;
        }
        return result;
    }
}    // namespace moris
