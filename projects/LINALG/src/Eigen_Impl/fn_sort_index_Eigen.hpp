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
    template< class E >
    Matrix< DDUMat > sort_index( const Eigen::MatrixBase< E >& a, const char* direction = "ascend" )
    {
        MORIS_ERROR( std::strcmp( direction, "ascend" ) == 0 || std::strcmp( direction, "descend" ) == 0, "sort_index: invalid direction." );
        auto             indices = linalg_internal::eigen_sorted_indices( a.derived(), std::strcmp( direction, "descend" ) == 0 );
        Matrix< DDUMat > result( indices.size(), 1 );
        for ( size_t i = 0; i < indices.size(); ++i )
        {
            result( i ) = indices[ i ];
        }
        return result;
    }
}    // namespace moris
