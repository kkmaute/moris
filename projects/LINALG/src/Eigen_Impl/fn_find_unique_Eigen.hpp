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
    template< class E >
    auto find_unique( const Eigen::MatrixBase< E >& a )
    {
        auto                indices = linalg_internal::eigen_sorted_indices( a.derived() );
        std::vector< uint > chosen;
        for ( auto i : indices )
        {
            if ( chosen.empty() || a.derived().coeff( i ) != a.derived().coeff( chosen.back() ) )
            {
                chosen.push_back( i );
            }
        }
        // Armadillo returns the selected original indices in increasing index order.
        std::sort( chosen.begin(), chosen.end() );
        Eigen::Matrix< uint, Eigen::Dynamic, Eigen::Dynamic > result( chosen.size(), 1 );
        for ( size_t i = 0; i < chosen.size(); ++i )
        {
            result( i ) = chosen[ i ];
        }
        return result;
    }
}    // namespace moris
