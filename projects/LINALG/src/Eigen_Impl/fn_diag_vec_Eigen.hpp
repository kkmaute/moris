/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"

namespace moris
{
    template< class E >
    auto diag_vec( const Eigen::MatrixBase< E >& a, size_t k = 0 ) -> decltype( a.diagonal( static_cast< Eigen::Index >( k ) ) )
    {
        return a.diagonal( static_cast< Eigen::Index >( k ) );
    }
}    // namespace moris
