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
    auto sum( const Eigen::MatrixBase< E >& a ) -> decltype( a.sum() )
    {
        return a.sum();
    }

    template< class E >
    auto sum( const Eigen::MatrixBase< E >& a, uint dim )
    {
        return linalg_internal::eigen_reduction< E, false >( a.derived(), dim );
    }
}    // namespace moris
