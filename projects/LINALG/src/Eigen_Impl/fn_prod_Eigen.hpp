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
    auto prod( const Eigen::MatrixBase< E >& a )
    {
        return linalg_internal::eigen_reduction< E, true >( a.derived(), 0 );
    }
}    // namespace moris
