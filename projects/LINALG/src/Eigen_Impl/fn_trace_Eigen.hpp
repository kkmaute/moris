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
    auto trace( const Eigen::MatrixBase< E >& a )
    {
        return a.trace();
    }
}    // namespace moris
