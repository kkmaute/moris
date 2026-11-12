/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#pragma once

#include <Eigen/Dense>
#include "cl_Matrix.hpp"

namespace moris::linalg_internal
{
    template< class M >
    void eigen_eye( size_t rows, size_t cols, Matrix< M >& a )
    {
        a.set_size( rows, cols, 0 );
        a.matrix_data().setIdentity();
    }
    inline auto eye( size_t rows, size_t cols )
    {
        return Eigen::MatrixXd::Identity( rows, cols );
    }
}    // namespace moris::linalg_internal
