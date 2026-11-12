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
    template< class Q, class R, class E >
    void qr( Matrix< Q >& q, Matrix< R >& r, const Eigen::MatrixBase< E >& a )
    {
        using Plain = linalg_internal::eigen_plain< E >;
        Eigen::HouseholderQR< Plain > solver( a );
        Plain                         fullQ = solver.householderQ() * Plain::Identity( a.rows(), a.rows() );
        // The QR factor already stores R; avoid another dense product.
        Plain fullR = solver.matrixQR().template triangularView< Eigen::Upper >();
        q           = fullQ;
        r           = fullR;
    }
}    // namespace moris
