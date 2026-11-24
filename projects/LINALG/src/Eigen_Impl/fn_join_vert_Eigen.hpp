/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include "cl_Matrix.hpp"
#include "fn_Eigen_Utilities.hpp"

namespace moris
{
    template< class A, class B >
    auto join_vert( const Eigen::MatrixBase< A >& aA, const Eigen::MatrixBase< B >& aB )
    {
        MORIS_ERROR( aA.cols() == aB.cols() || ( aA.rows() == 0 && aA.cols() == 0 )
                             || ( aB.rows() == 0 && aB.cols() == 0 ),
                "join_vert: column counts differ." );

        linalg_internal::eigen_plain< A > tResult( aA.rows() + aB.rows(), std::max( aA.cols(), aB.cols() ) );
        if ( aA.size() != 0 )
        {
            tResult.topRows( aA.rows() ) = aA;
        }
        if ( aB.size() != 0 )
        {
            tResult.bottomRows( aB.rows() ) = aB;
        }
        return tResult;
    }
}    // namespace moris
