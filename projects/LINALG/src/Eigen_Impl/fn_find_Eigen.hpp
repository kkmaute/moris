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
    auto find( const Eigen::DenseBase< E >& aMatrix, uint aCount = 0 )
    {
        Eigen::Matrix< uint, Eigen::Dynamic, Eigen::Dynamic > tResult( aMatrix.size(), 1 );
        Eigen::Index                                          tCount = 0;
        for ( Eigen::Index j = 0; j < aMatrix.cols(); ++j )
        {
            for ( Eigen::Index i = 0; i < aMatrix.rows(); ++i )
            {
                if ( aMatrix( i, j ) != typename E::Scalar( 0 ) && ( !aCount || tCount < aCount ) )
                {
                    tResult( tCount++ ) = i + j * aMatrix.rows();
                }
            }
        }
        tResult.conservativeResize( tCount, 1 );
        return tResult;
    }
}    // namespace moris
