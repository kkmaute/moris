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
    template< class A, class B >
    auto histc( const Eigen::MatrixBase< A >& aMatrix, const Eigen::MatrixBase< B >& aEdges )
    {
        const bool                                            tIsRowVector = aMatrix.rows() == 1;
        Eigen::Matrix< uint, Eigen::Dynamic, Eigen::Dynamic > tCounts      = Eigen::Matrix< uint, Eigen::Dynamic, Eigen::Dynamic >::Zero( tIsRowVector ? 1 : aEdges.size(), tIsRowVector ? aEdges.size() : aMatrix.cols() );
        for ( Eigen::Index k = 1; k < aEdges.size(); ++k )
        {
            MORIS_ERROR( aEdges.derived().coeff( k ) >= aEdges.derived().coeff( k - 1 ), "histc requires sorted bin edges." );
        }
        if ( !aEdges.size() )
        {
            return tCounts;
        }
        for ( Eigen::Index j = 0; j < ( tIsRowVector ? 1 : aMatrix.cols() ); ++j )
        {
            for ( Eigen::Index i = 0; i < ( tIsRowVector ? aMatrix.cols() : aMatrix.rows() ); ++i )
            {
                const auto   tValue = tIsRowVector ? aMatrix( 0, i ) : aMatrix( i, j );
                Eigen::Index tBin   = -1;
                // The final bin contains only values equal to the last edge.
                if ( tValue == aEdges.derived().coeff( aEdges.size() - 1 ) )
                {
                    tBin = aEdges.size() - 1;
                }
                else if ( tValue >= aEdges.derived().coeff( 0 ) && tValue < aEdges.derived().coeff( aEdges.size() - 1 ) )
                {
                    Eigen::Index tLow = 0, tHigh = aEdges.size() - 1;
                    while ( tLow + 1 < tHigh )
                    {
                        const auto tMiddle = ( tLow + tHigh ) / 2;
                        if ( tValue < aEdges.derived().coeff( tMiddle ) )
                        {
                            tHigh = tMiddle;
                        }
                        else
                        {
                            tLow = tMiddle;
                        }
                    }
                    tBin = tLow;
                }
                if ( tBin >= 0 )
                {
                    if ( tIsRowVector )
                    {
                        ++tCounts( 0, tBin );
                    }
                    else
                    {
                        ++tCounts( tBin, j );
                    }
                }
            }
        }
        return tCounts;
    }
}    // namespace moris
