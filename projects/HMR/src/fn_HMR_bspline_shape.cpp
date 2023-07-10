/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_bspline_shape.cpp
 *
 */

#include "fn_HMR_bspline_shape.hpp"
#include "cl_Matrix.hpp"

namespace moris::hmr
{
    
    // -----------------------------------------------------------------------------------------------------------------
    
    real bspline_shape(
            uint aOrder,
            uint aBasisNumber,
            real aXi )
    {
        // max number of entries in lookup table
        uint tSteps = 2 * ( aOrder + 1 );

        // temporary matrix that contains B-Spline segments
        Matrix< DDRMat > tDeltaXi( tSteps, 1, 0 );
        for ( uint i = 0; i < tSteps; ++i )
        {
            tDeltaXi( i ) = ( ( (real)i ) - ( (real)aOrder ) ) * 2.0 - 1.0;
        }

        // temporary matrix that contains evaluated values
        Matrix< DDRMat > tN( aOrder + 1, 1, 0 );

        // initialize zero order values
        for ( uint i = 0; i <= aOrder; ++i )
        {
            if ( tDeltaXi( i + aBasisNumber ) <= aXi && aXi < tDeltaXi( i + aBasisNumber + 1 ) )
            {
                tN( i ) = 1.0;
            }
        }

        // loop over all orders
        for ( uint r = 1; r <= aOrder; ++r )
        {
            // copy values of tN into old matrix
            Matrix< DDRMat > tNold( tN );

            // loop over all contributions
            for ( uint i = 0; i <= aOrder - r; ++i )
            {
                // help values
                real tA = aXi - tDeltaXi( i + aBasisNumber );
                real tB = tDeltaXi( i + aBasisNumber + r + 1 ) - aXi;

                tN( i ) = 0.5 * ( tA * tNold( i ) + tB * ( tNold( i + 1 ) ) ) / ( (real)r );
            }
        }

        // first value in entry is shape value
        return tN( 0 );
    }

    // -----------------------------------------------------------------------------------------------------------------

    real bspline_shape_extended(
            uint aOrder,
            uint aBasisNumber,
            real aXi )
    {
        switch ( aOrder )
        {
            // linear interpolation
            case 1:
            {
                // local ordering of basis function
                switch ( aBasisNumber )
                {
                    case 0:
                    {
                        return 0.5 * ( 1.0 - aXi );
                    }
                    case 1:
                    {
                        return 0.5 * ( 1.0 + aXi );
                    }
                    default:
                    {
                        MORIS_ERROR( false, "The specified local basis %u is not implemented", aBasisNumber );
                        return 0.0;
                    }
                }
            }

            default:
            {
                MORIS_ERROR( false, "The specified order %u is not implemented", aOrder );
                return 0.0;
            }
        }
    }

    // -----------------------------------------------------------------------------------------------------------------
    
}
