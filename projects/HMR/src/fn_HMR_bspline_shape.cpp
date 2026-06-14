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

    real
    bspline_shape(
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
        for ( uint iPolyOrder = 0; iPolyOrder <= aOrder; ++iPolyOrder )
        {
            if ( tDeltaXi( iPolyOrder + aBasisNumber ) <= aXi && aXi < tDeltaXi( iPolyOrder + aBasisNumber + 1 ) )
            {
                tN( iPolyOrder ) = 1.0;
            }
        }

        // loop over all orders
        for ( uint iPolyOrder = 1; iPolyOrder <= aOrder; ++iPolyOrder )
        {
            // copy values of tN into old matrix
            Matrix< DDRMat > tNold( tN );

            // loop over all contributions
            for ( uint i = 0; i <= aOrder - iPolyOrder; ++i )
            {
                // help values
                real tA = aXi - tDeltaXi( i + aBasisNumber );
                real tB = tDeltaXi( i + aBasisNumber + iPolyOrder + 1 ) - aXi;

                tN( i ) = 0.5 * ( tA * tNold( i ) + tB * ( tNold( i + 1 ) ) ) / ( (real)iPolyOrder );
            }
        }

        // first value in entry is shape value
        return tN( 0 );
    }

    // -----------------------------------------------------------------------------------------------------------------

    real
    bspline_shape_extended(
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
    //* NEW *
    // -----------------------------------------------------------------------------------------------------------------

    real
    expline_1D_constant(
            const moris_index aBasisNumber,
            const real aX ) // in range [0,1]
    {
        // explicit polynomials depending on position of basis function relative to knot span
        switch ( aBasisNumber ) 
        {
            case 0: // left-most basis function
            {
                return 1.0;
            }

            default: // not supported in knot span
            {
                return 0.0;
            }

        } // end switch: location of BF relative to knot span

    } // end function: expline_1D_constant()

    // -----------------------------------------------------------------------------------------------------------------

    real
    expline_1D_linear(
            const moris_index aBasisNumber,
            const real aX ) // in range [0,1]
    {
        // explicit polynomials depending on position of basis function relative to knot span
        switch ( aBasisNumber ) 
        {
            case 0: // left-most basis function
            {
                return 1.0 - aX;
            }

            case 1: // right-most basis function
            {
                return aX;
            }

            default: // not supported in knot span
            {
                return 0.0;
            }

        } // end switch: location of BF relative to knot span

    } // end function: expline_1D_linear()

    // -----------------------------------------------------------------------------------------------------------------

    real
    expline_1D_quadratic(
            const moris_index aBasisNumber,
            const real aX ) // in range [0,1]
    {
        // explicit polynomials depending on position of basis function relative to knot span
        switch ( aBasisNumber ) 
        {
            case 0: // left-most basis function
            {
                return 0.5 * ( aX - 1.0 ) * ( aX - 1.0 );
            }

            case 1:
            {
                return -aX * ( 0.5 * aX - 1.0 ) - 0.5 * ( aX + 1.0 ) * ( aX - 1.0 );
            }

            case 2: // right-most basis function
            {
                return 0.5 * aX * aX;
            }

            default: // not supported in knot span
            {
                return 0.0;
            }

        } // end switch: location of BF relative to knot span

    } // end function: expline_1D_quadratic()

    // -----------------------------------------------------------------------------------------------------------------

    real
    expline_1D_cubic(
            const moris_index aBasisNumber,
            const real aX ) // in range [0,1]
    {
        // explicit polynomials depending on position of basis function relative to knot span
        switch ( aBasisNumber ) 
        {
            case 0: // left-most basis function
            {
                return -0.5 * ( aX - 1.0 ) * ( 1.0 / 3.0 ) * ( aX - 1.0 ) *( aX - 1.0 );
            }

            case 1:
            {
                // (xi/3 - 2/3)*(xi*(xi/2 - 1) + (xi/2 + 1/2)*(xi - 1)) + (xi/2 - 1/2)*(xi/3 + 2/3)*(xi - 1) 
                return ( aX - 2.0 ) * ( aX * ( 0.5 * aX - 1.0 ) ) / 3.0
                        + 0.5 * ( aX + 1.0 ) * ( aX - 1.0 )
                        + 0.5 * ( aX - 1.0 ) * ( aX + 2.0 ) * ( aX - 1.0 ) / 3.0;
            }

            case 2:
            {
                // - (xi^2*(xi/3 - 1))/2 - (xi/3 + 1/3)*(xi*(xi/2 - 1) + (xi/2 + 1/2)*(xi - 1))
                return  - 0.5 * ( aX * aX * ( aX / 3.0 - 1.0 ) ) 
                        - ( aX + 1.0 ) * ( aX * ( 0.5 * aX - 1.0 ) ) / 3.0 
                        + 0.5 * ( aX + 1.0 ) * ( aX - 1.0 );
            }

            case 3: // right-most basis function
            {
                return aX * aX * aX / 6.0; // xi^3/6
            }

            default: // not supported in knot span
            {
                return 0.0;
            }

        } // end switch: location of BF relative to knot span

    } // end function: expline_1D_cubic()

    // -----------------------------------------------------------------------------------------------------------------

    // FIXME: this function should be templated against the number of dimensions and polynomial order, as these are known at compile time
    real
    expline_1D(
            const uint aOrder,
            const moris_index aBasisNumber, // ordered going from left-most to right-most basis function supported in know span
            const real aX ) // in range [0,1]
    {
        // explicit polynomials depending on order ...
        switch (aOrder) 
        {
            case 0: // constant
            {
                return expline_1D_constant( aBasisNumber, aX );
            }

            case 1: // linear
            {
                return expline_1D_linear( aBasisNumber, aX );
            }

            case 2: // quadratic
            {
                return expline_1D_quadratic( aBasisNumber, aX );
            }

            case 3: // cubic
            {
                return expline_1D_cubic( aBasisNumber, aX );
            }

            default:
                MORIS_ERROR( false, "expline_1D() - only implemented for orders p in {0,1,2,3}." );
                return 0.0;
        } 
    }

    // -----------------------------------------------------------------------------------------------------------------

    real
    eval_spline(
            const uint aNumDims, // d
            const uint aOrder,   // p
            Vector< moris_index > const & aRelativeIJK,
            Matrix< DDRMat > const & aXi ) // in range [-1,1]^d
    {
        // first dimension
        real tX = 0.5 * ( aXi( 0 ) + 1.0 ); // rescale to range [0,1] used in the expline functions
        real tVal = expline_1D( aOrder, aRelativeIJK( 0 ), tX );

        for ( uint iDim = 1; iDim < aNumDims; iDim++ )
        {
            tX = 0.5 * ( aXi( iDim ) + 1.0 ); // rescale to range [0,1] used in the expline functions
            tVal *= expline_1D( aOrder, aRelativeIJK( iDim ), tX );
        }

        return tVal;
    }

    // -----------------------------------------------------------------------------------------------------------------

    // FIXME: this is an inefficient way to perform the evaluation
    real
    eval( 
            const uint aNumDims,
            const uint aOrder,
            const uint aBfLevel,
            const luint* aBfIJK,
            const uint aElementLevel,
            const luint* aElementIJK,
            const Matrix< DDRMat > & aXi ) // in range [-1,1]^d
    {
        // sanity check input
        MORIS_ASSERT( aElementLevel >= aBfLevel, "HMR::eval() - Lagrange element coarser than B-spline BF. An evaluation should not happen in this setting." );

        // create modifiable copies
        Matrix< DDRMat > tXi( aXi );
        uint tCurrentLevel = aElementLevel;
        Vector< luint > tCurrentIJK( aNumDims );
        for ( uint iDim = 0; iDim < aNumDims; iDim++ )
        {
            tCurrentIJK( iDim ) = aElementIJK[ iDim ];
        }

        // perform the re-mapping to the BF-level
        while ( tCurrentLevel > aBfLevel )
        {
            // do this in 1D for each independent dimension
            for ( uint iDim = 0; iDim < aNumDims; iDim++ )
            {
                if ( tCurrentIJK( iDim ) % 2 == 0 ) // even I
                {
                    tCurrentIJK( iDim ) = tCurrentIJK( iDim ) / 2;
                    tXi( iDim ) = 0.5 * tXi( iDim ) - 0.5;
                }
                else // uneven I
                {
                    tCurrentIJK( iDim ) = ( tCurrentIJK( iDim ) - 1 ) / 2;
                    tXi( iDim ) = 0.5 * tXi( iDim ) + 0.5;
                }
            }

            // move up
            tCurrentLevel--;

        } // end while: 

        // get the relative IJk between BF and element
        Vector< moris_index > tRelativeIJK( aNumDims );
        moris_index tAuraCorrection = (moris_index)aOrder * ( std::pow( 2, (moris_index)aBfLevel ) - 1 );
        for ( uint iDim = 0; iDim < aNumDims; iDim++ )
        {
            tRelativeIJK( iDim ) = (moris_index)( aBfIJK[ iDim ] - tCurrentIJK( iDim ) ) - tAuraCorrection;
        }

        // evaluate and return the B-spline basis function
        return eval_spline( aNumDims, aOrder, tRelativeIJK, tXi );

    }

    // -----------------------------------------------------------------------------------------------------------------


}    // namespace moris::hmr
