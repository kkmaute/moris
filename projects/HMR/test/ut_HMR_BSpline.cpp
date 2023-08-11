/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_T_Matrix_Private.cpp
 *
 */

#include <catch.hpp>
#include "fn_HMR_bspline_shape.hpp"
#include "fn_norm.hpp"

namespace moris::hmr
{
    TEST_CASE( "HMR B-spline Shape Function", "[moris], [hmr], [bspline]" )
    {
        for ( uint tOrder = 1; tOrder <= 5; tOrder++ )
        {
            // points where the function is tested
            Matrix< DDRMat > tXi = {{-1, -0.5, 0, 0.5, 1}};

            // container for solution
            Matrix< DDRMat > tSolution;

            // fill solution vector with precomputed values
            switch ( tOrder )
            {
                case 1:
                {
                    Matrix< DDRMat > tSolution1 = {{ 1, 0 }, { 0.75, 0.25 }, { 0.5, 0.5 }, { 0.25, 0.75 }, { 0, 1 }};
                    tSolution = tSolution1;
                    break;
                }
                case 2:
                {
                    Matrix< DDRMat > tSolution3 = {{ 0.5, 0.5, 0 }, { 0.28125, 0.6875, 0.03125 },
                            { 0.125, 0.75, 0.125 }, { 0.03125, 0.6875, 0.28125 }, { 0, 0.5, 0.5 }};

                    tSolution = tSolution3;
                    break;
                }
                case 3:
                {
                    Matrix< DDRMat > tSolution3 = {{ 0.5, 2, 0.5, 0 }, { 0.2109375, 1.8359375, 0.9453125, 0.0078125 },
                            { 0.0625, 1.4375, 1.4375, 0.0625 }, { 0.0078125, 0.9453125, 1.8359375, 0.2109375 },
                            { 0, 0.5, 2, 0.5 }};
                    tSolution = tSolution3 * ( 1. / 3. );
                    break;
                }
                case 4:
                {
                    Matrix< DDRMat > tSolution4 = {{ 0.125, 1.375, 1.375, 0.125, 0 },
                            { 0.03955078125, 0.974609375, 1.6826171875, 0.302734375, 0.00048828125 },
                            { 0.0078125, 0.59375, 1.796875, 0.59375, 0.0078125 },
                            { 0.00048828125, 0.302734375, 1.6826171875, 0.974609375, 0.03955078125 },
                            { 0, 0.125, 1.375, 1.375, 0.125 }};

                    tSolution = tSolution4 * ( 1. / 3. );
                    break;
                }
                case 5:
                {
                    Matrix< DDRMat > tSolution5 = {{ 0.025, 0.65, 1.65, 0.65, 0.025, 0 },
                            { 0.0059326171875, 0.3747314453125, 1.558935546875, 0.984228515625, 0.0761474609375, 0.0000244140625 },
                            { 0.00078125, 0.18515625, 1.3140625, 1.3140625, 0.18515625, 0.00078125 },
                            { 2.44140625e-05, 0.0761474609375, 0.984228515625, 1.558935546875, 0.3747314453125, 0.0059326171875 },
                            { 0, 0.025, 0.65, 1.65, 0.65, 0.025 }};

                    tSolution = tSolution5 * ( 1. / 3. );
                    break;
                }
            }

            // loop over all contributing functions
            for ( uint iBasisNumber = 0; iBasisNumber < tOrder + 1; iBasisNumber++ )
            {
                // reset error
                Matrix< DDRMat > tError( 5, 1, 0 );

                for ( uint iTestPoint = 0; iTestPoint < 5; iTestPoint++ )
                {
                    // save error into matrix
                    tError( iTestPoint ) = bspline_shape( tOrder, iBasisNumber, tXi( iTestPoint ) ) - tSolution( iTestPoint, iBasisNumber );
                }

                // test solution
                REQUIRE ( norm( tError ) < 1e-12 );
            }
        }
    }
}
