/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Geometry_Interpolator_Test.cpp
 *
 */

#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr

using namespace moris;
using namespace fem;

TEST_CASE( "Geometry_Interpolator", "[moris],[fem],[GeoInterpolator]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    SECTION( "Geometry Interpolator : 1D space time" )
    {
        // space and time geometry interpolator
        //------------------------------------------------------------------------------
        // create a bar2 space element
        Matrix< DDRMat > tXHat( 2, 1 );
        tXHat( 0, 0 ) = -1.0;
        tXHat( 1, 0 ) =  3.0;

        // create a bar2 time element
        Matrix< DDRMat > tTHat( 2, 1 );
        tTHat( 0 ) = 0.0;
        tTHat( 1 ) = 5.0;

        // define the quad4 space element in the param space xi
        Matrix< DDRMat > tXiHat( 2, 1 );
        tXiHat( 0, 0 ) = -1.0;
        tXiHat( 1, 0 ) =  1.0;

        //define the line time element in the param space tau
        Matrix< DDRMat > tTauHat( 2, 1 );
        tTauHat( 0 ) = -1.0;
        tTauHat( 1 ) =  1.0;

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeoInterpRule(
                mtk::Geometry_Type::LINE,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeoInterpolator( tGeoInterpRule );

        //set the coefficients xHat, tHat
        tGeoInterpolator.set_coeff( tXHat, tTHat );

        //set the coefficients xHat, tHat
        tGeoInterpolator.set_space_param_coeff( tXiHat );
        tGeoInterpolator.set_time_param_coeff( tTauHat );

        // create an evaluation point xi, tau
        Matrix< DDRMat > tXi( 1, 1, 0.35 );
        Matrix< DDRMat > tTau( 1, 1, 0.70 );

        tGeoInterpolator.set_space( tXi );
        tGeoInterpolator.set_time( tTau );

        // check methods of the geometry interpolator
        //------------------------------------------------------------------------------
        // check space interpolation order
        REQUIRE( tGeoInterpolator.get_space_interpolation_order() == mtk::Interpolation_Order::LINEAR );

        // check time interpolation order
        REQUIRE( tGeoInterpolator.get_time_interpolation_order() == mtk::Interpolation_Order::LINEAR );

        // check space interpolation type
        REQUIRE( tGeoInterpolator.get_space_interpolation_type() == mtk::Interpolation_Type::LAGRANGE );

        // check time interpolation type
        REQUIRE( tGeoInterpolator.get_time_interpolation_type() == mtk::Interpolation_Type::LAGRANGE );

        // check number of space dimensions
        REQUIRE( tGeoInterpolator.get_number_of_space_dimensions() == 1 );

        // check number of time dimensions
        REQUIRE( tGeoInterpolator.get_number_of_time_dimensions() == 1 );

        // check number of space bases
        REQUIRE( tGeoInterpolator.get_number_of_space_bases() == 2 );

        // check number of time bases
        REQUIRE( tGeoInterpolator.get_number_of_time_bases() == 2 );

        // check space coefficients
        Matrix< DDRMat > tXHatCheck = tGeoInterpolator.get_space_coeff();
        bool tXHatCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_dimensions(); i++)
        {
            tXHatCheckBool = tXHatCheckBool && ( std::abs( tXHatCheck( i ) - tXHat( i ) ) < tEpsilon );
        }
        REQUIRE( tXHatCheckBool );

        // check time coefficients
        Matrix< DDRMat > tTHatCheck = tGeoInterpolator.get_time_coeff();
        bool tTHatCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            tTHatCheckBool = tTHatCheckBool && ( std::abs( tTHatCheck( i ) - tTHat( i ) ) < tEpsilon );
        }
        REQUIRE( tTHatCheckBool );

        // check space jacobian
        Matrix< DDRMat > tdNdXi = tGeoInterpolator.dNdXi();
        Matrix< DDRMat > tSpaceJt = tGeoInterpolator.space_jacobian();
        Matrix< DDRMat > tSpaceJtMatlab( 1, 1, 2.0 );

        bool tSpaceJCheck = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_dimensions(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_space_dimensions(); j++)
            {
                tSpaceJCheck = tSpaceJCheck && ( std::abs( tSpaceJtMatlab( i, j ) - tSpaceJt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceJCheck );

        // check time jacobian
        Matrix< DDRMat > tdNdTau = tGeoInterpolator.dNdTau();
        Matrix< DDRMat > tTimeJt = tGeoInterpolator.time_jacobian();
        Matrix< DDRMat > tTimeJtMatlab( 1, 1, 2.5 );

        bool tTimeJCheck = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_time_dimensions(); j++)
            {
                tTimeJCheck = tTimeJCheck && ( std::abs( tTimeJtMatlab( i, j ) - tTimeJt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeJCheck );

        // check space jacobian and matrices for second derivatives
        Matrix< DDRMat > tSpaceKt, tSpaceLt;
        Matrix< DDRMat > td2NdXi2 = tGeoInterpolator.d2NdXi2();
        tGeoInterpolator.space_jacobian_and_matrices_for_second_derivatives( tSpaceJt,
                tSpaceKt,
                tSpaceLt,
                tdNdXi,
                td2NdXi2 );
        Matrix< DDRMat > tSpaceKtMatlab( 1, 1, 0.0 );
        Matrix< DDRMat > tSpaceLtMatlab( 1, 1, 4.0 );

        bool tSpaceKCheck = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_dimensions(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_space_dimensions(); j++)
            {
                tSpaceKCheck = tSpaceKCheck && ( std::abs( tSpaceKtMatlab( i, j ) - tSpaceKt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceKCheck );

        bool tSpaceLCheck = true;
        for ( uint i = 0; i < 1; i++)
        {
            for ( uint j = 0; j < 1; j++)
            {
                tSpaceLCheck = tSpaceLCheck && ( std::abs( tSpaceLtMatlab( i, j ) - tSpaceLt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceLCheck );

        // check time jacobian and matrices for second derivatives
        Matrix< DDRMat > tTimeKt, tTimeLt;
        Matrix< DDRMat > td2NdTau2 = tGeoInterpolator.d2NdTau2();
        tGeoInterpolator.time_jacobian_and_matrices_for_second_derivatives( tTimeJt,
                tTimeKt,
                tTimeLt,
                tdNdTau,
                td2NdTau2 );
        Matrix< DDRMat > tTimeKtMatlab( 1, 1, 0.0 );
        Matrix< DDRMat > tTimeLtMatlab( 1, 1, 6.25 );

        bool tTimeKCheck = true;
        for ( uint i = 0; i < 1; i++)
        {
            for ( uint j = 0; j < 1; j++)
            {
                tTimeKCheck = tTimeKCheck && ( std::abs( tTimeKtMatlab( i, j ) - tTimeKt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeKCheck );

        bool tTimeLCheck = true;
        for ( uint i = 0; i < 1; i++)
        {
            for ( uint j = 0; j < 1; j++)
            {
                tTimeLCheck = tTimeLCheck && ( std::abs( tTimeLtMatlab( i, j ) - tTimeLt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeLCheck );

        // check space interpolation
        Matrix< DDRMat > tx = tGeoInterpolator.valx();
        Matrix< DDRMat > txMatlab( 1, 1, 1.7 );
        bool txCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_dimensions(); i++)
        {
            txCheckBool = txCheckBool && ( std::abs( txMatlab( i ) - tx( i ) ) < tEpsilon );
        }
        REQUIRE( txCheckBool );

        // check time interpolation
        Matrix< DDRMat > tt = tGeoInterpolator.valt();
        Matrix< DDRMat > ttMatlab( 1, 1, 4.25 );
        bool ttCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            ttCheckBool = ttCheckBool && ( std::abs( ttMatlab( i ) - tt( i ) ) < tEpsilon );
        }
        REQUIRE( ttCheckBool );
    }

    SECTION( "Geometry Interpolator : 2D space time" )
    {
        // space and time geometry interpolator
        //------------------------------------------------------------------------------
        // define a quad4 space element in the phys space x
        Matrix< DDRMat > tXHat( 4, 2 );
        tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
        tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
        tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
        tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;

        //define a line time element in the phys space t
        Matrix< DDRMat > tTHat( 2, 1 );
        tTHat( 0 ) = 0.0;
        tTHat( 1 ) = 5.0;

        // define the quad4 space element in the param space xi
        Matrix< DDRMat > tXiHat( 4, 2 );
        tXiHat( 0, 0 ) = -1.0; tXiHat( 0, 1 ) = -1.0;
        tXiHat( 1, 0 ) =  1.0; tXiHat( 1, 1 ) = -1.0;
        tXiHat( 2, 0 ) =  1.0; tXiHat( 2, 1 ) =  1.0;
        tXiHat( 3, 0 ) = -1.0; tXiHat( 3, 1 ) =  1.0;

        //define the line time element in the param space tau
        Matrix< DDRMat > tTauHat( 2, 1 );
        tTauHat( 0 ) = -1.0;
        tTauHat( 1 ) =  1.0;

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeoInterpRule( mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeoInterpolator( tGeoInterpRule, mtk::CellShape::GENERAL, true );

        //set the coefficients xHat, tHat
        tGeoInterpolator.set_space_coeff( tXHat );
        tGeoInterpolator.set_time_coeff( tTHat );

        //set the coefficients xiHat, tauHat
        tGeoInterpolator.set_space_param_coeff( tXiHat );
        tGeoInterpolator.set_time_param_coeff( tTauHat );

        // create an evaluation point xi, tau
        Matrix< DDRMat > tXi( 2, 1 );
        tXi( 0, 0 ) =  0.35;
        tXi( 1, 0 ) = -0.25;
        Matrix< DDRMat > tTau( 1, 1, 0.70 );
        tGeoInterpolator.set_space( tXi );
        tGeoInterpolator.set_time( tTau );

        // check methods of the geometry interpolator
        //------------------------------------------------------------------------------
        // check space interpolation order
        REQUIRE( tGeoInterpolator.get_space_interpolation_order() == mtk::Interpolation_Order::LINEAR );

        // check time interpolation order
        REQUIRE( tGeoInterpolator.get_time_interpolation_order() == mtk::Interpolation_Order::LINEAR );

        // check space interpolation type
        REQUIRE( tGeoInterpolator.get_space_interpolation_type() == mtk::Interpolation_Type::LAGRANGE );

        // check time interpolation type
        REQUIRE( tGeoInterpolator.get_time_interpolation_type() == mtk::Interpolation_Type::LAGRANGE );

        // check number of space dimensions
        REQUIRE( tGeoInterpolator.get_number_of_space_dimensions() == 2 );

        // check number of time dimensions
        REQUIRE( tGeoInterpolator.get_number_of_time_dimensions() == 1 );

        // check number of space bases
        REQUIRE( tGeoInterpolator.get_number_of_space_bases() == 4 );

        // check number of time bases
        REQUIRE( tGeoInterpolator.get_number_of_time_bases() == 2 );

        // check space coefficients
        Matrix< DDRMat > tXHatCheck = tGeoInterpolator.get_space_coeff();
        bool tXHatCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_bases(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_space_dimensions(); j++)
            {
                tXHatCheckBool = tXHatCheckBool && ( std::abs( tXHatCheck( i, j ) - tXHat( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tXHatCheckBool );

        // check time coefficients
        Matrix< DDRMat > tTHatCheck = tGeoInterpolator.get_time_coeff();
        bool tTHatCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            tTHatCheckBool = tTHatCheckBool && ( std::abs( tTHatCheck( i ) - tTHat( i ) ) < tEpsilon );
        }
        REQUIRE( tTHatCheckBool );

        // check space jacobian
        Matrix< DDRMat > tdNdXi = tGeoInterpolator.dNdXi();
        Matrix< DDRMat > tSpaceJt = tGeoInterpolator.space_jacobian();
        Matrix< DDRMat > tSpaceJtMatlab = {{ 1.593750000000000, 0.531250000000000 },
                { 0.668750000000000, 1.456250000000000}};

        bool tSpaceJCheck = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_dimensions(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_space_dimensions(); j++)
            {
                tSpaceJCheck = tSpaceJCheck && ( std::abs( tSpaceJtMatlab( i, j ) - tSpaceJt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceJCheck );

        // check time jacobian
        Matrix< DDRMat > tdNdTau = tGeoInterpolator.dNdTau();
        Matrix< DDRMat > tTimeJt = tGeoInterpolator.time_jacobian();
        Matrix< DDRMat > tTimeJtMatlab( 1, 1, 2.5 );

        bool tTimeJCheck = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_time_dimensions(); j++)
            {
                tTimeJCheck = tTimeJCheck && ( std::abs( tTimeJtMatlab( i, j ) - tTimeJt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeJCheck );

        // check space jacobian and matrices for second derivatives
        Matrix< DDRMat > tSpaceKt, tSpaceLt;
        Matrix< DDRMat > td2NdXi2 = tGeoInterpolator.d2NdXi2();
        tGeoInterpolator.space_jacobian_and_matrices_for_second_derivatives( tSpaceJt,
                tSpaceKt,
                tSpaceLt,
                tdNdXi,
                td2NdXi2 );
        Matrix< DDRMat > tSpaceKtMatlab = {
                { 0.000000000000000,  0.000000000000000 },
                { 0.000000000000000,  0.000000000000000 },
                { 0.125000000000000, -0.125000000000000 }
        };
        Matrix< DDRMat > tSpaceLtMatlab = {
                { 2.540039062500000,   0.282226562500000,   1.693359375000000 },
                { 0.447226562500000,   2.120664062500000,   1.947734375000000 },
                { 1.065820312500000,   0.773632812500000,   2.676171875000000 }
        };

        bool tSpaceKCheck = true;
        for ( uint i = 0; i < 3; i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_space_dimensions(); j++)
            {
                tSpaceKCheck = tSpaceKCheck && ( std::abs( tSpaceKtMatlab( i, j ) - tSpaceKt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceKCheck );

        bool tSpaceLCheck = true;
        for ( uint i = 0; i < 3; i++)
        {
            for ( uint j = 0; j < 3; j++)
            {
                tSpaceLCheck = tSpaceLCheck && ( std::abs( tSpaceLtMatlab( i, j ) - tSpaceLt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceLCheck );

        // check time jacobian and matrices for second derivatives
        Matrix< DDRMat > tTimeKt, tTimeLt;
        Matrix< DDRMat > td2NdTau2 = tGeoInterpolator.d2NdTau2();
        tGeoInterpolator.time_jacobian_and_matrices_for_second_derivatives( tTimeJt,
                tTimeKt,
                tTimeLt,
                tdNdTau,
                td2NdTau2 );
        Matrix< DDRMat > tTimeKtMatlab( 1, 1, 0.0 );
        Matrix< DDRMat > tTimeLtMatlab( 1, 1, 6.25 );

        bool tTimeKCheck = true;
        for ( uint i = 0; i < 1; i++)
        {
            for ( uint j = 0; j < 1; j++)
            {
                tTimeKCheck = tTimeKCheck && ( std::abs( tTimeKtMatlab( i, j ) - tTimeKt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeKCheck );

        bool tTimeLCheck = true;
        for ( uint i = 0; i < 1; i++)
        {
            for ( uint j = 0; j < 1; j++)
            {
                tTimeLCheck = tTimeLCheck && ( std::abs( tTimeLtMatlab( i, j ) - tTimeLt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeLCheck );

        // check space interpolation
        Matrix< DDRMat > tx = tGeoInterpolator.valx();
        Matrix< DDRMat > txMatlab = {
                { 2.526562500000000,   1.935937500000000 },
        };

        bool txCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_dimensions(); i++)
        {
            txCheckBool = txCheckBool && ( std::abs( txMatlab( i ) - tx( i ) ) < tEpsilon );
        }
        REQUIRE( txCheckBool );

        // check time interpolation
        Matrix< DDRMat > tt = tGeoInterpolator.valt();
        Matrix< DDRMat > ttMatlab( 1, 1, 4.25 );
        bool ttCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            ttCheckBool = ttCheckBool && ( std::abs( ttMatlab( i ) - tt( i ) ) < tEpsilon );
        }
        REQUIRE( ttCheckBool );

        // check the parametric coordinates of a space side
        Vector< moris_index > tAllSpaceOrdinals = { 0, 1, 2, 3 };

        // loop over the sides
        for( uint iSide = 0; iSide < tAllSpaceOrdinals.size(); iSide++ )
        {
            // get the treated side ordinal
            moris_index tSpaceOrdinal = tAllSpaceOrdinals( iSide );

            // get side coords for check
            Matrix< DDRMat > tSpaceParamCoords;
            Matrix< DDRMat > tSpacePhysCoords;
            switch( tSpaceOrdinal )
            {
                case 0 :
                    tSpaceParamCoords = {{ -1.0, -1.0 }, { 1.0, -1.0 }};
                    tSpacePhysCoords  = {{  0.0,  0.0 }, { 3.0,  1.25 }};
                    break;
                case 1 :
                    tSpaceParamCoords = {{ 1.0, -1.0 }, { 1.0, 1.0 }};
                    tSpacePhysCoords  = {{ 3.0, 1.25 }, { 4.5, 4.0 }};
                    break;
                case 2 :
                    tSpaceParamCoords = {{ 1.0, 1.0 }, { -1.0, 1.0 }};
                    tSpacePhysCoords  = {{ 4.5, 4.0 }, { 1.0, 3.25 }};
                    break;
                case 3 :
                    tSpaceParamCoords = {{ -1.0, 1.0 }, { -1.0, -1.0 }};
                    tSpacePhysCoords  = {{ 1.0, 3.25 }, { 0.0, 0.0 }};
                    break;
                default:
                    MORIS_ASSERT( false, " wrong side ordinal ");
                    break;
            }
        }
    }

    SECTION( "Geometry Interpolator : 3D space time" )
    {
        // space and time geometry interpolator
        //------------------------------------------------------------------------------
        // define a hex8 space element in the phys space x
        Matrix< DDRMat > tXHat( 8, 3 );
        tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;  tXHat( 0, 2 ) = 0.0;
        tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25; tXHat( 1, 2 ) = 0.0;
        tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;  tXHat( 2, 2 ) = 0.0;
        tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25; tXHat( 3, 2 ) = 0.0;
        tXHat( 4, 0 ) = 1.0; tXHat( 4, 1 ) = 1.0;  tXHat( 4, 2 ) = 1.0;
        tXHat( 5, 0 ) = 4.0; tXHat( 5, 1 ) = 2.25; tXHat( 5, 2 ) = 1.5;
        tXHat( 6, 0 ) = 5.5; tXHat( 6, 1 ) = 5.0;  tXHat( 6, 2 ) = 2.0;
        tXHat( 7, 0 ) = 2.0; tXHat( 7, 1 ) = 4.25; tXHat( 7, 2 ) = 3.0;

        //define a line time element in the phys space t
        Matrix< DDRMat > tTHat( 2, 1 );
        tTHat( 0 ) = 0.0;
        tTHat( 1 ) = 5.0;

        // define the quad4 space element in the param space xi
        Matrix< DDRMat > tXiHat( 3, 8 );
        tXiHat( 0, 0 ) = -1.000000; tXiHat( 1, 0 ) = -1.000000; tXiHat( 2, 0 ) = -1.000000;
        tXiHat( 0, 1 ) =  1.000000; tXiHat( 1, 1 ) = -1.000000; tXiHat( 2, 1 ) = -1.000000;
        tXiHat( 0, 2 ) =  1.000000; tXiHat( 1, 2 ) =  1.000000; tXiHat( 2, 2 ) = -1.000000;
        tXiHat( 0, 3 ) = -1.000000; tXiHat( 1, 3 ) =  1.000000; tXiHat( 2, 3 ) = -1.000000;
        tXiHat( 0, 4 ) = -1.000000; tXiHat( 1, 4 ) = -1.000000; tXiHat( 2, 4 ) =  1.000000;
        tXiHat( 0, 5 ) =  1.000000; tXiHat( 1, 5 ) = -1.000000; tXiHat( 2, 5 ) =  1.000000;
        tXiHat( 0, 6 ) =  1.000000; tXiHat( 1, 6 ) =  1.000000; tXiHat( 2, 6 ) =  1.000000;
        tXiHat( 0, 7 ) = -1.000000; tXiHat( 1, 7 ) =  1.000000; tXiHat( 2, 7 ) =  1.000000;
        tXiHat = trans( tXiHat );

        //define the line time element in the param space tau
        Matrix< DDRMat > tTauHat( 2, 1 );
        tTauHat( 0 ) = -1.0;
        tTauHat( 1 ) =  1.0;

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeomInterpRule(
                mtk::Geometry_Type::HEX,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeoInterpolator( tGeomInterpRule );

        //set the coefficients xHat, tHat
        tGeoInterpolator.set_coeff( tXHat, tTHat );

        //set the coefficients xiHat, tauHat
        tGeoInterpolator.set_space_param_coeff( tXiHat );
        tGeoInterpolator.set_time_param_coeff( tTauHat );

        // create an evaluation point xi, tau
        Matrix< DDRMat > tXi( 3, 1 );
        tXi( 0, 0 ) =  0.35; tXi( 1, 0 ) = -0.25; tXi( 2, 0 ) = 0.74;
        Matrix< DDRMat > tTau( 1, 1, 0.70);

        tGeoInterpolator.set_space( tXi );
        tGeoInterpolator.set_time( tTau );

        // check methods of the geometry interpolator
        //------------------------------------------------------------------------------
        // check space interpolation order
        REQUIRE( tGeoInterpolator.get_space_interpolation_order() == mtk::Interpolation_Order::LINEAR );

        // check time interpolation order
        REQUIRE( tGeoInterpolator.get_time_interpolation_order() == mtk::Interpolation_Order::LINEAR );

        // check space interpolation type
        REQUIRE( tGeoInterpolator.get_space_interpolation_type() == mtk::Interpolation_Type::LAGRANGE );

        // check time interpolation type
        REQUIRE( tGeoInterpolator.get_time_interpolation_type() == mtk::Interpolation_Type::LAGRANGE );

        // check number of space dimensions
        REQUIRE( tGeoInterpolator.get_number_of_space_dimensions() == 3 );

        // check number of time dimensions
        REQUIRE( tGeoInterpolator.get_number_of_time_dimensions() == 1 );

        // check number of space bases
        REQUIRE( tGeoInterpolator.get_number_of_space_bases() == 8 );

        // check number of time bases
        REQUIRE( tGeoInterpolator.get_number_of_time_bases() == 2 );

        // check space coefficients
        Matrix< DDRMat > tXHatCheck = tGeoInterpolator.get_space_coeff();
        bool tXHatCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_bases(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_space_dimensions(); j++)
            {
                tXHatCheckBool = tXHatCheckBool && ( std::abs( tXHatCheck( i, j ) - tXHat( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tXHatCheckBool );

        // check time coefficients
        Matrix< DDRMat > tTHatCheck = tGeoInterpolator.get_time_coeff();
        bool tTHatCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            tTHatCheckBool = tTHatCheckBool && ( std::abs( tTHatCheck( i ) - tTHat( i ) ) < tEpsilon );
        }
        REQUIRE( tTHatCheckBool );

        // check space jacobian
        Matrix< DDRMat > tdNdXi = tGeoInterpolator.dNdXi();
        Matrix< DDRMat > tSpaceJt = tGeoInterpolator.space_jacobian();
        Matrix< DDRMat > tSpaceJtMatlab = {
                { 1.593750000000000,   0.531250000000000,  -0.027187500000000},
                { 0.668750000000000,   1.456250000000000,   0.429562500000000},
                { 0.500000000000000,   0.500000000000000,   0.853906250000000}};

        bool tSpaceJCheck = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_dimensions(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_space_dimensions(); j++)
            {
                tSpaceJCheck = tSpaceJCheck && ( std::abs( tSpaceJtMatlab( i, j ) - tSpaceJt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceJCheck );

        // check time jacobian
        Matrix< DDRMat > tdNdTau = tGeoInterpolator.dNdTau();
        Matrix< DDRMat > tTimeJt = tGeoInterpolator.time_jacobian();
        Matrix< DDRMat > tTimeJtMatlab( 1, 1, 2.5 );

        bool tTimeJCheck = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_time_dimensions(); j++)
            {
                tTimeJCheck = tTimeJCheck && ( std::abs( tTimeJtMatlab( i, j ) - tTimeJt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeJCheck );

        // check space jacobian and matrices for second derivatives
        Matrix< DDRMat > tSpaceKt, tSpaceLt;
        Matrix< DDRMat > td2NdXi2 = tGeoInterpolator.d2NdXi2();
        tGeoInterpolator.space_jacobian_and_matrices_for_second_derivatives( tSpaceJt,
                tSpaceKt,
                tSpaceLt,
                tdNdXi,
                td2NdXi2 );
        Matrix< DDRMat > tSpaceKtMatlab = {
                { 0,                   0,                    0 },
                { 0,                   0,                    0 },
                { 0,                   0,                    0 },
                { 0,                  -0.000000000000000,    0.246875000000000 },
                { 0,                   0,                   -0.015625000000000 },
                { 0.125000000000000,  -0.125000000000000,   -0.326250000000000}};
        Matrix< DDRMat > tSpaceLtMatlab = {
                { 2.540039062500000,   0.282226562500000,   0.000739160156250,  -0.028886718750000,  -0.086660156250000,   1.693359375000000 },
                { 0.447226562500000,   2.120664062500000,   0.184523941406250,   1.251100781250000,   0.574539843750000,   1.947734375000001 },
                { 0.250000000000000,   0.250000000000000,   0.729155883789063,   0.853906250000000,   0.853906250000000,   0.500000000000000 },
                { 0.334375000000000,   0.728125000000000,   0.366806103515625,   1.458282226562500,   0.785831054687500,   1.062500000000000 },
                { 0.796875000000000,   0.265625000000000,  -0.023215576171875,   0.440043945312500,   1.347319335937500,   1.062500000000000 },
                { 1.065820312500000,   0.773632812500000,  -0.011678730468750,   0.188613281250000,   0.666433593750000,   2.676171875000000 }
        };

        bool tSpaceKCheck = true;
        for ( uint i = 0; i < 3; i++)
        {
            for ( uint j = 0; j < tGeoInterpolator.get_number_of_space_dimensions(); j++)
            {
                tSpaceKCheck = tSpaceKCheck && ( std::abs( tSpaceKtMatlab( i, j ) - tSpaceKt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceKCheck );

        bool tSpaceLCheck = true;
        for ( uint i = 0; i < 3; i++)
        {
            for ( uint j = 0; j < 3; j++)
            {
                tSpaceLCheck = tSpaceLCheck && ( std::abs( tSpaceLtMatlab( i, j ) - tSpaceLt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tSpaceLCheck );

        // check time jacobian and matrices for second derivatives
        Matrix< DDRMat > tTimeKt, tTimeLt;
        Matrix< DDRMat > td2NdTau2 = tGeoInterpolator.d2NdTau2();
        tGeoInterpolator.time_jacobian_and_matrices_for_second_derivatives( tTimeJt,
                tTimeKt,
                tTimeLt,
                tdNdTau,
                td2NdTau2 );
        Matrix< DDRMat > tTimeKtMatlab( 1, 1, 0.0 );
        Matrix< DDRMat > tTimeLtMatlab( 1, 1, 6.25 );

        bool tTimeKCheck = true;
        for ( uint i = 0; i < 1; i++)
        {
            for ( uint j = 0; j < 1; j++)
            {
                tTimeKCheck = tTimeKCheck && ( std::abs( tTimeKtMatlab( i, j ) - tTimeKt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeKCheck );

        bool tTimeLCheck = true;
        for ( uint i = 0; i < 1; i++)
        {
            for ( uint j = 0; j < 1; j++)
            {
                tTimeLCheck = tTimeLCheck && ( std::abs( tTimeLtMatlab( i, j ) - tTimeLt( i, j ) ) < tEpsilon );
            }
        }
        REQUIRE( tTimeLCheck );

        // check space interpolation
        Matrix< DDRMat > tx = tGeoInterpolator.valx();
        Matrix< DDRMat > txMatlab = {
                { 3.396562500000000,   2.805937500000000,   1.485796875000000 },
        };

        bool txCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_space_dimensions(); i++)
        {
            txCheckBool = txCheckBool && ( std::abs( txMatlab( i ) - tx( i ) ) < tEpsilon );
        }
        REQUIRE( txCheckBool );

        // check time interpolation
        Matrix< DDRMat > tt = tGeoInterpolator.valt();
        Matrix< DDRMat > ttMatlab( 1, 1, 4.25 );
        bool ttCheckBool = true;
        for ( uint i = 0; i < tGeoInterpolator.get_number_of_time_dimensions(); i++)
        {
            ttCheckBool = ttCheckBool && ( std::abs( ttMatlab( i ) - tt( i ) ) < tEpsilon );
        }
        REQUIRE( ttCheckBool );

        // check the parametric coordinates of a space side
        Vector< moris_index > tAllSpaceOrdinals = { 0, 1, 2, 3, 4, 5 };

        // loop over the sides
        for( uint iSide = 0; iSide < tAllSpaceOrdinals.size(); iSide++ )
        {
            // get the treated side ordinal
            moris_index tSpaceOrdinal = tAllSpaceOrdinals( iSide );

            // get side coords for check
            Matrix< DDRMat > tSpaceParamCoords;
            Matrix< DDRMat > tSpacePhysCoords;
            switch( tSpaceOrdinal )
            {
                case 0 :
                {
                    tSpaceParamCoords = {{ -1.0, -1.0, -1.0 }, { 1.0, -1.0, -1.0 }, { 1.0, -1.0, 1.0 }, { -1.0, -1.0, 1.0 }};
                    tSpacePhysCoords = {{ 0.0, 0.0, 0.0 }, { 3.0, 1.25, 0.0 }, { 4.0, 2.25, 1.5 }, { 1.0, 1.0, 1.0 }};
                    break;
                }
                case 1 :
                {
                    tSpaceParamCoords = {{ 1.0, -1.0, -1.0 }, { 1.0, 1.0, -1.0 }, { 1.0, 1.0, 1.0 }, { 1.0, -1.0, 1.0 }};
                    tSpacePhysCoords  = {{ 3.0, 1.25, 0.0 }, { 4.5, 4.0, 0.0 }, { 5.5, 5.0, 2.0 }, { 4.0, 2.25, 1.5}};
                    break;
                }
                case 2 :
                {
                    tSpaceParamCoords = {{ 1.0, 1.0, -1.0 }, { -1.0, 1.0, -1.0 }, { -1.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }};
                    tSpacePhysCoords  = {{ 4.5, 4.0, 0.0 }, { 1.0, 3.25, 0.0 }, { 2.0, 4.25, 3.0 }, { 5.5, 5.0, 2.0 }};
                    break;
                }
                case 3 :
                {
                    tSpaceParamCoords = {{ -1.0, -1.0, -1.0 }, { -1.0, -1.0, 1.0 }, { -1.0, 1.0, 1.0 }, { -1.0, 1.0, -1.0 }};
                    tSpacePhysCoords  = {{ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 }, { 2.0, 4.25, 3.0 }, { 1.0, 3.25, 0.0 }};
                    break;
                }
                case 4 :
                {
                    tSpaceParamCoords = {{ -1.0, -1.0, -1.0 }, { -1.0, 1.0, -1.0 }, { 1.0, 1.0, -1.0 }, { 1.0, -1.0, -1.0 }};
                    tSpacePhysCoords  = {{ 0.0, 0.0, 0.0 }, { 1.0, 3.25, 0.0 }, { 4.5, 4.0, 0.0 }, { 3.0, 1.25, 0.0 }};
                    break;
                }
                case 5 :
                {
                    tSpaceParamCoords = {{ -1.0, -1.0, 1.0 }, { 1.0, -1.0, 1.0 }, { 1.0, 1.0, 1.0 }, { -1.0, 1.0, 1.0 }};
                    tSpacePhysCoords  = {{ 1.0, 1.0, 1.0 }, { 4.0, 2.25, 1.5 }, { 5.5, 5.0, 2.0 }, { 2.0, 4.25, 3.0 }};
                    break;
                }
                default:
                {
                    MORIS_ASSERT( false, " wrong side ordinal "); break;
                }
            }
        }
    }

    //------------------------------------------------------------------------------
}
