/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Integrator_Test.cpp
 *
 */

#include "catch.hpp"
#include "cl_MTK_Integrator.hpp" //FEM//INT//src
#include "fn_FEM_Check.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "Integrator", "[moris],[fem],[Integrator]" )
{

    SECTION( "Integrator : Space bar - Time bar" )
    {
        // space time integrator for space bar and time bar
        //------------------------------------------------------------------------------
        // create a space time integration rule
        mtk::Integration_Rule tFieldIntegRule( mtk::Geometry_Type::LINE,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::BAR_5,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::BAR_5);

        // create an integrator
        mtk::Integrator tFieldIntegrator( tFieldIntegRule );

        // space HEX2x2x2 for comparison
        //------------------------------------------------------------------------------
        // create a space integration rule
        mtk::Integration_Order tSpaceIntOrder  = mtk::Integration_Order::QUAD_5x5;
        mtk::Integration_Rule tSpaceIntegRule( mtk::Geometry_Type::QUAD,
                                          mtk::Integration_Type::GAUSS,
                                          tSpaceIntOrder,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::BAR_1);
        // create an integrator
        mtk::Integrator tSpaceIntegrator( tSpaceIntegRule );

        // check integrator
        //------------------------------------------------------------------------------
        // define an epsilon environment
        double tEpsilon = 1E-12;

        // check the number of points
        REQUIRE( tFieldIntegrator.get_number_of_points() == tSpaceIntegrator.get_number_of_points());

        // switch points and their weights, since they are in a different order in space time
        Matrix< DDRMat > tSpaceIntegPoints;
        tSpaceIntegrator.get_points( tSpaceIntegPoints );
        Matrix< DDRMat > tSpaceIntegWeights;
        tSpaceIntegrator.get_weights( tSpaceIntegWeights );
        Matrix< DDRMat > tSpaceIntegPoints2  = tSpaceIntegPoints;
        Matrix< DDRMat > tSpaceIntegWeights2 = tSpaceIntegWeights;
        switch ( tSpaceIntOrder )
        {
            case ( mtk::Integration_Order::QUAD_2x2 ) :
            {
                tSpaceIntegPoints2( 0, 2 ) = tSpaceIntegPoints( 0, 3 );
                tSpaceIntegPoints2( 1, 2 ) = tSpaceIntegPoints( 1, 3 );
                tSpaceIntegPoints2( 0, 3 ) = tSpaceIntegPoints( 0, 2 );
                tSpaceIntegPoints2( 1, 3 ) = tSpaceIntegPoints( 1, 2 );

                tSpaceIntegWeights2( 2 ) = tSpaceIntegWeights( 3 );
                tSpaceIntegWeights2( 3 ) = tSpaceIntegWeights( 2 );
                break;
             }
            case ( mtk::Integration_Order::QUAD_3x3 ) :
            {
                tSpaceIntegPoints2( 0, 1 ) = tSpaceIntegPoints( 0, 4 );
                tSpaceIntegPoints2( 1, 1 ) = tSpaceIntegPoints( 1, 4 );
                tSpaceIntegPoints2( 0, 2 ) = tSpaceIntegPoints( 0, 1 );
                tSpaceIntegPoints2( 1, 2 ) = tSpaceIntegPoints( 1, 1 );
                tSpaceIntegPoints2( 0, 3 ) = tSpaceIntegPoints( 0, 7 );
                tSpaceIntegPoints2( 1, 3 ) = tSpaceIntegPoints( 1, 7 );
                tSpaceIntegPoints2( 0, 4 ) = tSpaceIntegPoints( 0, 8 );
                tSpaceIntegPoints2( 1, 4 ) = tSpaceIntegPoints( 1, 8 );
                tSpaceIntegPoints2( 0, 6 ) = tSpaceIntegPoints( 0, 3 );
                tSpaceIntegPoints2( 1, 6 ) = tSpaceIntegPoints( 1, 3 );
                tSpaceIntegPoints2( 0, 7 ) = tSpaceIntegPoints( 0, 6 );
                tSpaceIntegPoints2( 1, 7 ) = tSpaceIntegPoints( 1, 6 );
                tSpaceIntegPoints2( 0, 8 ) = tSpaceIntegPoints( 0, 2 );
                tSpaceIntegPoints2( 1, 8 ) = tSpaceIntegPoints( 1, 2 );

                tSpaceIntegWeights2( 1 ) = tSpaceIntegWeights( 4 );
                tSpaceIntegWeights2( 2 ) = tSpaceIntegWeights( 1 );
                tSpaceIntegWeights2( 3 ) = tSpaceIntegWeights( 7 );
                tSpaceIntegWeights2( 4 ) = tSpaceIntegWeights( 8 );
                tSpaceIntegWeights2( 6 ) = tSpaceIntegWeights( 3 );
                tSpaceIntegWeights2( 7 ) = tSpaceIntegWeights( 6 );
                tSpaceIntegWeights2( 8 ) = tSpaceIntegWeights( 2 );
                break;
            }
            case ( mtk::Integration_Order::QUAD_4x4 ) :
            {
                break;
            }
            case ( mtk::Integration_Order::QUAD_5x5 ) :
            {
                break;
            }
            default :
            {
                MORIS_ERROR( false, "Unknown integration order.");
                break;
            }
        }

        // check the points coordinates
        bool tCheckPoints = true;
        Matrix< DDRMat > tFieldIntegPoints;
        tFieldIntegrator.get_points( tFieldIntegPoints );

        for ( uint i = 0; i < tFieldIntegrator.get_number_of_points(); i++)
        {
            for( uint j = 0; j < 2; j++ )
            {
                tCheckPoints = tCheckPoints && ( std::abs( tFieldIntegPoints( j, i ) - tSpaceIntegPoints2( j, i ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckPoints );

        // check the points weights
        bool tCheckWeights = true;
        Matrix< DDRMat > tFieldIntegWeights;
        tFieldIntegrator.get_weights( tFieldIntegWeights );

        for ( uint i = 0; i < tFieldIntegrator.get_number_of_points(); i++)
        {
            tCheckWeights = tCheckWeights && ( std::abs( tFieldIntegWeights( i ) - tSpaceIntegWeights2( i )/2 ) < tEpsilon );
        }
        REQUIRE( tCheckWeights );
    }

    SECTION( "Integrator : Space quad - Time bar" )
    {
        // space time integrator for space QUAD2x2 and time bar2
        //------------------------------------------------------------------------------
        // create a space time integration rule
        mtk::Integration_Rule tFieldIntegRule( mtk::Geometry_Type::QUAD,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::QUAD_3x3,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::BAR_3 );

        // create an integrator
        mtk::Integrator tFieldIntegrator( tFieldIntegRule );

        // space HEX2x2x2 for comparison
        //------------------------------------------------------------------------------
        // create a space integration rule
        mtk::Integration_Rule tSpaceIntegRule( mtk::Geometry_Type::HEX,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::HEX_3x3x3,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::BAR_1 );
        // create an integrator
        mtk::Integrator tSpaceIntegrator( tSpaceIntegRule );

        // check integrator
        //------------------------------------------------------------------------------
        // define an epsilon environment
        double tEpsilon = 1E-12;

        // check the number of points
        REQUIRE( tFieldIntegrator.get_number_of_points() == tSpaceIntegrator.get_number_of_points());

        // check the points coordinates
        bool tCheckPoints = true;
        Matrix< DDRMat > tFieldIntegPoints;
        tFieldIntegrator.get_points( tFieldIntegPoints );
        Matrix< DDRMat > tSpaceIntegPoints;
        tSpaceIntegrator.get_points( tSpaceIntegPoints );
        for ( uint i = 0; i < tFieldIntegrator.get_number_of_points(); i++)
        {
            for( uint j = 0; j < 3; j++ )
            {
                tCheckPoints = tCheckPoints && ( std::abs( tFieldIntegPoints( j, i ) - tSpaceIntegPoints( j, i ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckPoints );

        // check the points weights
        bool tCheckWeights = true;
        Matrix< DDRMat > tFieldIntegWeights;
        tFieldIntegrator.get_weights( tFieldIntegWeights );
        Matrix< DDRMat > tSpaceIntegWeights;
        tSpaceIntegrator.get_weights( tSpaceIntegWeights );
        for ( uint i = 0; i < tFieldIntegrator.get_number_of_points(); i++)
        {
            tCheckWeights = tCheckWeights && ( std::abs( tFieldIntegWeights( i ) - tSpaceIntegWeights( i )/2 ) < tEpsilon );
        }
        REQUIRE( tCheckWeights );
    }

    SECTION( "Integrator : Space hex - Time bar" )
    {
        // space time integrator for space HEX2x2x2 and time bar2
        //------------------------------------------------------------------------------
        // create a space time integration rule
        mtk::Integration_Rule tFieldIntegRule( mtk::Geometry_Type::HEX,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::HEX_2x2x2,
                                          mtk::Integration_Type::GAUSS,
                                          mtk::Integration_Order::BAR_2);

        // create an integrator
        mtk::Integrator tFieldIntegrator( tFieldIntegRule );

        // check integrator
        //------------------------------------------------------------------------------

        // check the number of points
        REQUIRE( tFieldIntegrator.get_number_of_points() == 16);

        // check the points coords?
        Matrix< DDRMat > tFieldIntegPoints;
        tFieldIntegrator.get_points( tFieldIntegPoints );

        // check the points weights?
        Matrix< DDRMat > tFieldIntegWeights;
        tFieldIntegrator.get_weights( tFieldIntegWeights );
    }
}
