#include <string>
#include <catch.hpp>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp" //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "Integrator", "[moris],[fem],[Integrator]" )
{

    SECTION( "Integrator : Space bar - Time bar" )
    {
        // space time integrator for space bar and time bar
        //------------------------------------------------------------------------------
        // create a space time integration rule
        Integration_Rule tFieldIntegRule( mtk::Geometry_Type::LINE,
                                          Integration_Type::GAUSS,
                                          Integration_Order::BAR_5,
                                          Integration_Type::GAUSS,
                                          Integration_Order::BAR_5);

        // create an integrator
        Integrator tFieldIntegrator( tFieldIntegRule );

        // space HEX2x2x2 for comparison
        //------------------------------------------------------------------------------
        // create a space integration rule
        Integration_Order tSpaceIntOrder  = Integration_Order::QUAD_5x5;
        Integration_Rule tSpaceIntegRule( mtk::Geometry_Type::QUAD,
                                          Integration_Type::GAUSS,
                                          tSpaceIntOrder);
        // create an integrator
        Integrator tSpaceIntegrator( tSpaceIntegRule );

        // check integrator
        //------------------------------------------------------------------------------
        // define an epsilon environment
        double tEpsilon = 1E-12;

        // check the number of dimensions
        REQUIRE( tFieldIntegrator.get_number_of_dimensions() == tSpaceIntegrator.get_number_of_dimensions());

        // check the number of points
        REQUIRE( tFieldIntegrator.get_number_of_points() == tSpaceIntegrator.get_number_of_points());


        // switch points and their weights, since they are in a different order in space time
        Matrix< DDRMat > tSpaceIntegPoint   = tSpaceIntegrator.get_points();
        Matrix< DDRMat > tSpaceIntegWeight  = tSpaceIntegrator.get_weights();
        Matrix< DDRMat > tSpaceIntegPoint2  = tSpaceIntegPoint;
        Matrix< DDRMat > tSpaceIntegWeight2 = tSpaceIntegWeight;
        switch ( tSpaceIntOrder )
        {
            case ( Integration_Order::QUAD_2x2 ) :
            {
                tSpaceIntegPoint2( 0, 2 ) = tSpaceIntegPoint( 0, 3 );
                tSpaceIntegPoint2( 1, 2 ) = tSpaceIntegPoint( 1, 3 );
                tSpaceIntegPoint2( 0, 3 ) = tSpaceIntegPoint( 0, 2 );
                tSpaceIntegPoint2( 1, 3 ) = tSpaceIntegPoint( 1, 2 );

                tSpaceIntegWeight2( 2 ) = tSpaceIntegWeight( 3 );
                tSpaceIntegWeight2( 3 ) = tSpaceIntegWeight( 2 );
                break;
             }
            case ( Integration_Order::QUAD_3x3 ) :
            {
                tSpaceIntegPoint2( 0, 1 ) = tSpaceIntegPoint( 0, 4 );
                tSpaceIntegPoint2( 1, 1 ) = tSpaceIntegPoint( 1, 4 );
                tSpaceIntegPoint2( 0, 2 ) = tSpaceIntegPoint( 0, 1 );
                tSpaceIntegPoint2( 1, 2 ) = tSpaceIntegPoint( 1, 1 );
                tSpaceIntegPoint2( 0, 3 ) = tSpaceIntegPoint( 0, 7 );
                tSpaceIntegPoint2( 1, 3 ) = tSpaceIntegPoint( 1, 7 );
                tSpaceIntegPoint2( 0, 4 ) = tSpaceIntegPoint( 0, 8 );
                tSpaceIntegPoint2( 1, 4 ) = tSpaceIntegPoint( 1, 8 );
                tSpaceIntegPoint2( 0, 6 ) = tSpaceIntegPoint( 0, 3 );
                tSpaceIntegPoint2( 1, 6 ) = tSpaceIntegPoint( 1, 3 );
                tSpaceIntegPoint2( 0, 7 ) = tSpaceIntegPoint( 0, 6 );
                tSpaceIntegPoint2( 1, 7 ) = tSpaceIntegPoint( 1, 6 );
                tSpaceIntegPoint2( 0, 8 ) = tSpaceIntegPoint( 0, 2 );
                tSpaceIntegPoint2( 1, 8 ) = tSpaceIntegPoint( 1, 2 );

                tSpaceIntegWeight2( 1 ) = tSpaceIntegWeight( 4 );
                tSpaceIntegWeight2( 2 ) = tSpaceIntegWeight( 1 );
                tSpaceIntegWeight2( 3 ) = tSpaceIntegWeight( 7 );
                tSpaceIntegWeight2( 4 ) = tSpaceIntegWeight( 8 );
                tSpaceIntegWeight2( 6 ) = tSpaceIntegWeight( 3 );
                tSpaceIntegWeight2( 7 ) = tSpaceIntegWeight( 6 );
                tSpaceIntegWeight2( 8 ) = tSpaceIntegWeight( 2 );
                break;
            }
            case ( Integration_Order::QUAD_4x4 ) :
            {
                break;
            }
            case ( Integration_Order::QUAD_5x5 ) :
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
        Matrix< DDRMat > tSpaceTimeIntegPoint = tFieldIntegrator.get_points();

        for ( uint i = 0; i < tFieldIntegrator.get_number_of_points(); i++)
        {
            for( uint j = 0; j < tFieldIntegrator.get_number_of_dimensions(); j++ )
            {
                tCheckPoints = tCheckPoints && ( ( tSpaceTimeIntegPoint( j, i ) - tSpaceIntegPoint2( j, i )) < tEpsilon );
            }
        }
        REQUIRE( tCheckPoints );

        // check the points weights
        bool tCheckWeights = true;
        Matrix< DDRMat > tSpaceTimeIntegWeight = tFieldIntegrator.get_weights();

        for ( uint i = 0; i < tFieldIntegrator.get_number_of_points(); i++)
        {
            tCheckWeights = tCheckWeights && ( ( tSpaceTimeIntegWeight( i ) - tSpaceIntegWeight2( i )) < tEpsilon );
        }
        REQUIRE( tCheckWeights );
    }


    SECTION( "Integrator : Space quad - Time bar" )
    {
        // space time integrator for space QUAD2x2 and time bar2
        //------------------------------------------------------------------------------
        // create a space time integration rule
        Integration_Rule tFieldIntegRule( mtk::Geometry_Type::QUAD,
                                          Integration_Type::GAUSS,
                                          Integration_Order::QUAD_2x2,
                                          Integration_Type::GAUSS,
                                          Integration_Order::BAR_2);

        // create an integrator
        Integrator tFieldIntegrator( tFieldIntegRule );

        // space HEX2x2x2 for comparison
        //------------------------------------------------------------------------------
        // create a space integration rule
        Integration_Rule tSpaceIntegRule( mtk::Geometry_Type::HEX,
                                          Integration_Type::GAUSS,
                                          Integration_Order::HEX_2x2x2);
        // create an integrator
        Integrator tSpaceIntegrator( tSpaceIntegRule );

        // check integrator
        //------------------------------------------------------------------------------
        // define an epsilon environment
        double tEpsilon = 1E-12;

        // check the number of dimensions
        REQUIRE( tFieldIntegrator.get_number_of_dimensions() == tSpaceIntegrator.get_number_of_dimensions());

        // check the number of points
        REQUIRE( tFieldIntegrator.get_number_of_points() == tSpaceIntegrator.get_number_of_points());

        // check the points coordinates
        bool tCheckPoints = true;
        Matrix< DDRMat > tSpaceTimeIntegPoint = tFieldIntegrator.get_points();
        Matrix< DDRMat > tSpaceIntegPoint     = tSpaceIntegrator.get_points();
        for ( uint i = 0; i < tFieldIntegrator.get_number_of_points(); i++)
        {
            for( uint j = 0; j < tFieldIntegrator.get_number_of_dimensions(); j++ )
            {
                tCheckPoints = tCheckPoints && ( ( tSpaceTimeIntegPoint( j, i ) - tSpaceIntegPoint( j, i )) < tEpsilon );
            }
        }
        REQUIRE( tCheckPoints );

        // check the points weights
        bool tCheckWeights = true;
        Matrix< DDRMat > tSpaceTimeIntegWeight = tFieldIntegrator.get_weights();
        Matrix< DDRMat > tSpaceIntegWeight     = tSpaceIntegrator.get_weights();
        for ( uint i = 0; i < tFieldIntegrator.get_number_of_points(); i++)
        {
            tCheckWeights = tCheckWeights && ( ( tSpaceTimeIntegWeight( i ) - tSpaceIntegWeight( i )) < tEpsilon );
        }
        REQUIRE( tCheckWeights );

    }


    SECTION( "Integrator : Space hex - Time bar" )
    {
        // space time integrator for space HEX2x2x2 and time bar2
        //------------------------------------------------------------------------------
        // create a space time integration rule
        Integration_Rule tFieldIntegRule( mtk::Geometry_Type::HEX,
                                          Integration_Type::GAUSS,
                                          Integration_Order::HEX_2x2x2,
                                          Integration_Type::GAUSS,
                                          Integration_Order::BAR_2);

        // create an integrator
        Integrator tFieldIntegrator( tFieldIntegRule );

        // check integrator
        //------------------------------------------------------------------------------

        // check the number of dimensions
        REQUIRE( tFieldIntegrator.get_number_of_dimensions() == 4 );

        // check the number of points
        REQUIRE( tFieldIntegrator.get_number_of_points() == 16);

        // check the points coords?
        Matrix< DDRMat > tSpaceTimeIntegPoints = tFieldIntegrator.get_points();

        // check the points weights?
        Matrix< DDRMat > tSpaceTimeIntegWeight = tFieldIntegrator.get_weights();
    }
}


