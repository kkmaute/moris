/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Side_Geometry_Interpolation_Test.cpp
 *
 */

#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr
#include "cl_MTK_Integrator.hpp" //MTK/sr

using namespace moris;
using namespace fem;

TEST_CASE("Side_Geometry_Interpolation : QUAD4 - QUAD9 - QUAD16", "[moris],[fem],[SideGeoInterp_QUAD]")
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // get side surface for comparison
    Matrix< DDRMat > tSideSurfaceExact( 4, 1 );
    tSideSurfaceExact( 0 ) = 3.250000000000000;
    tSideSurfaceExact( 1 ) = 3.132491021535417;
    tSideSurfaceExact( 2 ) = 3.579455265819088;
    tSideSurfaceExact( 3 ) = 3.400367627183861;

    // get side normals for comparison
    Matrix< DDRMat > tSideNormalsExact( 2, 4 );
    tSideNormalsExact( 0, 0 ) =  0.3846153846153846;
    tSideNormalsExact( 1, 0 ) = -0.9230769230769231;
    tSideNormalsExact( 0, 1 ) =  0.8778955729143844;
    tSideNormalsExact( 1, 1 ) = -0.4788521306805733;
    tSideNormalsExact( 0, 2 ) = -0.2095290887308734;
    tSideNormalsExact( 1, 2 ) =  0.9778024140774094;
    tSideNormalsExact( 0, 3 ) = -0.9557790087219500;
    tSideNormalsExact( 1, 3 ) =  0.2940858488375231;

    // define a list of side ordinals
    Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2, 3 };

    // loop over the interpolation order
    for( uint iOrder = 0; iOrder < 3; iOrder++ )
    {
        // define a QUAD of order tSpaceInterpolationOrder
        Matrix< DDRMat > tXHat;
        Matrix< DDRMat > tXiHat;
        mtk::Interpolation_Order tSpaceInterpolationOrder;
        Matrix< DDSMat > tSideNodes;

        // loop over the interpolation orders
        switch( iOrder )
        {
            case 0 :
            {
                // define a QUAD4 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                tXHat = {{ 0.0, 0.0 }, { 3.0, 1.25 }, { 4.5, 4.0 }, { 1.0, 3.25} };
                // define a QUAD4 in the parametric space, i.e. space param coordinates xiHat
                tXiHat = {{ -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 }};
                // mapping of nodes per side for HEX8
                tSideNodes = {{ 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }};
                break;
            }

            case 1 :
            {
                // define a QUAD9 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;
                tXHat = {{ 0.0,   0.0   }, { 3.0,  1.25  }, { 4.5, 4.0 }, { 1.0, 3.25},
                         { 1.5,   0.625 }, { 3.75, 2.625 },
                         { 2.75,  3.625 }, { 0.5,  1.625 },
                         { 2.125, 2.125 }};
                // define a QUAD9 in the parametric space, i.e. space param coordinates xiHat
                tXiHat.set_size( 2, 9 );
                tXiHat( 0, 0 ) = -1.000000; tXiHat( 1, 0 ) = -1.000000;
                tXiHat( 0, 1 ) =  1.000000; tXiHat( 1, 1 ) = -1.000000;
                tXiHat( 0, 2 ) =  1.000000; tXiHat( 1, 2 ) =  1.000000;
                tXiHat( 0, 3 ) = -1.000000; tXiHat( 1, 3 ) =  1.000000;
                tXiHat( 0, 4 ) =  0.000000; tXiHat( 1, 4 ) = -1.000000;
                tXiHat( 0, 5 ) =  1.000000; tXiHat( 1, 5 ) =  0.000000;
                tXiHat( 0, 6 ) =  0.000000; tXiHat( 1, 6 ) =  1.000000;
                tXiHat( 0, 7 ) = -1.000000; tXiHat( 1, 7 ) =  0.000000;
                tXiHat( 0, 8 ) =  0.000000; tXiHat( 1, 8 ) =  0.000000;
                tXiHat = trans( tXiHat );
                tSideNodes = { { 0, 1, 4 }, { 1, 2, 5 }, { 2, 3, 6 }, { 3, 0, 7 } };
                break;
            }
            case 2 :
            {
                // define a QUAD16 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::CUBIC;
                tXHat = {{0.000000000000000,  0.000000000000000}, {3.000000000000000,  1.250000000000000},
                         {4.500000000000000,  4.000000000000000}, {1.000000000000000,  3.250000000000000},
                         {1.000000000000000,  0.416666666666666}, {2.000000000000000,  0.833333333333333},
                         {3.500000000000000,  2.166666666666667}, {4.000000000000000,  3.083333333333333},
                         {3.333333333333333,  3.750000000000000}, {2.166666666666667,  3.500000000000000},
                         {0.666666666666666,  2.166666666666667}, {0.333333333333333,  1.083333333333333},
                         {1.388888888888889,  1.444444444444445}, {2.444444444444445,  1.805555555555556},
                         {2.888888888888889,  2.777777777777778}, {1.777777777777778,  2.472222222222222}};
                // define a QUAD16 in the parametric space, i.e. space param coordinates xiHat
                tXiHat.set_size( 2, 16 ); real c = 1.0/3.0;
                tXiHat( 0,  0 ) = -1.000000; tXiHat( 1,  0 ) = -1.000000;
                tXiHat( 0,  1 ) =  1.000000; tXiHat( 1,  1 ) = -1.000000;
                tXiHat( 0,  2 ) =  1.000000; tXiHat( 1,  2 ) =  1.000000;
                tXiHat( 0,  3 ) = -1.000000; tXiHat( 1,  3 ) =  1.000000;
                tXiHat( 0,  4 ) = -c;        tXiHat( 1,  4 ) = -1.000000;
                tXiHat( 0,  5 ) =  c;        tXiHat( 1,  5 ) = -1.000000;
                tXiHat( 0,  6 ) =  1.000000; tXiHat( 1,  6 ) = -c;
                tXiHat( 0,  7 ) =  1.000000; tXiHat( 1,  7 ) =  c;
                tXiHat( 0,  8 ) =  c;        tXiHat( 1,  8 ) =  1.000000;
                tXiHat( 0,  9 ) = -c;        tXiHat( 1,  9 ) =  1.000000;
                tXiHat( 0, 10 ) = -1.000000; tXiHat( 1, 10 ) =  c;
                tXiHat( 0, 11 ) = -1.000000; tXiHat( 1, 11 ) = -c;
                tXiHat( 0, 12 ) = -c;        tXiHat( 1, 12 ) = -c;
                tXiHat( 0, 13 ) =  c;        tXiHat( 1, 13 ) = -c;
                tXiHat( 0, 14 ) =  c;        tXiHat( 1, 14 ) =  c;
                tXiHat( 0, 15 ) = -c;        tXiHat( 1, 15 ) =  c;
                tXiHat = trans( tXiHat );
                tSideNodes = { { 0, 1,  4,  5 }, { 1, 2,  6,  7 }, { 2, 3,  8,  9 }, { 3, 0, 10, 11 } };
                break;
            }
        }

        //create a line time element
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};
        Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

        // side geometry type and space ordinal
        mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::LINE;

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tSideGeoInterpRule(
                tSideGeometryType,
                                                mtk::Interpolation_Type::LAGRANGE,
                                                tSpaceInterpolationOrder,
                                                mtk::Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tSideGeoInterp( tSideGeoInterpRule, mtk::CellShape::GENERAL, true );

        // create a side integration
        mtk::Integration_Rule tSideIntegRule( tSideGeometryType,
                                         mtk::Integration_Type::GAUSS,
                                         mtk::Integration_Order::BAR_5,
                                         mtk::Integration_Type::GAUSS,
                                         mtk::Integration_Order::BAR_1 );

        // create a side integrator
        mtk::Integrator tSideIntegrator( tSideIntegRule );

        //get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tSideIntegPoints;
        tSideIntegrator.get_points( tSideIntegPoints );
        Matrix< DDRMat > tSideIntegWeights;
        tSideIntegrator.get_weights( tSideIntegWeights );

        // booleans for checks
        bool tSideSurfaceCheck = true;

        // loop over the element sides
        for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
        {
            // get treated side ordinal
            moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

            // get the node ids associated to the follower side ordinal
            Matrix< DDSMat > tSideOrdNodes = tSideNodes.get_row( tSideOrdinal );

            // phys coords and param coords in IP param space for the follower side
            uint tNumSideNodes = tSideOrdNodes.numel();
            Matrix< DDRMat > tSideXHat( tNumSideNodes, 2 );
            Matrix< DDRMat > tSideXiHat( tNumSideNodes, 2 );
            for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
            {
                tSideXHat.get_row( iNode )  = tXHat.get_row( tSideOrdNodes( iNode ) );
                tSideXiHat.get_row( iNode ) = tXiHat.get_row( tSideOrdNodes( iNode ) );
            }
            //set the coefficients xHat, tHat
            tSideGeoInterp.set_space_coeff( tSideXHat );
            tSideGeoInterp.set_time_coeff( tTHat );

            // set the param coefficients xiHat, tauHat
            tSideGeoInterp.set_space_param_coeff( tSideXiHat );
            tSideGeoInterp.set_time_param_coeff( tTauHat );

            // init the side surface
            real tSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tLocalSideIntegPoint = tSideIntegPoints.get_column( iGP );

                // set the treated integration point location in the surface ref space for the geometry interp
                tSideGeoInterp.set_space_time( tLocalSideIntegPoint );

                // evaluate surfDetJ
                real tSurfDetJ = tSideGeoInterp.det_J();

                // add contribution to the surface
                tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );

                // get normal
                Matrix< DDRMat > tSideNormal;
                tSideGeoInterp.get_normal( tSideNormal );

                // check the normal
                bool tSideNormalCheck = true;
                for( uint i = 0; i < 2; i++ )
                {
                    tSideNormalCheck = tSideNormalCheck && ( std::abs( tSideNormal( i ) - tSideNormalsExact( i, iSide ) ) < tEpsilon );
                }
                REQUIRE( tSideNormalCheck );
            }

            // check the surface value
            tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );
        }
        // check surfDetJ and surfNormal values
        REQUIRE( tSideSurfaceCheck );

//        // time side geometry type and space ordinal------------------------------------
//        //------------------------------------------------------------------------------
//
//        Cell< moris::moris_index > tListOfTimeOrdinals = { 0, 1 };
//
//        // create a side integration
//        mtk::Integration_Rule tTimeSideIntegRule( mtk::Geometry_Type::QUAD,
//                                             mtk::Integration_Type::GAUSS,
//                                             mtk::Integration_Order::QUAD_3x3,
//                                             mtk::Integration_Type::GAUSS,
//                                             mtk::Integration_Order::BAR_1 );
//
//        // create a side integrator
//        mtk::Integrator tTimeSideIntegrator( tTimeSideIntegRule );
//
//        //get number of integration points, integration points and weights
//        uint             tNumOfTimeIntegPoints = tTimeSideIntegrator.get_number_of_points();
//        Matrix< DDRMat > tTimeSideIntegPoints  = tTimeSideIntegrator.get_points();
//        Matrix< DDRMat > tTimeSideIntegWeights = tTimeSideIntegrator.get_weights();
//
//        // boolean for surface check
//        bool tTimeSurfaceCheck = true;
//
//        // loop over the time side
//        for ( uint iTimeSide = 0; iTimeSide < tListOfTimeOrdinals.size(); iTimeSide++ )
//        {
//            // get the time ordinal
//            moris_index tTimeSideOrdinal = tListOfTimeOrdinals( iTimeSide );
//
//            // build and get the side time param coords
//            tElemGeoInterp.build_time_side_time_param_coeff( tTimeSideOrdinal );
//
//            // build and get the side time phys coords
//            tElemGeoInterp.build_time_side_time_phys_coeff( tTimeSideOrdinal );
//
//            // init the side surface
//            real tTimeSideSurface = 0;
//
//            // loop over the integration points
//            for ( uint iGP = 0; iGP < tNumOfTimeIntegPoints; iGP++ )
//            {
//                // get the treated integration point location in the surface ref space
//                Matrix< DDRMat > tTimeSideIntegPointI = tTimeSideIntegPoints.get_column( iGP );
//
//                // get the treated integration point location in the volume ref space
//                Matrix< DDRMat > tTimeVolIntegPointI = tElemGeoInterp.time_surf_val( tTimeSideIntegPointI );
//
//                // evaluate timesurfDetJ
//                real tSurfDetJ = tElemGeoInterp.time_surf_det_J( tTimeSideIntegPointI );
//
//                // add contribution to the surface
//                tTimeSideSurface = tTimeSideSurface + tSurfDetJ * tTimeSideIntegWeights( iGP );
//            }
//            // check the surface value
//            tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 8.5 ) < tEpsilon );
//        }
//        REQUIRE( tTimeSurfaceCheck );
    }
}

//------------------------------------------------------------------------------
TEST_CASE( "Side_Geometry_Interpolation : TRI3 - TRI6 - TRI10", "[moris],[fem],[SideGeoInterp_TRI]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // side surface for comparison
    Matrix< DDRMat > tSideSurfaceExact( 3, 1, 2.0 );
    tSideSurfaceExact( 1 ) = 2 * std::sqrt( 2 );

    // Side normals for comparison
    Matrix< DDRMat > tSideNormalsExact( 2, 3 );
    tSideNormalsExact( 0, 0 ) =  0.0;
    tSideNormalsExact( 1, 0 ) = -1.0;
    tSideNormalsExact( 0, 1 ) =  0.707106781186;
    tSideNormalsExact( 1, 1 ) =  0.707106781186; ;
    tSideNormalsExact( 0, 2 ) = -1.0;
    tSideNormalsExact( 1, 2 ) =  0.0;

    // loop over the interpolation order
    for( uint iOrder = 0; iOrder < 3; iOrder++ )
    {
        // define a TRI of order tSpaceInterpolationOrder
        Matrix< DDRMat > tXHat;
        Matrix< DDRMat > tXiHat;
        mtk::Interpolation_Order tSpaceInterpolationOrder;
        Matrix< DDSMat > tSideNodes;

        switch( iOrder )
        {
            case ( 0 ):
            {
                // define a TRI4 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                tXHat = {{ 0.0, 0.0 },
                         { 1.0, 0.0 },
                         { 0.0, 1.0 }};
                tXiHat = {{ 1.0, 0.0, 0.0 },
                          { 0.0, 1.0, 0.0 },
                          { 0.0, 0.0, 1.0 }};
                tSideNodes = {{ 0, 1 }, { 1, 2 }, { 2, 0 }};
                break;
            }
            case ( 1 ):
            {
                // define a TRI6 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;
                real t12 = 1.0/2.0;
                tXHat = {{ 0.0, 0.0 },
                         { 1.0, 0.0 },
                         { 0.0, 1.0 },
                         { t12, 0.0 },
                         { t12, t12 },
                         { 0.0, t12 }};
                tXiHat = {{ 1.000000000000000, 0.000000000000000, 0.000000000000000,
                            0.500000000000000, 0.000000000000000, 0.500000000000000 },
                          { 0.000000000000000, 1.000000000000000, 0.000000000000000,
                            0.500000000000000, 0.500000000000000, 0.000000000000000 },
                          { 0.000000000000000, 0.000000000000000, 1.000000000000000,
                            0.000000000000000, 0.500000000000000, 0.500000000000000 }};
                tXiHat = trans( tXiHat );
                tSideNodes = {{ 0, 1, 3 }, { 1, 2, 4 }, { 2, 0, 5 }};
                break;
            }
            case ( 2 ):
            {
                // define a TRI10 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::CUBIC;
                real t13 = 1.0/3.0;
                real t23 = 2.0/3.0;
                tXHat = {{ 0.0, 0.0 },
                         { 1.0, 0.0 },
                         { 0.0, 1.0 },
                         { t13, 0.0 },
                         { t23, 0.0 },
                         { t23, t13 },
                         { t13, t23 },
                         { 0.0, t23 },
                         { 0.0, t13 },
                         { t13, t13 }};
                tXiHat = {{ 1.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t23, t13 },
                          { 0.0, 1.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, t13 },
                          { 0.0, 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, t13 }};
                tXiHat = trans( tXiHat );
                tSideNodes = {{ 0, 1, 3, 4 }, { 1, 2, 5, 6 }, { 2, 0, 7, 8 }};
                break;
            }
        }

        // define a LINE2 in the physical space, i.e. time coordinates tHat
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 2.0 }};
        Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 } };

        // side geometry type and space ordinal
        mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::LINE;

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tSideGeoInterpRule( tSideGeometryType,
                                               mtk::Interpolation_Type::LAGRANGE,
                                               tSpaceInterpolationOrder,
                                               mtk::Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::LINEAR );

       // create a space and time geometry interpolator
       Geometry_Interpolator tSideGeoInterp( tSideGeoInterpRule, mtk::CellShape::GENERAL, true );

        // side space ordinal
        Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2 };

        // create a side integration
        mtk::Integration_Rule tSideIntegrationRule( tSideGeometryType,
                                               mtk::Integration_Type::GAUSS,
                                               mtk::Integration_Order::BAR_5,
                                               mtk::Integration_Type::GAUSS,
                                               mtk::Integration_Order::BAR_1 );

        // create a side integrator
        mtk::Integrator tSideIntegrator( tSideIntegrationRule );

        //get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tSideIntegPoints;
        tSideIntegrator.get_points( tSideIntegPoints );
        Matrix< DDRMat > tSideIntegWeights;
        tSideIntegrator.get_weights( tSideIntegWeights );

        // booleans for checks
        bool tSideSurfaceCheck = true;

        // loop over the sides
        for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
        {
            // get treated side ordinal
            moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

            // get the node ids associated to the follower side ordinal
            Matrix< DDSMat > tSideOrdNodes = tSideNodes.get_row( tSideOrdinal );

            // phys coords and param coords in IP param space for the follower side
            uint tNumSideNodes = tSideOrdNodes.numel();
            Matrix< DDRMat > tSideXHat( tNumSideNodes, 2 );
            Matrix< DDRMat > tSideXiHat( tNumSideNodes, 3 );
            for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
            {
                tSideXHat.get_row( iNode )  = tXHat.get_row( tSideOrdNodes( iNode ) );
                tSideXiHat.get_row( iNode ) = tXiHat.get_row( tSideOrdNodes( iNode ) );
            }

            //set the coefficients xHat, tHat
            tSideGeoInterp.set_space_coeff( tSideXHat );
            tSideGeoInterp.set_time_coeff( tTHat );

            // set the param coefficients xiHat, tauHat
            tSideGeoInterp.set_space_param_coeff( tSideXiHat );
            tSideGeoInterp.set_time_param_coeff( tTauHat );

            // init the side surface and the normal
            real tSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tLocalSideIntegPoint = tSideIntegPoints.get_column( iGP );

                // set the treated integration point location in the surface ref space for the geometry interp
                tSideGeoInterp.set_space_time( tLocalSideIntegPoint );

                // evaluate surfDetJ and normal
                real tSurfDetJ = tSideGeoInterp.det_J();

                // add contribution to the surface
                tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );

                // get normal
                Matrix< DDRMat > tSideNormal;
                tSideGeoInterp.get_normal( tSideNormal );

                // check the normal
                bool tSideNormalCheck = true;
                for( uint i = 0; i < 2; i++ )
                {
                    tSideNormalCheck = tSideNormalCheck && ( std::abs( tSideNormal( i ) - tSideNormalsExact( i, iSide ) ) < tEpsilon );
                }
                REQUIRE( tSideNormalCheck );
            }

            // check the surface value
            tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );
        }

        REQUIRE( tSideSurfaceCheck );

    }
}

//------------------------------------------------------------------------------
TEST_CASE( "Side_Geometry_Interpolation :  TET4 - TET10 - TET20", "[moris],[fem],[SideGeoInterp_TET]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

     // Side surface for comparison
     Matrix< DDRMat >   tSideSurfaceExact( 4, 1, 1.0 );
     tSideSurfaceExact( 1 ) = std::sqrt( 3 );

     // Side normals for comparison
     Matrix< DDRMat > tSideNormalsExact( 3, 4 );
     tSideNormalsExact( 0, 0 ) = -1.0;
     tSideNormalsExact( 1, 0 ) =  0.0;
     tSideNormalsExact( 2, 0 ) =  0.0;
     tSideNormalsExact( 0, 1 ) =  0.5773502691896;
     tSideNormalsExact( 1, 1 ) = -0.5773502691896;
     tSideNormalsExact( 2, 1 ) =  0.5773502691896;
     tSideNormalsExact( 0, 2 ) =  0.0;
     tSideNormalsExact( 1, 2 ) =  1.0;
     tSideNormalsExact( 2, 2 ) =  0.0;
     tSideNormalsExact( 0, 3 ) =  0.0;
     tSideNormalsExact( 1, 3 ) =  0.0;
     tSideNormalsExact( 2, 3 ) = -1.0;

     // switch on the TET interpolation order
     for ( uint iOrder = 0; iOrder < 3; iOrder++ )
     {
         // define a TET of order tSpaceInterpolationOrder
         Matrix< DDRMat > tXHat;
         Matrix< DDRMat > tXiHat;
         mtk::Interpolation_Order tSpaceInterpolationOrder;
         Matrix< DDSMat > tSideNodes;

         switch( iOrder )
         {
             case ( 0 ):
             {
                 // define a TET4 in the physical space, i.e. space coordinates xHat
                 tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                 tXHat = {{ 0.0,  0.0, 0.0 },
                          { 0.0, -1.0, 0.0 },
                          { 1.0,  0.0, 0.0 },
                          { 0.0,  0.0, 1.0 }};
                 tXiHat = {{ 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 },
                           { 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 },
                           { 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000 },
                           { 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000 }};
                 tXiHat = trans( tXiHat );
                 tSideNodes = {{ 0, 1, 3 }, { 1, 2, 3 }, { 0, 3, 2 }, { 0, 2, 1 }};
                 break;
             }
             case ( 1 ):
             {
                 // define a TET10 in the physical space, i.e. space coordinates xHat
                 tSpaceInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;
                 real t12 = 1.0/2.0;
                 tXHat = {{ 0.0,  0.0, 0.0 },
                          { 0.0, -1.0, 0.0 },
                          { 1.0,  0.0, 0.0 },
                          { 0.0,  0.0, 1.0 },
                          { 0.0, -t12, 0.0 },
                          { t12, -t12, 0.0 },
                          { t12,  0.0, 0.0 },
                          { 0.0,  0.0, t12 },
                          { 0.0, -t12, t12 },
                          { t12,  0.0, t12 } };
                 tXiHat ={{ 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000,
                            0.000000000000000, 0.500000000000000, 0.500000000000000, 0.000000000000000, 0.000000000000000 },
                          { 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000,
                            0.500000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000, 0.000000000000000 },
                          { 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000,
                            0.500000000000000, 0.500000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000 },
                          { 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000,
                            0.000000000000000, 0.000000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000 }};
                 tXiHat = trans( tXiHat );
                 tSideNodes = {{ 0, 1, 3, 4, 8, 7 }, { 1, 2, 3, 5, 9, 8 }, { 0, 3, 2, 7, 9, 6 }, { 0, 2, 1, 6, 5, 4 }};
                 break;
             }
             case ( 2 ):
             {
                 // define a TET20 in the physical space, i.e. space coordinates xHat
                 tSpaceInterpolationOrder = mtk::Interpolation_Order::CUBIC;
                 real t13 = 1.0/3.0;
                 real t23 = 2.0/3.0;
                 tXHat = {{ 0.0,  0.0, 0.0 },
                          { 0.0, -1.0, 0.0 },
                          { 1.0,  0.0, 0.0 },
                          { 0.0,  0.0, 1.0 },
                          { 0.0, -t13, 0.0 },
                          { 0.0, -t23, 0.0 },
                          { t13, -t23, 0.0 },
                          { t23, -t13, 0.0 },
                          { t13,  0.0, 0.0 },
                          { t23,  0.0, 0.0 },
                          { 0.0,  0.0, t13 },
                          { 0.0,  0.0, t23 },
                          { 0.0, -t23, t13 },
                          { 0.0, -t13, t23 },
                          { t23,  0.0, t13 },
                          { t13,  0.0, t23 },
                          { t13, -t13, 0.0 },
                          { 0.0, -t13, t13 },
                          { t13, -t13, t13 },
                          { t13,  0.0, t13 } };
                 tXiHat = {{ 1.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t23, t13, t23, t13, 0.0, 0.0, 0.0, 0.0, t13, t13, 0.0, t13 },
                            { 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t13, t13, 0.0 },
                            { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, 0.0, 0.0, 0.0, 0.0, t23, t13, t13, 0.0, t13, t13 },
                            { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, t13, t23, 0.0, t13, t13, t13 }};
                 tXiHat = trans( tXiHat );
                 tSideNodes = {{ 0, 1, 3,  4,  5, 12, 13, 11, 10, 17 }, { 1, 2, 3,  6,  7, 14, 15, 13, 12, 18 },
                               { 0, 3, 2, 10, 11, 15, 14,  9,  8, 19 }, { 0, 2, 1,  8,  9,  7,  6,  5,  4, 16 }};
                 break;
             }
         }

         // define a LINE2 in the physical space, i.e. time coordinates tHat
         Matrix< DDRMat > tTHat = {{ 0.0 }, { 2.0 }};
         Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 } };

         // side geometry type
         mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::TRI;

         // create a space and time geometry interpolation rule
         mtk::Interpolation_Rule tSideGeoInterpRule( tSideGeometryType,
                                                mtk::Interpolation_Type::LAGRANGE,
                                                tSpaceInterpolationOrder,
                                                mtk::Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

         // create a space and time geometry interpolator
         Geometry_Interpolator tSideGeoInterp( tSideGeoInterpRule, mtk::CellShape::GENERAL, true );

         // side space ordinal
         Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2, 3 };

         // create a side integration
         mtk::Integration_Rule tSideIntegrationRule( tSideGeometryType,
                                                mtk::Integration_Type::GAUSS,
                                                mtk::Integration_Order::TRI_6,
                                                mtk::Integration_Type::GAUSS,
                                                mtk::Integration_Order::BAR_1 );

         // create a side integrator
         mtk::Integrator tSideIntegrator( tSideIntegrationRule );

         //get number of integration points, integration points and weights
         uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
         Matrix< DDRMat > tSideIntegPoints;
         tSideIntegrator.get_points( tSideIntegPoints );
         Matrix< DDRMat > tSideIntegWeights;
         tSideIntegrator.get_weights( tSideIntegWeights );

         // booleans for checks
         bool tSideSurfaceCheck = true;

         // loop over the element sides
         for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
         {
             // get treated side ordinal
             moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

             // get the node ids associated to the follower side ordinal
             Matrix< DDSMat > tSideOrdNodes = tSideNodes.get_row( tSideOrdinal );

             // phys coords and param coords in IP param space for the follower side
             uint tNumSideNodes = tSideOrdNodes.numel();
             Matrix< DDRMat > tSideXHat( tNumSideNodes, 3 );
             Matrix< DDRMat > tSideXiHat( tNumSideNodes, 4 );
             for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
             {
                 tSideXHat.get_row( iNode )  = tXHat.get_row( tSideOrdNodes( iNode ) );
                 tSideXiHat.get_row( iNode ) = tXiHat.get_row( tSideOrdNodes( iNode ) );
             }

             //set the coefficients xHat, tHat
             tSideGeoInterp.set_space_coeff( tSideXHat );
             tSideGeoInterp.set_time_coeff( tTHat );

             // set the param coefficients xiHat, tauHat
             tSideGeoInterp.set_space_param_coeff( tSideXiHat );
             tSideGeoInterp.set_time_param_coeff( tTauHat );

             // init the side surface and the normal
             real tSideSurface = 0;

             // loop over the integration points
             for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
             {
                 // get the treated integration point location in the surface ref space
                 Matrix< DDRMat > tLocalSideIntegPoint = tSideIntegPoints.get_column( iGP );

                 // set the treated integration point location in the surface ref space for the geometry interp
                 tSideGeoInterp.set_space_time( tLocalSideIntegPoint );

                 // evaluate surfDetJ
                 real tSurfDetJ = tSideGeoInterp.det_J();

                 // add contribution to the surface
                 tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );

                 // get normal
                 Matrix< DDRMat > tSideNormal;
                 tSideGeoInterp.get_normal( tSideNormal );

                 // check the normal
                 bool tSideNormalCheck = true;
                 for( uint i = 0; i < 3; i++ )
                 {
                     tSideNormalCheck = tSideNormalCheck && ( std::abs( tSideNormal( i ) - tSideNormalsExact( i, iSide ) ) < tEpsilon );
                 }
                 REQUIRE( tSideNormalCheck );
             }

             // check the surface value
             tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );
         }
         REQUIRE( tSideSurfaceCheck );
    }
}

TEST_CASE( "Side_Geometry_Interpolation : HEX8", "[moris],[fem],[SideGeoInterp_HEX]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // Side surface for comparison
    Matrix< DDRMat > tSideSurfaceExact( 6, 1, 1.0 );

    // Side normals for comparison
    Matrix< DDRMat > tSideNormalsExact( 3, 6 );
    tSideNormalsExact( 0, 0 ) =  0.0;
    tSideNormalsExact( 1, 0 ) = -1.0;
    tSideNormalsExact( 2, 0 ) =  0.0;
    tSideNormalsExact( 0, 1 ) =  1.0;
    tSideNormalsExact( 1, 1 ) =  0.0;
    tSideNormalsExact( 2, 1 ) =  0.0;
    tSideNormalsExact( 0, 2 ) =  0.0;
    tSideNormalsExact( 1, 2 ) =  1.0;
    tSideNormalsExact( 2, 2 ) =  0.0;
    tSideNormalsExact( 0, 3 ) = -1.0;
    tSideNormalsExact( 1, 3 ) =  0.0;
    tSideNormalsExact( 2, 3 ) =  0.0;
    tSideNormalsExact( 0, 4 ) =  0.0;
    tSideNormalsExact( 1, 4 ) =  0.0;
    tSideNormalsExact( 2, 4 ) = -1.0;
    tSideNormalsExact( 0, 5 ) =  0.0;
    tSideNormalsExact( 1, 5 ) =  0.0;
    tSideNormalsExact( 2, 5 ) =  1.0;

    // loop over the interpolation order
    for( uint iOrder = 0; iOrder < 1; iOrder++ )
    {
        // define a TRI of order tSpaceInterpolationOrder
        Matrix< DDRMat > tXHat;
        Matrix< DDRMat > tXiHat;
        mtk::Interpolation_Order tSpaceInterpolationOrder;
        Matrix< DDSMat > tSideNodes;

        switch( iOrder )
        {
            case ( 0 ):
            {
                // define a HEX8 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                tXHat = { { 0.0, 0.0, 0.0 },
                          { 1.0, 0.0, 0.0 },
                          { 1.0, 1.0, 0.0 },
                          { 0.0, 1.0, 0.0 },
                          { 0.0, 0.0, 1.0 },
                          { 1.0, 0.0, 1.0 },
                          { 1.0, 1.0, 1.0 },
                          { 0.0, 1.0, 1.0 }};
                tXiHat = {{ -1.0, -1.0, -1.0 },
                          {  1.0, -1.0, -1.0 },
                          {  1.0,  1.0, -1.0 },
                          { -1.0,  1.0, -1.0 },
                          { -1.0, -1.0,  1.0 },
                          {  1.0, -1.0,  1.0 },
                          {  1.0,  1.0,  1.0 },
                          { -1.0,  1.0,  1.0 }};
                tSideNodes = { { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 },
                               { 0, 4, 7, 3 }, { 0, 3, 2, 1 }, { 4, 5, 6, 7 } };
                break;
            }
            // FIXME for QUADRATIC and CUBIC
        }

        //create a line time element
        Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 } };
        Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 } };

        // side geometry type and space ordinal
        mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::QUAD;

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tSideGeoInterpRule( tSideGeometryType,
                                               mtk::Interpolation_Type::LAGRANGE,
                                               tSpaceInterpolationOrder,
                                               mtk::Interpolation_Type::LAGRANGE,
                                               mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tSideGeoInterp( tSideGeoInterpRule, mtk::CellShape::GENERAL, true );

        // side geometry type and space ordinal
        Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2, 3, 4, 5 };

        // create a side integration
        mtk::Integration_Rule tSideIntegrationRule( tSideGeometryType,
                                               mtk::Integration_Type::GAUSS,
                                               mtk::Integration_Order::QUAD_3x3,
                                               mtk::Integration_Type::GAUSS,
                                               mtk::Integration_Order::BAR_1 );

        // create a side integrator
        mtk::Integrator tSideIntegrator( tSideIntegrationRule );

        //get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tSideIntegPoints;
        tSideIntegrator.get_points( tSideIntegPoints );
        Matrix< DDRMat > tSideIntegWeights;
        tSideIntegrator.get_weights( tSideIntegWeights );

        // booleans for checks
        bool tSideSurfaceCheck = true;

        // loop over the element sides
        for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
        {
            // get treated side ordinal
            moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

            // get the node ids associated to the follower side ordinal
            Matrix< DDSMat > tSideOrdNodes = tSideNodes.get_row( tSideOrdinal );

            // phys coords and param coords in IP param space for the follower side
            uint tNumSideNodes = tSideOrdNodes.numel();
            Matrix< DDRMat > tSideXHat( tNumSideNodes, 3 );
            Matrix< DDRMat > tSideXiHat( tNumSideNodes, 3 );
            for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
            {
                tSideXHat.get_row( iNode )  = tXHat.get_row( tSideOrdNodes( iNode ) );
                tSideXiHat.get_row( iNode ) = tXiHat.get_row( tSideOrdNodes( iNode ) );
            }

            //set the coefficients xHat, tHat
            tSideGeoInterp.set_space_coeff( tSideXHat );
            tSideGeoInterp.set_time_coeff( tTHat );

            // set the param coefficients xiHat, tauHat
            tSideGeoInterp.set_space_param_coeff( tSideXiHat );
            tSideGeoInterp.set_time_param_coeff( tTauHat );

            // init the side surface
            real tSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tLocalSideIntegPoint = tSideIntegPoints.get_column( iGP );

                // set the treated integration point location in the surface ref space for the geometry interp
                tSideGeoInterp.set_space_time( tLocalSideIntegPoint );

                // evaluate surfDetJ
                real tSurfDetJ = tSideGeoInterp.det_J();

                // add contribution to the surface
                tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );

                // get normal
                Matrix< DDRMat > tSideNormal;
                tSideGeoInterp.get_normal( tSideNormal );

                // check the normal
                bool tSideNormalCheck = true;
                for( uint i = 0; i < 3; i++ )
                {
                    tSideNormalCheck = tSideNormalCheck && ( std::abs( tSideNormal( i ) - tSideNormalsExact( i, iSide ) ) < tEpsilon );
                }
                REQUIRE( tSideNormalCheck );
            }

            // check the surface value
            tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );

        }
        REQUIRE( tSideSurfaceCheck );
    }
}
