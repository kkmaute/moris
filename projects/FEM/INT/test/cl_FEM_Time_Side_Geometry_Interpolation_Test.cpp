/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Time_Side_Geometry_Interpolation_Test.cpp
 *
 */

#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr
#include "cl_MTK_Integrator.hpp" //MTK/sr

using namespace moris;
using namespace fem;

TEST_CASE("Time_Side_Geometry_Interpolation : QUAD4 - QUAD9 - QUAD16", "[moris],[fem],[TimeSideGeoInterp_QUAD]")
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // loop over the interpolation order
    for( uint iOrder = 0; iOrder < 3; iOrder++ )
    {
        // define a QUAD of order tSpaceInterpolationOrder
        Matrix< DDRMat > tXHat;
        Matrix< DDRMat > tXiHat;
        mtk::Interpolation_Order tSpaceInterpolationOrder;

        // loop over the interpolation orders
        switch( iOrder )
        {
            case ( 0 ):
            {
                // define a QUAD4 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                tXHat = {{ 0.0, 0.0 }, { 3.0, 1.25 }, { 4.5, 4.0 }, { 1.0, 3.25} };
                // define a QUAD4 in the parametric space, i.e. space param coordinates xiHat
                tXiHat = {{ -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 }};
                break;
            }

            case ( 1 ):
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
                break;
            }
            case ( 2 ):
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
                break;
            }
        }

        //create a line time element
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};
        Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 }};

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tSideGeoInterpRule( mtk::Geometry_Type::QUAD,
                                               mtk::Interpolation_Type::LAGRANGE,
                                               tSpaceInterpolationOrder,
                                               mtk::Geometry_Type::POINT,
                                               mtk::Interpolation_Type::CONSTANT,
                                               mtk::Interpolation_Order::UNDEFINED );

        // create a space and time geometry interpolator
        Geometry_Interpolator tSideGeoInterp( tSideGeoInterpRule, mtk::CellShape::GENERAL, false, true );

        // time side geometry type and space ordinal------------------------------------
        Cell< moris::moris_index > tListOfTimeOrdinals = { 0, 1 };

        // create a side integration
        mtk::Integration_Rule tTimeSideIntegRule( mtk::Geometry_Type::QUAD,
                                             mtk::Integration_Type::GAUSS,
                                             mtk::Integration_Order::QUAD_3x3,
                                             mtk::Geometry_Type::POINT,
                                             mtk::Integration_Type::GAUSS,
                                             mtk::Integration_Order::UNDEFINED );

        // create a side integrator
        mtk::Integrator tTimeSideIntegrator( tTimeSideIntegRule );

        //get number of integration points, integration points and weights
        uint             tNumOfTimeIntegPoints = tTimeSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tTimeSideIntegPoints;
        tTimeSideIntegrator.get_points( tTimeSideIntegPoints );
        Matrix< DDRMat > tTimeSideIntegWeights;
        tTimeSideIntegrator.get_weights( tTimeSideIntegWeights );

        // boolean for surface check
        bool tTimeSurfaceCheck = true;

        // loop over the time side
        for ( uint iTimeSide = 0; iTimeSide < tListOfTimeOrdinals.size(); iTimeSide++ )
        {
            // get the time ordinal
            moris_index tTimeSideOrdinal = tListOfTimeOrdinals( iTimeSide );

            //set the coefficients xHat, tHat
           tSideGeoInterp.set_space_coeff( tXHat );
           tSideGeoInterp.set_time_coeff( {{ tTHat( tTimeSideOrdinal ) }} );

           // set the param coefficients xiHat, tauHat
           tSideGeoInterp.set_space_param_coeff( tXiHat );
           tSideGeoInterp.set_time_param_coeff( {{ tTauHat( tTimeSideOrdinal )}} );

            // init the side surface
            real tTimeSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfTimeIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tTimeSideIntegPointI = tTimeSideIntegPoints.get_column( iGP );

                // set the treated integration point location in the surface ref space for the geometry interp
                tSideGeoInterp.set_space_time( tTimeSideIntegPointI );

                // evaluate timesurfDetJ
                real tSurfDetJ = tSideGeoInterp.det_J();

                // add contribution to the surface
                tTimeSideSurface = tTimeSideSurface + tSurfDetJ * tTimeSideIntegWeights( iGP );
            }
            // check the surface value
            tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 8.5 ) < tEpsilon );
        }
        REQUIRE( tTimeSurfaceCheck );
    }
}

//------------------------------------------------------------------------------
TEST_CASE( "Time_Side_Geometry_Interpolation : TRI3 - TRI6 - TRI10", "[moris],[fem],[TimeSideGeoInterp_TRI]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // loop over the interpolation order
    for( uint iOrder = 0; iOrder < 3; iOrder++ )
    {
        // define a TRI of order tSpaceInterpolationOrder
        Matrix< DDRMat > tXHat;
        Matrix< DDRMat > tXiHat;
        mtk::Interpolation_Order tSpaceInterpolationOrder;

        switch( iOrder )
        {
            case ( 0 ):
            {
                // define a TRI4 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                tXHat = {{ 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }};
                tXiHat = {{ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 }};
                break;
            }
            case ( 1 ):
            {
                // define a TRI6 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;
                real t12 = 1.0/2.0;
                tXHat = {{ 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 },
                         { t12, 0.0 }, { t12, t12 }, { 0.0, t12 }};
                tXiHat = {{ 1.000000000000000, 0.000000000000000, 0.000000000000000,
                            0.500000000000000, 0.000000000000000, 0.500000000000000 },
                          { 0.000000000000000, 1.000000000000000, 0.000000000000000,
                            0.500000000000000, 0.500000000000000, 0.000000000000000 },
                          { 0.000000000000000, 0.000000000000000, 1.000000000000000,
                            0.000000000000000, 0.500000000000000, 0.500000000000000 }};
                tXiHat = trans( tXiHat );
                break;
            }
            case ( 2 ):
            {
                // define a TRI10 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::CUBIC;
                real t13 = 1.0/3.0;
                real t23 = 2.0/3.0;
                tXHat = {{ 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 },
                         { t13, 0.0 }, { t23, 0.0 },
                         { t23, t13 }, { t13, t23 },
                         { 0.0, t23 }, { 0.0, t13 },
                         { t13, t13 }};
                tXiHat = {{ 1.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t23, t13 },
                          { 0.0, 1.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, t13 },
                          { 0.0, 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, t13 }};
                tXiHat = trans( tXiHat );
                break;
            }
        }

        // define a LINE2 in the physical space, i.e. time coordinates tHat
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 2.0 }};
        Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 } };

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tSideGeoInterpRule( mtk::Geometry_Type::TRI,
                                               mtk::Interpolation_Type::LAGRANGE,
                                               tSpaceInterpolationOrder,
                                               mtk::Geometry_Type::POINT,
                                               mtk::Interpolation_Type::CONSTANT,
                                               mtk::Interpolation_Order::UNDEFINED );

       // create a space and time geometry interpolator
       Geometry_Interpolator tSideGeoInterp( tSideGeoInterpRule, mtk::CellShape::GENERAL, false, true );

       // time side geometry type and space ordinal
       Cell< moris::moris_index > tListOfTimeOrdinals = { 0, 1 };

       // create a side integration
       mtk::Integration_Rule tTimeSideIntegRule( mtk::Geometry_Type::TRI,
                                            mtk::Integration_Type::GAUSS,
                                            mtk::Integration_Order::TRI_6,
                                            mtk::Geometry_Type::POINT,
                                            mtk::Integration_Type::GAUSS,
                                            mtk::Integration_Order::UNDEFINED );

       // create a side integrator
       mtk::Integrator tTimeSideIntegrator( tTimeSideIntegRule );

       //get number of integration points, integration points and weights
       uint             tNumOfTimeIntegPoints = tTimeSideIntegrator.get_number_of_points();
       Matrix< DDRMat > tTimeSideIntegPoints;
       tTimeSideIntegrator.get_points( tTimeSideIntegPoints );
       Matrix< DDRMat > tTimeSideIntegWeights;
       tTimeSideIntegrator.get_weights( tTimeSideIntegWeights );

       // boolean for surface check
       bool tTimeSurfaceCheck = true;

       // loop over the time side
       for ( uint iTimeSide = 0; iTimeSide < tListOfTimeOrdinals.size(); iTimeSide++ )
       {
           // get the time ordinal
           moris_index tTimeSideOrdinal = tListOfTimeOrdinals( iTimeSide );

           //set the coefficients xHat, tHat
          tSideGeoInterp.set_space_coeff( tXHat );
          tSideGeoInterp.set_time_coeff( {{ tTHat( tTimeSideOrdinal ) }} );

          // set the param coefficients xiHat, tauHat
          tSideGeoInterp.set_space_param_coeff( tXiHat );
          tSideGeoInterp.set_time_param_coeff( {{ tTauHat( tTimeSideOrdinal )}} );

           // init the side surface
           real tTimeSideSurface = 0;

           // loop over the integration points
           for ( uint iGP = 0; iGP < tNumOfTimeIntegPoints; iGP++ )
           {
               // get the treated integration point location in the surface ref space
               Matrix< DDRMat > tTimeSideIntegPointI = tTimeSideIntegPoints.get_column( iGP );

               // set the treated integration point location in the surface ref space for the geometry interp
               tSideGeoInterp.set_space_time( tTimeSideIntegPointI );

               // evaluate timesurfDetJ
               real tSurfDetJ = tSideGeoInterp.det_J();

               // add contribution to the surface
               tTimeSideSurface = tTimeSideSurface + tSurfDetJ * tTimeSideIntegWeights( iGP );
           }
           // check the surface value
           tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 0.5 ) < tEpsilon );
           REQUIRE( tTimeSurfaceCheck );
       }
    }
}

//------------------------------------------------------------------------------
TEST_CASE( "Time_Side_Geometry_Interpolation :  TET4 - TET10 - TET20", "[moris],[fem],[TimeSideGeoInterp_TET]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

     // switch on the TET interpolation order
     for ( uint iOrder = 0; iOrder < 3; iOrder++ )
     {
         // define a TET of order tSpaceInterpolationOrder
         Matrix< DDRMat > tXHat;
         Matrix< DDRMat > tXiHat;
         mtk::Interpolation_Order tSpaceInterpolationOrder;

         switch( iOrder )
         {
             case ( 0 ):
             {
                 // define a TET4 in the physical space, i.e. space coordinates xHat
                 tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                 tXHat = {{ 0.0,  0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 1.0,  0.0, 0.0 }, { 0.0,  0.0, 1.0 }};
                 tXiHat = {{ 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000 },
                           { 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000 },
                           { 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000 },
                           { 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000 }};
                 tXiHat = trans( tXiHat );
                 break;
             }
             case ( 1 ):
             {
                 // define a TET10 in the physical space, i.e. space coordinates xHat
                 tSpaceInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;
                 real t12 = 1.0/2.0;
                 tXHat = {{ 0.0,  0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 1.0,  0.0, 0.0 }, { 0.0,  0.0, 1.0 },
                          { 0.0, -t12, 0.0 }, { t12, -t12, 0.0 },
                          { t12,  0.0, 0.0 }, { 0.0,  0.0, t12 },
                          { 0.0, -t12, t12 }, { t12,  0.0, t12 } };
                 tXiHat ={{ 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000,
                            0.000000000000000, 0.500000000000000, 0.500000000000000, 0.000000000000000, 0.000000000000000 },
                          { 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000,
                            0.500000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000, 0.000000000000000 },
                          { 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000, 0.000000000000000,
                            0.500000000000000, 0.500000000000000, 0.000000000000000, 0.000000000000000, 0.500000000000000 },
                          { 0.000000000000000, 0.000000000000000, 0.000000000000000, 1.000000000000000, 0.000000000000000,
                            0.000000000000000, 0.000000000000000, 0.500000000000000, 0.500000000000000, 0.500000000000000 }};
                 tXiHat = trans( tXiHat );
                 break;
             }
             case ( 2 ):
             {
                 // define a TET20 in the physical space, i.e. space coordinates xHat
                 tSpaceInterpolationOrder = mtk::Interpolation_Order::CUBIC;
                 real t13 = 1.0/3.0;
                 real t23 = 2.0/3.0;
                 tXHat = {{ 0.0,  0.0, 0.0 }, { 0.0, -1.0, 0.0 }, { 1.0,  0.0, 0.0 }, { 0.0,  0.0, 1.0 },
                          { 0.0, -t13, 0.0 }, { 0.0, -t23, 0.0 },
                          { t13, -t23, 0.0 }, { t23, -t13, 0.0 },
                          { t13,  0.0, 0.0 }, { t23,  0.0, 0.0 },
                          { 0.0,  0.0, t13 }, { 0.0,  0.0, t23 },
                          { 0.0, -t23, t13 }, { 0.0, -t13, t23 },
                          { t23,  0.0, t13 }, { t13,  0.0, t23 },
                          { t13, -t13, 0.0 }, { 0.0, -t13, t13 },
                          { t13, -t13, t13 }, { t13,  0.0, t13 } };
                 tXiHat = {{ 1.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t23, t13, t23, t13, 0.0, 0.0, 0.0, 0.0, t13, t13, 0.0, t13 },
                            { 0.0, 1.0, 0.0, 0.0, t13, t23, t23, t13, 0.0, 0.0, 0.0, 0.0, t23, t13, 0.0, 0.0, t13, t13, t13, 0.0 },
                            { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, 0.0, 0.0, 0.0, 0.0, t23, t13, t13, 0.0, t13, t13 },
                            { 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, t13, t23, t13, t23, t13, t23, 0.0, t13, t13, t13 }};
                 tXiHat = trans( tXiHat );
                 break;
             }
         }

         // define a LINE2 in the physical space, i.e. time coordinates tHat
         Matrix< DDRMat > tTHat = {{ 0.0 }, { 2.0 }};
         Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 } };

         // create a space and time geometry interpolation rule
         mtk::Interpolation_Rule tSideGeoInterpRule( mtk::Geometry_Type::TET,
                                                mtk::Interpolation_Type::LAGRANGE,
                                                tSpaceInterpolationOrder,
                                                mtk::Geometry_Type::POINT,
                                                mtk::Interpolation_Type::CONSTANT,
                                                mtk::Interpolation_Order::UNDEFINED );

         // create a space and time geometry interpolator
         Geometry_Interpolator tSideGeoInterp( tSideGeoInterpRule, mtk::CellShape::GENERAL, false, true );

         // time side geometry type and space ordinal
         Cell< moris::moris_index > tListOfTimeOrdinals = { 0, 1 };

         // create a side integration
         mtk::Integration_Rule tTimeSideIntegRule( mtk::Geometry_Type::TET,
                                              mtk::Integration_Type::GAUSS,
                                              mtk::Integration_Order::TET_5,
                                              mtk::Geometry_Type::POINT,
                                              mtk::Integration_Type::GAUSS,
                                              mtk::Integration_Order::UNDEFINED );

         // create a side integrator
         mtk::Integrator tTimeSideIntegrator( tTimeSideIntegRule );

         //get number of integration points, integration points and weights
         uint             tNumOfTimeIntegPoints = tTimeSideIntegrator.get_number_of_points();
         Matrix< DDRMat > tTimeSideIntegPoints;
         tTimeSideIntegrator.get_points( tTimeSideIntegPoints );
         Matrix< DDRMat > tTimeSideIntegWeights;
         tTimeSideIntegrator.get_weights( tTimeSideIntegWeights );

         // loop over the time side
         for ( uint iTimeSide = 0; iTimeSide < tListOfTimeOrdinals.size(); iTimeSide++ )
         {
             // get the time ordinal
             moris_index tTimeSideOrdinal = tListOfTimeOrdinals( iTimeSide );

             //set the coefficients xHat, tHat
            tSideGeoInterp.set_space_coeff( tXHat );
            tSideGeoInterp.set_time_coeff( {{ tTHat( tTimeSideOrdinal ) }} );

            // set the param coefficients xiHat, tauHat
            tSideGeoInterp.set_space_param_coeff( tXiHat );
            tSideGeoInterp.set_time_param_coeff( {{ tTauHat( tTimeSideOrdinal )}} );

             // init the side surface
             real tTimeSideSurface = 0;

             // boolean for surface check
             bool tTimeSurfaceCheck = true;

             // loop over the integration points
             for ( uint iGP = 0; iGP < tNumOfTimeIntegPoints; iGP++ )
             {
                 // get the treated integration point location in the surface ref space
                 Matrix< DDRMat > tTimeSideIntegPointI = tTimeSideIntegPoints.get_column( iGP );

                 // set the treated integration point location in the surface ref space for the geometry interp
                 tSideGeoInterp.set_space_time( tTimeSideIntegPointI );

                 // evaluate timesurfDetJ
                 real tSurfDetJ = tSideGeoInterp.det_J();

                 // add contribution to the surface
                 tTimeSideSurface = tTimeSideSurface + tSurfDetJ * tTimeSideIntegWeights( iGP );
             }
             // check the surface value
             tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 1.0/6.0 ) < tEpsilon );
             REQUIRE( tTimeSurfaceCheck );
         }
    }
}

TEST_CASE( "Time_Side_Geometry_Interpolation : HEX8", "[moris],[fem],[TimeSideGeoInterp_HEX]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // loop over the interpolation order
    for( uint iOrder = 0; iOrder < 1; iOrder++ )
    {
        // define a TRI of order tSpaceInterpolationOrder
        Matrix< DDRMat > tXHat;
        Matrix< DDRMat > tXiHat;
        mtk::Interpolation_Order tSpaceInterpolationOrder;

        switch( iOrder )
        {
            case ( 0 ):
            {
                // define a HEX8 in the physical space, i.e. space coordinates xHat
                tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                tXHat = { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
                          { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 }};
                tXiHat = {{ -1.0, -1.0, -1.0 }, {  1.0, -1.0, -1.0 }, {  1.0,  1.0, -1.0 }, { -1.0,  1.0, -1.0 },
                          { -1.0, -1.0,  1.0 }, {  1.0, -1.0,  1.0 }, {  1.0,  1.0,  1.0 }, { -1.0,  1.0,  1.0 }};
                break;
            }
            // FIXME for QUADRATIC and CUBIC
        }

        //create a line time element
        Matrix< DDRMat > tTHat   = {{  0.0 }, { 1.0 } };
        Matrix< DDRMat > tTauHat = {{ -1.0 }, { 1.0 } };

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tSideGeoInterpRule( mtk::Geometry_Type::HEX,
                                               mtk::Interpolation_Type::LAGRANGE,
                                               tSpaceInterpolationOrder,
                                               mtk::Geometry_Type::POINT,
                                               mtk::Interpolation_Type::CONSTANT,
                                               mtk::Interpolation_Order::UNDEFINED );

        // create a space and time geometry interpolator
        Geometry_Interpolator tSideGeoInterp( tSideGeoInterpRule, mtk::CellShape::GENERAL, false, true );

        // time side geometry type and space ordinal
        Cell< moris::moris_index > tListOfTimeOrdinals = { 0, 1 };

        // create a side integration
        mtk::Integration_Rule tTimeSideIntegRule( mtk::Geometry_Type::HEX,
                                             mtk::Integration_Type::GAUSS,
                                             mtk::Integration_Order::HEX_3x3x3,
                                             mtk::Geometry_Type::POINT,
                                             mtk::Integration_Type::GAUSS,
                                             mtk::Integration_Order::UNDEFINED );

        // create a side integrator
        mtk::Integrator tTimeSideIntegrator( tTimeSideIntegRule );

        //get number of integration points, integration points and weights
        uint             tNumOfTimeIntegPoints = tTimeSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tTimeSideIntegPoints;
        tTimeSideIntegrator.get_points( tTimeSideIntegPoints );
        Matrix< DDRMat > tTimeSideIntegWeights;
        tTimeSideIntegrator.get_weights( tTimeSideIntegWeights );

        // boolean for surface check
        bool tTimeSurfaceCheck = true;

        // loop over the time side
        for ( uint iTimeSide = 0; iTimeSide < tListOfTimeOrdinals.size(); iTimeSide++ )
        {
            // get the time ordinal
            moris_index tTimeSideOrdinal = tListOfTimeOrdinals( iTimeSide );

            //set the coefficients xHat, tHat
           tSideGeoInterp.set_space_coeff( tXHat );
           tSideGeoInterp.set_time_coeff( {{ tTHat( tTimeSideOrdinal ) }} );

           // set the param coefficients xiHat, tauHat
           tSideGeoInterp.set_space_param_coeff( tXiHat );
           tSideGeoInterp.set_time_param_coeff( {{ tTauHat( tTimeSideOrdinal )}} );

            // init the side surface
            real tTimeSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfTimeIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tTimeSideIntegPointI = tTimeSideIntegPoints.get_column( iGP );

                // set the treated integration point location in the surface ref space for the geometry interp
                tSideGeoInterp.set_space_time( tTimeSideIntegPointI );

                // evaluate timesurfDetJ
                real tSurfDetJ = tSideGeoInterp.det_J();

                // add contribution to the surface
                tTimeSideSurface = tTimeSideSurface + tSurfDetJ * tTimeSideIntegWeights( iGP );
            }
            // check the surface value
            tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 1.0 ) < tEpsilon );
            REQUIRE( tTimeSurfaceCheck );
        }
    }
}
