#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr
#include "cl_FEM_Integrator.hpp" //FEM/INT/sr

using namespace moris;
using namespace fem;

TEST_CASE( "Side_Geometry_Interpolation", "[moris],[fem],[SideGeoInterp]" )
{
    // define an epsilon environment
    double tEpsilon = 1E-12;

    SECTION( "Geometry Interpolator : 2D space time" )
    {
        // space and time geometry interpolator
        //------------------------------------------------------------------------------
        // create a quad4 space element
        Matrix< DDRMat > tXHat = {{ 0.0, 0.0 },
                                  { 3.0, 1.25 },
                                  { 4.5, 4.0 },
                                  { 1.0, 3.25} };

        //create a line time element
        Matrix< DDRMat > tTHat = {{ 0.0 },
                                  { 1.0 }};

        // get side surface for comparison
        Matrix< IndexMat > tSideListOfNodes = {{ 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }};
        Matrix< DDRMat >   tSideSurfaceExact( 4, 1, 0.0 );
        for( uint iSide = 0; iSide < 4; iSide++ )
        {
            Matrix< IndexMat > tTreatedNodes = tSideListOfNodes.get_row( iSide );
            tSideSurfaceExact( iSide ) = norm( tXHat.get_row( tTreatedNodes( 0 ) ) - tXHat.get_row( tTreatedNodes( 1 ) ) )
                                       * ( tTHat( 1 ) - tTHat( 0 ) );
        }

        // create a space and time geometry interpolation rule
        Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

        //set the coefficients xHat, tHat
        tGeomInterpolator.set_coeff( tXHat, tTHat );

        // side geometry type and space ordinal-----------------------------------------
        //------------------------------------------------------------------------------

        mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::LINE;
        Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2, 3};

        // create a side integration
        Integration_Rule tSideIntegrationRule( tSideGeometryType,
                                               Integration_Type::GAUSS,
                                               Integration_Order::BAR_3,
                                               Integration_Type::GAUSS,
                                               Integration_Order::BAR_1 );

        // create a side integrator
        Integrator tSideIntegrator( tSideIntegrationRule );

        //get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tSideIntegPoints  = tSideIntegrator.get_points();
        Matrix< DDRMat > tSideIntegWeights = tSideIntegrator.get_weights();


        bool tSideSurfaceCheck = true;

        for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
        {

            // get treated side ordinal
            moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

            // init the side surface
            real tSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tSideIntegPointI = tSideIntegPoints.get_column( iGP );

                // get the treated integration point location in the volume ref space
                Matrix< DDRMat > tVolIntegPointI = tGeomInterpolator.surf_val( tSideIntegPointI,
                                                                               tSideOrdinal );

                // evaluate surfDetJ and normal
                real tSurfDetJ;
                Matrix< DDRMat > tNormal;
                tGeomInterpolator.surf_det_J( tSurfDetJ, tNormal,
                                              tSideIntegPointI, tSideOrdinal );
                //print( tNormal, "tNormal" );

                // add contribution to the surface
                tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );
            }

            tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );
            //std::cout<<tSideSurface<<std::endl;
        }
        REQUIRE( tSideSurfaceCheck );

        // time side geometry type and space ordinal------------------------------------
        //------------------------------------------------------------------------------

        Matrix< IndexMat > tListOfTimeOrdinals = {{ 0 }, { 1 }};

        // create a side integration
        Integration_Rule tTimeSideIntegrationRule( mtk::Geometry_Type::QUAD,
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::QUAD_3x3,
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_1 );

        // create a side integrator
        Integrator tTimeSideIntegrator( tTimeSideIntegrationRule );

        //get number of integration points, integration points and weights
        uint             tNumOfTimeIntegPoints = tTimeSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tTimeSideIntegPoints  = tTimeSideIntegrator.get_points();
        Matrix< DDRMat > tTimeSideIntegWeights = tTimeSideIntegrator.get_weights();

        // init surface check
        bool tTimeSurfaceCheck = true;

        // loop over the time side
        for ( uint iTimeSide = 0; iTimeSide < tListOfTimeOrdinals.numel(); iTimeSide++ )
        {
            // get the time ordinal
            moris_index tTimeSideOrdinal = tListOfTimeOrdinals( iTimeSide );

            // init the side surface
            real tTimeSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfTimeIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tTimeSideIntegPointI = tTimeSideIntegPoints.get_column( iGP );

                // get the treated integration point location in the volume ref space
                Matrix< DDRMat > tTimeVolIntegPointI = tGeomInterpolator.time_surf_val( tTimeSideIntegPointI,
                                                                                        tTimeSideOrdinal );

                // evaluate surfDetJ and normal
                real tSurfDetJ;
                tGeomInterpolator.time_surf_det_J( tSurfDetJ,
                                                   tTimeSideIntegPointI,
                                                   tTimeSideOrdinal );

                // add contribution to the surface
                tTimeSideSurface = tTimeSideSurface + tSurfDetJ * tTimeSideIntegWeights( iGP );
            }
            // check the surface value
            tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 8.5 ) < tEpsilon );
            //std::cout<<tTimeSideSurface<<std::endl;
        }
        REQUIRE( tTimeSurfaceCheck );

    }

    SECTION( "Geometry Interpolator : 3D space time" )
    {
        // space and time geometry interpolator
        //------------------------------------------------------------------------------
        // create a quad4 space element
//        Matrix< DDRMat > tXHat = { { 0.0, 0.0,  0.0 },
//                                   { 3.0, 1.25, 0.0 },
//                                   { 4.5, 4.0,  0.0 },
//                                   { 1.0, 3.25, 0.0 },
//                                   { 1.0, 1.0,  1.0 },
//                                   { 4.0, 2.25, 1.5 },
//                                   { 5.5, 5.0,  2.0 },
//                                   { 2.0, 4.25, 3.0 }};

    	Matrix< DDRMat > tXHat = { { 0.0, 0.0, 0.0 },
    	                           { 1.0, 0.0, 0.0 },
    	                           { 1.0, 1.0, 0.0 },
    	                           { 0.0, 1.0, 0.0 },
    	                           { 0.0, 0.0, 1.0 },
    	                           { 1.0, 0.0, 1.0 },
    	                           { 1.0, 1.0, 1.0 },
    	                           { 0.0, 1.0, 1.0 }};

        //create a line time element
        Matrix< DDRMat > tTHat = { { 0.0 },
                                   { 1.0 } };

        // get side surface for comparison
        Matrix< DDRMat >   tSideSurfaceExact( 6, 1, 1.0 );

        // create a space and time geometry interpolation rule
        Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

        //set the coefficients xHat, tHat
        tGeomInterpolator.set_coeff( tXHat, tTHat );

        // side geometry type and space ordinal
        mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::QUAD;
        Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2, 3, 4, 5 };

        // create a side integration
        Integration_Rule tSideIntegrationRule( tSideGeometryType,
                                               Integration_Type::GAUSS,
                                               Integration_Order::QUAD_3x3,
                                               Integration_Type::GAUSS,
                                               Integration_Order::BAR_1 );

        // create a side integrator
        Integrator tSideIntegrator( tSideIntegrationRule );

        //get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tSideIntegPoints  = tSideIntegrator.get_points();
        Matrix< DDRMat > tSideIntegWeights = tSideIntegrator.get_weights();

        bool tSideSurfaceCheck = true;

        for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
        {

            // get treated side ordinal
            moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

            // init the side surface
            real tSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tSideIntegPointI = tSideIntegPoints.get_column( iGP );

                // get the treated integration point location in the volume ref space
                Matrix< DDRMat > tVolIntegPointI = tGeomInterpolator.surf_val( tSideIntegPointI,
                                                                               tSideOrdinal );

                // evaluate surfDetJ and normal
                real tSurfDetJ;
                Matrix< DDRMat > tNormal;
                tGeomInterpolator.surf_det_J( tSurfDetJ, tNormal,
                                              tSideIntegPointI, tSideOrdinal );
                //print( tNormal, "tNormal" );

                // add contribution to the surface
                tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );
            }

            tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );
            //std::cout<<tSideSurface<<std::endl;
        }
        REQUIRE( tSideSurfaceCheck );

        // time side geometry type and space ordinal------------------------------------
        //------------------------------------------------------------------------------

        Matrix< IndexMat > tListOfTimeOrdinals = {{ 0 }, { 1 }};

        // create a side integration
        Integration_Rule tTimeSideIntegrationRule( mtk::Geometry_Type::HEX,
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::HEX_3x3x3,
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_1 );

        // create a side integrator
        Integrator tTimeSideIntegrator( tTimeSideIntegrationRule );

        //get number of integration points, integration points and weights
        uint             tNumOfTimeIntegPoints = tTimeSideIntegrator.get_number_of_points();
        Matrix< DDRMat > tTimeSideIntegPoints  = tTimeSideIntegrator.get_points();
        Matrix< DDRMat > tTimeSideIntegWeights = tTimeSideIntegrator.get_weights();
        //print( tSideIntegWeights, "tSideIntegWeights" );

        // init surface check
        bool tTimeSurfaceCheck = true;

        // loop over the time side
        for ( uint iTimeSide = 0; iTimeSide < tListOfTimeOrdinals.numel(); iTimeSide++ )
        {
            // get the time ordinal
            moris_index tTimeSideOrdinal = tListOfTimeOrdinals( iTimeSide );

            // init the side surface
            real tTimeSideSurface = 0;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfTimeIntegPoints; iGP++ )
            {
                // get the treated integration point location in the surface ref space
                Matrix< DDRMat > tTimeSideIntegPointI = tTimeSideIntegPoints.get_column( iGP );

                // get the treated integration point location in the volume ref space
                Matrix< DDRMat > tTimeVolIntegPointI = tGeomInterpolator.time_surf_val( tTimeSideIntegPointI,
                                                                                        tTimeSideOrdinal );

                // evaluate surfDetJ and normal
                real tSurfDetJ;
                tGeomInterpolator.time_surf_det_J( tSurfDetJ,
                                                   tTimeSideIntegPointI,
                                                   tTimeSideOrdinal );

                // add contribution to the surface
                tTimeSideSurface = tTimeSideSurface + tSurfDetJ * tTimeSideIntegWeights( iGP );
            }
            // check the surface value
            tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 1.0 ) < tEpsilon );
            //std::cout<<tTimeSideSurface<<std::endl;
        }
        REQUIRE( tTimeSurfaceCheck );

    }
//------------------------------------------------------------------------------
}


