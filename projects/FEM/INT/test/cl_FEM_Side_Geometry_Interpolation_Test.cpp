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
        moris::moris_index tSideOrdinal = 3;

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
        //print( tSideIntegWeights, "tSideIntegWeights" );

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
        //std::cout<<tSideSurface<<std::endl;

        // time side geometry type and space ordinal------------------------------------
        //------------------------------------------------------------------------------

        moris::moris_index tTimeSideOrdinal = 0;

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
        //print( tSideIntegWeights, "tSideIntegWeights" );

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
        //std::cout<<tTimeSideSurface<<std::endl;

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
    	Matrix< DDRMat > tXHat = { { 0.0, 0.0,  0.0 },
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
        moris::moris_index tSideOrdinal = 0;

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
        //std::cout<<tSideSurface<<std::endl;

        // time side geometry type and space ordinal------------------------------------
        //------------------------------------------------------------------------------

        moris::moris_index tTimeSideOrdinal = 0;

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
            std::cout<<tSurfDetJ<<std::endl;
        }
        std::cout<<tTimeSideSurface<<std::endl;


    }
//------------------------------------------------------------------------------
}


