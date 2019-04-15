#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/sr
#include "cl_FEM_Integrator.hpp" //FEM/INT/sr

using namespace moris;
using namespace fem;

TEST_CASE( "Side_Geometry_Interpolation", "[moris],[fem],[SideGeoInterp]" )
{
    // define an epsilon environment
    double tEpsilon = 1E-12;

    SECTION( "Side Geometry Interpolation : QUAD4 - QUAD9 - QUAD16 " )
    {
        // get side surface for comparison
        Matrix< DDRMat > tSideSurfaceExact( 4, 1 );
        tSideSurfaceExact( 0 ) = 3.250000000000000;
        tSideSurfaceExact( 1 ) = 3.132491021535417;
        tSideSurfaceExact( 2 ) = 3.579455265819088;
        tSideSurfaceExact( 3 ) = 3.400367627183861;

        // get side normals for comparison
        Matrix< DDRMat > tSideNormalsExact( 2, 4, 0.0 );
        tSideNormalsExact( 0, 0 ) =  0.3846153846153846;
        tSideNormalsExact( 1, 0 ) = -0.9230769230769231;
        tSideNormalsExact( 0, 1 ) =  0.8778955729143844;
        tSideNormalsExact( 1, 1 ) = -0.4788521306805733;
        tSideNormalsExact( 0, 2 ) = -0.2095290887308734;
        tSideNormalsExact( 1, 2 ) =  0.9778024140774094;
        tSideNormalsExact( 0, 3 ) = -0.9557790087219500;
        tSideNormalsExact( 1, 3 ) =  0.2940858488375231;

        // loop over the interpolation order
        for( uint iOrder = 0; iOrder < 3; iOrder++ )
        {
            // define a TRI of order tSpaceInterpolationOrder
            Matrix< DDRMat > tXHat;
            mtk::Interpolation_Order tSpaceInterpolationOrder;

            // loop over the interpolation orders
            switch( iOrder )
            {
                case ( 0 ):
                {
                    // define a QUAD4 in the physical space, i.e. space coordinates xHat
                    tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                    tXHat = {{ 0.0, 0.0 },
                             { 3.0, 1.25 },
                             { 4.5, 4.0 },
                             { 1.0, 3.25} };
                    break;
                }

                case ( 1 ):
                {
                    // define a QUAD9 in the physical space, i.e. space coordinates xHat
                    tSpaceInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;
                    tXHat = {{ 0.0, 0.0 },
                             { 3.0, 1.25 },
                             { 4.5, 4.0 },
                             { 1.0, 3.25},
                             { 1.5, 0.625 },
                             { 3.75, 2.625 },
                             { 2.75, 3.625 },
                             { 0.5, 1.625 },
                             { 2.125, 2.125 }};
                    break;
                }
                case ( 2 ):
                {
                    // define a QUAD16 in the physical space, i.e. space coordinates xHat
                    tSpaceInterpolationOrder = mtk::Interpolation_Order::CUBIC;
                    tXHat = {{0.000000000000000e+00,  0.000000000000000e+00},
                             {3.000000000000000e+00,  1.250000000000000e+00},
                             {4.500000000000000e+00,  4.000000000000000e+00},
                             {1.000000000000000e+00,  3.250000000000000e+00},
                             {1.000000000000000,  4.166666666666667e-01},
                             {2.000000000000000,  8.333333333333334e-01},
                             {3.500000000000000,  2.166666666666667},
                             {4.000000000000000,  3.083333333333333},
                             {3.333333333333333,  3.750000000000000},
                             {2.166666666666667,  3.500000000000000},
                             {6.666666666666666e-01,  2.166666666666667},
                             {3.333333333333333e-01,  1.083333333333333},
                             {1.388888888888889,  1.444444444444445},
                             {2.444444444444445,  1.805555555555556},
                             {2.888888888888889,  2.777777777777778},
                             {1.777777777777778,  2.472222222222222}};
                    break;
                }
            }

            //create a line time element
            Matrix< DDRMat > tTHat = {{ 0.0 },
                                      { 1.0 }};

            // create a space and time geometry interpolation rule
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                                Interpolation_Type::LAGRANGE,
                                                tSpaceInterpolationOrder,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

            // create a space and time geometry interpolator
            Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

            //set the coefficients xHat, tHat
            tGeomInterpolator.set_coeff( tXHat, tTHat );

            // side geometry type and space ordinal
            mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::LINE;
            Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2, 3 };

            // create a side integration
            Integration_Rule tSideIntegrationRule( tSideGeometryType,
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_5,
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_1 );

            // create a side integrator
            Integrator tSideIntegrator( tSideIntegrationRule );

            //get number of integration points, integration points and weights
            uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
            Matrix< DDRMat > tSideIntegPoints  = tSideIntegrator.get_points();
            Matrix< DDRMat > tSideIntegWeights = tSideIntegrator.get_weights();

            // booleans for checks
            bool tSideSurfaceCheck = true;
            bool tSideNormalCheck  = true;

            // loop over the element sides
            for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
            {
                // get treated side ordinal
                moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

                // init the side surface
                real tSideSurface = 0;
                Matrix< DDRMat > tNormal;

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
                    tGeomInterpolator.surf_det_J( tSurfDetJ, tNormal,
                                                  tSideIntegPointI, tSideOrdinal );

                    // add contribution to the surface
                    tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );
                }

                // check the surface value
                tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );

                // check the normals
                for( uint iNorm = 0; iNorm < 2; iNorm++ )
                {
                    tSideNormalCheck = tSideNormalCheck && ( std::abs( tNormal( iNorm ) - tSideNormalsExact( iNorm, iSide ) ) < tEpsilon );
                }
            }
            REQUIRE( tSideSurfaceCheck );
            REQUIRE( tSideNormalCheck );

            // time side geometry type and space ordinal------------------------------------
            //------------------------------------------------------------------------------

            Matrix< IndexMat > tListOfTimeOrdinals = { 0, 1 };

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

            // boolean for surface check
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
            }
            REQUIRE( tTimeSurfaceCheck );
        }
    }

//------------------------------------------------------------------------------

    SECTION( "Side Geometry Interpolation : HEX8 " )
    {
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
            mtk::Interpolation_Order tSpaceInterpolationOrder;

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
                    break;
                }
                // FIXME for QUADRATIC and CUBIC
            }

            //create a line time element
            Matrix< DDRMat > tTHat = { { 0.0 },
                                       { 1.0 } };

            // create a space and time geometry interpolation rule
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                                                Interpolation_Type::LAGRANGE,
                                                tSpaceInterpolationOrder,
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

            // booleans for checks
            bool tSideSurfaceCheck = true;
            bool tSideNormalCheck  = true;

            // loop over the element sides
            for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
            {
                // get treated side ordinal
                moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

                // init the side surface
                real tSideSurface = 0;
                Matrix< DDRMat > tNormal;

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
                    tGeomInterpolator.surf_det_J( tSurfDetJ, tNormal,
                                                  tSideIntegPointI, tSideOrdinal );

                    // add contribution to the surface
                    tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );
                }

                // check the surface value
                tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );

                // check the normals
                for( uint iNorm = 0; iNorm < 3; iNorm++ )
                {
                    tSideNormalCheck = tSideNormalCheck && ( std::abs( tNormal( iNorm ) - tSideNormalsExact( iNorm, iSide ) ) < tEpsilon );
                }
            }
            REQUIRE( tSideSurfaceCheck );
            REQUIRE( tSideNormalCheck );

            // time side geometry type and space ordinal------------------------------------
            //------------------------------------------------------------------------------

            Matrix< IndexMat > tListOfTimeOrdinals = { 0, 1 };

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
            }
            REQUIRE( tTimeSurfaceCheck );
        }
    }
//------------------------------------------------------------------------------

    SECTION( "Side Geometry Interpolation : TRI3 - TRI6 - TRI10 " )
    {
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
            mtk::Interpolation_Order tSpaceInterpolationOrder;

            switch( iOrder )
            {
                case ( 0 ):
                {
                    // define a TRI4 in the physical space, i.e. space coordinates xHat
                    tSpaceInterpolationOrder = mtk::Interpolation_Order::LINEAR;
                    tXHat = {{ 0.0, 0.0 },
                             { 1.0, 0.0 },
                             { 0.0, 1.0 }};
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
                    break;
                }
            }

            // define a LINE2 in the physical space, i.e. time coordinates tHat
            Matrix< DDRMat > tTHat = {{ 0.0 },
                                      { 2.0 }};

            // create a space and time geometry interpolation rule
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::TRI,
                                                Interpolation_Type::LAGRANGE,
                                                tSpaceInterpolationOrder,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

            // create a space and time geometry interpolator
            Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

            //set the coefficients xHat, tHat
            tGeomInterpolator.set_coeff( tXHat, tTHat );

            // side geometry type and space ordinal
            mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::LINE;
            Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2 };

            // create a side integration
            Integration_Rule tSideIntegrationRule( tSideGeometryType,
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_5,
                                                   Integration_Type::GAUSS,
                                                   Integration_Order::BAR_1 );

            // create a side integrator
            Integrator tSideIntegrator( tSideIntegrationRule );

            //get number of integration points, integration points and weights
            uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
            Matrix< DDRMat > tSideIntegPoints  = tSideIntegrator.get_points();
            Matrix< DDRMat > tSideIntegWeights = tSideIntegrator.get_weights();

            // booleans for checks
            bool tSideSurfaceCheck = true;
            bool tSideNormalCheck  = true;

            // loop over the sides
            for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
            {
                // get treated side ordinal
                moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

                // init the side surface and the normal
                real tSideSurface = 0;
                Matrix< DDRMat > tNormal;

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
                    tGeomInterpolator.surf_det_J( tSurfDetJ, tNormal,
                                                  tSideIntegPointI, tSideOrdinal );

                    // add contribution to the surface
                    tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );
                }

                // check the surface value
                tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );

                // check the normals
                for( uint iNorm = 0; iNorm < 2; iNorm++ )
                {
                    tSideNormalCheck = tSideNormalCheck && ( std::abs( tNormal( iNorm ) - tSideNormalsExact( iNorm, iSide ) ) < tEpsilon );
                }
            }

            REQUIRE( tSideSurfaceCheck );
            REQUIRE( tSideNormalCheck );

            // time side geometry type and space ordinal------------------------------------
            //------------------------------------------------------------------------------

            Matrix< IndexMat > tListOfTimeOrdinals = { 0, 1 };

            // create a side integration
            Integration_Rule tTimeSideIntegrationRule( mtk::Geometry_Type::TRI,
                                                       Integration_Type::GAUSS,
                                                       Integration_Order::TRI_6,
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
                tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 1.0 ) < tEpsilon );
            }
            REQUIRE( tTimeSurfaceCheck );
        }
    }

//------------------------------------------------------------------------------

    SECTION( "Side Geometry Interpolation : TET4 - TET10 - TET20 " )
     {
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
             mtk::Interpolation_Order tSpaceInterpolationOrder;

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
                     break;
                 }
             }

             // define a LINE2 in the physical space, i.e. time coordinates tHat
             Matrix< DDRMat > tTHat = {{ 0.0 },
                                       { 2.0 }};

             // create a space and time geometry interpolation rule
             Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::TET,
                                                 Interpolation_Type::LAGRANGE,
                                                 tSpaceInterpolationOrder,
                                                 Interpolation_Type::LAGRANGE,
                                                 mtk::Interpolation_Order::LINEAR );

             // create a space and time geometry interpolator
             Geometry_Interpolator tGeomInterpolator( tGeomInterpRule, true );

             //set the coefficients xHat, tHat
             tGeomInterpolator.set_coeff( tXHat, tTHat );

             // side geometry type and space ordinal
             mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::TRI;
             Cell< moris::moris_index > tListOfSideOrdinals = { 0, 1, 2, 3 };

             // create a side integration
             Integration_Rule tSideIntegrationRule( tSideGeometryType,
                                                    Integration_Type::GAUSS,
                                                    Integration_Order::TRI_6,
                                                    Integration_Type::GAUSS,
                                                    Integration_Order::BAR_1 );

             // create a side integrator
             Integrator tSideIntegrator( tSideIntegrationRule );

             //get number of integration points, integration points and weights
             uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
             Matrix< DDRMat > tSideIntegPoints  = tSideIntegrator.get_points();
             Matrix< DDRMat > tSideIntegWeights = tSideIntegrator.get_weights();

             // booleans for checks
             bool tSideSurfaceCheck = true;
             bool tSideNormalCheck  = true;

             // loop over the element sides
             for( uint iSide = 0; iSide < tListOfSideOrdinals.size(); iSide++ )
             {
                 // get treated side ordinal
                 moris_index tSideOrdinal = tListOfSideOrdinals( iSide );

                 // init the side surface and the normal
                 real tSideSurface = 0;
                 Matrix< DDRMat > tNormal;

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
                     tGeomInterpolator.surf_det_J( tSurfDetJ, tNormal,
                                                   tSideIntegPointI, tSideOrdinal );

                     // add contribution to the surface
                     tSideSurface = tSideSurface + tSurfDetJ * tSideIntegWeights( iGP );
                 }

                 // check the surface value
                 tSideSurfaceCheck = tSideSurfaceCheck && ( std::abs( tSideSurface - tSideSurfaceExact( iSide ) ) < tEpsilon );

                 // check the normals
                 for( uint iNorm = 0; iNorm < 3; iNorm++ )
                 {
                     tSideNormalCheck = tSideNormalCheck && ( std::abs( tNormal( iNorm ) - tSideNormalsExact( iNorm, iSide ) ) < tEpsilon );
                 }
             }
             REQUIRE( tSideSurfaceCheck );
             REQUIRE( tSideNormalCheck );


         // time side geometry type and space ordinal------------------------------------
         //------------------------------------------------------------------------------

         Matrix< IndexMat > tListOfTimeOrdinals = { 0, 1 };

         // create a side integration
         Integration_Rule tTimeSideIntegrationRule( mtk::Geometry_Type::TET,
                                                    Integration_Type::GAUSS,
                                                    Integration_Order::TET_5,
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
             tTimeSurfaceCheck = tTimeSurfaceCheck && ( std::abs( tTimeSideSurface - 1.0/3.0 ) < tEpsilon );
         }
         REQUIRE( tTimeSurfaceCheck );
         }
    }
}

