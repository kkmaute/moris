#include "catch.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"    //FEM/INT/src
#include "cl_MTK_Integrator.hpp"               //MTK/src
#include "fn_vectorize.hpp"                    //LINALG/src

using namespace moris;
using namespace fem;

TEST_CASE( "Intergration_Mesh", "[moris],[fem],[IntegMesh]" )
{
    // define an epsilon environment
    double tEpsilon = 1E-12;

    SECTION( "Interpolation mesh QUAD4 - Integration mesh QUAD4 " )
    {
        // define an interpolation mesh
        //------------------------------------------------------------------------------
        // define a QUAD4 space element, i.e. space coordinates xHat
        Matrix< DDRMat > tXHatIP = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };

        // define a line time element, i.e. time coordinates tHat
        Matrix< DDRMat > tTHatIP = { { 0.0 }, { 1.0 }, { 0.5 } };

        // the QUAD4 interpolation element in space and time param coordinates xiHat., tauHat
        Matrix< DDRMat > tXiHatIP  = { { -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 } };
        Matrix< DDRMat > tTauHatIP = { { -1.0 }, { 1.0 }, { 0.0 } };

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeoInterpIPRule( mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeoInterpIP( tGeoInterpIPRule );

        // set the coefficients xHat, tHat
        tGeoInterpIP.set_space_coeff( tXHatIP );
        tGeoInterpIP.set_time_coeff( tTHatIP );

        // set the coefficients xHat, tHat
        tGeoInterpIP.set_space_param_coeff( tXiHatIP );
        tGeoInterpIP.set_time_param_coeff( tTauHatIP );

        // space and time geometry interpolations
        mtk::Interpolation_Function_Base* tSpaceInterpolation = tGeoInterpIPRule.create_space_interpolation_function();
        mtk::Interpolation_Function_Base* tTimeInterpolation  = tGeoInterpIPRule.create_time_interpolation_function();

        uint tNumTimeBases     = tTimeInterpolation->get_number_of_bases();
        uint tNumTimeDim       = tTimeInterpolation->get_number_of_dimensions();
        uint tNumParamSpaceDim = tSpaceInterpolation->get_number_of_param_dimensions();

        // define an integration mesh
        //------------------------------------------------------------------------------
        // integration mesh geometry type
        mtk::Geometry_Type tGeoTypeIG = mtk::Geometry_Type::QUAD;

        // define a QUAD4 integration element, i.e. space param coordinates xiHat
        Matrix< DDRMat > tXiHatIG = { { -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 } };

        Matrix< DDRMat > tTauHatIG = { { -1.0 }, { 1.0 }, { 0.0 } };

        // the QUAD4 integration element in space physical coordinates xHat
        Matrix< DDRMat > tXHatIG = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
        Matrix< DDRMat > tTHatIG = { { 0.0 }, { 1.0 }, { 0.5 } };

        // integration mesh interpolation rule
        mtk::Interpolation_Rule tGeoInterpIGRule( tGeoTypeIG,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC );

        // create a space and time geometry interpolator fot the integration element
        Geometry_Interpolator tGeoInterpIG( tGeoInterpIGRule );

        // set the coefficients xHat, tHat
        tGeoInterpIG.set_space_coeff( tXHatIG );
        tGeoInterpIG.set_time_coeff( tTHatIG );

        // set the coefficients xiHat, tauHat
        tGeoInterpIG.set_space_param_coeff( tXiHatIG );
        tGeoInterpIG.set_time_param_coeff( tTauHatIG );

        // space geometry interpolations for integration mesh
        mtk::Interpolation_Function_Base* tIntegSpaceInterpolation = tGeoInterpIGRule.create_space_interpolation_function();
        uint                              tIntegNumSpaceBases      = tIntegSpaceInterpolation->get_number_of_bases();
        uint                              tIntegNumParamSpaceDim   = tIntegSpaceInterpolation->get_number_of_param_dimensions();

        // create a integration rule
        mtk::Integration_Rule tIntegrationRule( tGeoTypeIG,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::QUAD_3x3,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_3 );

        // create a side integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        // get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegPoints;
        tIntegrator.get_points( tIntegPoints );
        Matrix< DDRMat > tIntegWeights;
        tIntegrator.get_weights( tIntegWeights );

        // boolean for surface check
        bool tSurfaceCheck  = true;
        bool tSurfaceCheck2 = true;
        bool tMappedIPCheck = true;

        // init the surface of the integration mesh
        real tSurface  = 0;
        real tSurface2 = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // set the treated integration point location in the surface ref space for the geometry interp
            tGeoInterpIG.set_space_time( tIntegPointI );

            // get the treated integration point location in the interpolation space
            //------------------------------------------------------------------------------
            // unpack the space and time param coords of the integration point
            Matrix< DDRMat > tXi = tIntegPointI( { 0, tIntegNumParamSpaceDim - 1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, tIntegPointI( tIntegNumParamSpaceDim ) );

            Matrix< DDRMat > tNIntegSpace;
            Matrix< DDRMat > tNTime;

            // evaluate space interpolation shape functions at integration point
            tIntegSpaceInterpolation->eval_N( tXi, tNIntegSpace );

            // evaluate time interpolation shape functions at aParamPoint
            tTimeInterpolation->eval_N( tTau, tNTime );

            // build space time interpolation functions for integration mesh
            Matrix< DDRMat > tN = trans( vectorize( trans( tNIntegSpace ) * tNTime ) );

            // get the parametric coordinates of the integration mesh in interpolation mesh
            Matrix< DDRMat > tInterpParamCoords( tNumTimeDim + tNumParamSpaceDim, tNumTimeBases * tIntegNumSpaceBases );

            // get the time parametric coordinates
            Matrix< DDRMat > tTimeParamCoords;
            tTimeInterpolation->get_param_coords( tTimeParamCoords );

            // get a vector of ones
            Matrix< DDRMat > tOnes( 1, tIntegNumSpaceBases, 1.0 );

            // loop on the time bases
            for ( uint i = 0; i < tNumTimeBases; i++ )
            {
                // fill the space time parametric coordinates matrix with space coordinates
                tInterpParamCoords(
                        { 0, tIntegNumParamSpaceDim - 1 },
                        { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases - 1 } ) =
                        trans( tXiHatIG );

                // fill the space time parametric coordinates matrix with time coordinates
                tInterpParamCoords(
                        { tIntegNumParamSpaceDim, tIntegNumParamSpaceDim },
                        { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases - 1 } ) =
                        tTimeParamCoords( i ) * tOnes;
            }

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            Matrix< DDRMat > tRefIntegPointI = tInterpParamCoords * trans( tN );

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            const Matrix< DDRMat >& tRefIntegPointI2 = tGeoInterpIG.map_integration_point();

            for ( uint iCoords = 0; iCoords < 3; iCoords++ )
            {
                tMappedIPCheck = tMappedIPCheck && ( std::abs( tRefIntegPointI( iCoords ) - tRefIntegPointI2( iCoords ) ) < tEpsilon );
            }

            // evaluate detJ
            //------------------------------------------------------------------------------
            // set the global integration point for the IP geometry interp
            tGeoInterpIP.set_space_time( tRefIntegPointI );

            real tDetJ1 = tGeoInterpIP.det_J();

            // get the space jacobian
            Matrix< DDRMat > tdNSpacedXi;
            tIntegSpaceInterpolation->eval_dNdXi( tXi, tdNSpacedXi );
            Matrix< DDRMat > tSpaceJt = tdNSpacedXi * tXiHatIG;
            // get the time Jacobian
            Matrix< DDRMat > tdNTimedTau;
            tTimeInterpolation->eval_dNdXi( tTau, tdNTimedTau );
            Matrix< DDRMat > tTimeJt = tdNTimedTau * tTauHatIG;
            real             tDetJ2  = det( tSpaceJt ) * det( tTimeJt );

            real tDetJ = tDetJ1 * tDetJ2;

            // add contribution to the surface
            tSurface = tSurface + tDetJ * tIntegWeights( iGP );

            // add contribution to the surface from GeoInterpIG
            tSurface2 = tSurface2 + tGeoInterpIG.det_J() * tIntegWeights( iGP );
        }

        // check the integration point values
        REQUIRE( tMappedIPCheck );

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 1.0 ) < tEpsilon );

        REQUIRE( tSurfaceCheck );

        // check the surface value
        tSurfaceCheck2 = tSurfaceCheck2 && ( std::abs( tSurface2 - 1.0 ) < tEpsilon );

        REQUIRE( tSurfaceCheck2 );

        // clean up
        delete tSpaceInterpolation;
        delete tTimeInterpolation;
        delete tIntegSpaceInterpolation;
    }

    SECTION( "Interpolation mesh QUAD4 - Integration mesh TRI3 " )
    {
        // define an interpolation mesh
        //------------------------------------------------------------------------------
        // define a QUAD4 space element, i.e. space coordinates xHat
        Matrix< DDRMat > tXHatIP = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };

        // define a line time element, i.e. time coordinates tHat
        Matrix< DDRMat > tTHatIP = { { 0.0 }, { 1.0 }, { 0.5 } };

        // the QUAD4 interpolation element in space and time param coordinates xiHat., tauHat
        Matrix< DDRMat > tXiHatIP  = { { -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 } };
        Matrix< DDRMat > tTauHatIP = { { -1.0 }, { 1.0 }, { 0.0 } };

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeoInterpIPRule( mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeoInterpIP( tGeoInterpIPRule );

        // set the coefficients xHat, tHat
        tGeoInterpIP.set_space_coeff( tXHatIP );
        tGeoInterpIP.set_time_coeff( tTHatIP );

        // set the coefficients xHat, tHat
        tGeoInterpIP.set_space_param_coeff( tXiHatIP );
        tGeoInterpIP.set_time_param_coeff( tTauHatIP );

        // space and time geometry interpolations
        mtk::Interpolation_Function_Base* tSpaceInterpolation = tGeoInterpIPRule.create_space_interpolation_function();
        mtk::Interpolation_Function_Base* tTimeInterpolation  = tGeoInterpIPRule.create_time_interpolation_function();

        uint tNumTimeBases     = tTimeInterpolation->get_number_of_bases();
        uint tNumTimeDim       = tTimeInterpolation->get_number_of_dimensions();
        uint tNumParamSpaceDim = tSpaceInterpolation->get_number_of_param_dimensions();

        // define an integration mesh
        //------------------------------------------------------------------------------
        // integration mesh geometry type
        mtk::Geometry_Type tGeoTypeIG = mtk::Geometry_Type::TRI;

        // integration mesh interpolation rule
        mtk::Interpolation_Rule tGeoInterpIGRule( tGeoTypeIG,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC );

        // define a TRI3 integration element, i.e. space param coordinates xiHat
        Matrix< DDRMat > tXiHatIG  = { { -1.0, -1.0 }, { 0.0, -1.0 }, { -1.0, 1.0 } };
        Matrix< DDRMat > tTauHatIG = { { -1.0 }, { 1.0 }, { 0.0 } };

        // the TRI3 integration element in space physical coordinates xHat
        Matrix< DDRMat > tXHatIG = { { 0.0, 0.0 }, { 0.5, 0.0 }, { 0.0, 1.0 } };
        Matrix< DDRMat > tTHatIG = { { 0.0 }, { 1.0 }, { 0.5 } };

        // create a space and time geometry interpolator fot the integration element
        Geometry_Interpolator tGeoInterpIG( tGeoInterpIGRule, tGeoInterpIPRule );

        // set the coefficients xHat, tHat
        tGeoInterpIG.set_space_coeff( tXHatIG );
        tGeoInterpIG.set_time_coeff( tTHatIG );

        // set the coefficients xHat, tHat
        tGeoInterpIG.set_space_param_coeff( tXiHatIG );
        tGeoInterpIG.set_time_param_coeff( tTauHatIG );

        // space geometry interpolations for integration mesh
        mtk::Interpolation_Function_Base* tIntegSpaceInterpolation = tGeoInterpIGRule.create_space_interpolation_function();
        uint                              tIntegNumSpaceBases      = tIntegSpaceInterpolation->get_number_of_bases();
        uint                              tIntegNumParamSpaceDim   = tIntegSpaceInterpolation->get_number_of_param_dimensions();

        // create a integration rule
        mtk::Integration_Rule tIntegrationRule( tGeoTypeIG,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::TRI_3,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_3 );

        // create a side integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        // get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegPoints;
        tIntegrator.get_points( tIntegPoints );
        Matrix< DDRMat > tIntegWeights;
        tIntegrator.get_weights( tIntegWeights );

        // boolean for surface check
        bool tSurfaceCheck  = true;
        bool tSurfaceCheck2 = true;
        bool tMappedIPCheck = true;

        // init the surface of the integration mesh
        real tSurface  = 0;
        real tSurface2 = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // set the treated integration point location in the surface ref space for the geometry interp
            tGeoInterpIG.set_space_time( tIntegPointI );

            // get the treated integration point location in the interpolation space
            //------------------------------------------------------------------------------
            // unpack the space and time param coords of the integration point
            Matrix< DDRMat > tXi = tIntegPointI( { 0, tIntegNumParamSpaceDim - 1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, tIntegPointI( tIntegNumParamSpaceDim ) );

            Matrix< DDRMat > tNIntegSpace;
            Matrix< DDRMat > tNTime;

            // evaluate space interpolation shape functions at integration point
            tIntegSpaceInterpolation->eval_N( tXi, tNIntegSpace );

            // evaluate time interpolation shape functions at aParamPoint
            tTimeInterpolation->eval_N( tTau, tNTime );

            // build space time interpolation functions for integration mesh
            Matrix< DDRMat > tN = trans( vectorize( trans( tNIntegSpace ) * tNTime ) );

            // get the parametric coordinates of the integration mesh in interpolation mesh
            Matrix< DDRMat > tInterpParamCoords( tNumTimeDim + tNumParamSpaceDim, tNumTimeBases * tIntegNumSpaceBases );

            // get the time parametric coordinates
            Matrix< DDRMat > tTimeParamCoords;
            tTimeInterpolation->get_param_coords( tTimeParamCoords );

            // get a vector of ones
            Matrix< DDRMat > tOnes( 1, tIntegNumSpaceBases, 1.0 );

            // loop on the time bases
            for ( uint i = 0; i < tNumTimeBases; i++ )
            {
                // fill the space time parametric coordinates matrix with space coordinates
                tInterpParamCoords(
                        { 0, tNumParamSpaceDim - 1 }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases - 1 } ) =    //
                        trans( tXiHatIG );

                // fill the space time parametric coordinates matrix with time coordinates
                tInterpParamCoords(
                        { tNumParamSpaceDim, tNumParamSpaceDim }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases - 1 } ) =    //
                        tTimeParamCoords( i ) * tOnes;
            }

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            Matrix< DDRMat > tRefIntegPointI = tInterpParamCoords * trans( tN );

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            const Matrix< DDRMat >& tRefIntegPointI2 = tGeoInterpIG.map_integration_point();

            for ( uint iCoords = 0; iCoords < 3; iCoords++ )
            {
                tMappedIPCheck = tMappedIPCheck && ( std::abs( tRefIntegPointI( iCoords ) - tRefIntegPointI2( iCoords ) ) < tEpsilon );
            }

            // evaluate detJ
            //------------------------------------------------------------------------------
            // set the global integration point for the IP geometry interp
            tGeoInterpIP.set_space_time( tRefIntegPointI );

            real tDetJ1 = tGeoInterpIP.det_J();

            // get the space jacobian
            Matrix< DDRMat > tdNSpacedXi;
            tIntegSpaceInterpolation->eval_dNdXi( tXi, tdNSpacedXi );
            Matrix< DDRMat > tSpaceJt = tdNSpacedXi * tXiHatIG;
            // get the time Jacobian
            Matrix< DDRMat > tdNTimedTau;
            tTimeInterpolation->eval_dNdXi( tTau, tdNTimedTau );
            Matrix< DDRMat > tTimeJt = tdNTimedTau * tTauHatIG;

            Matrix< DDRMat > tSpaceJt2 = trans( tSpaceJt );
            real             detJSpace = det( tSpaceJt2 ) / 2.0;
            real             tDetJ2    = detJSpace * det( tTimeJt );

            real tDetJ = tDetJ1 * tDetJ2;

            // add contribution to the surface
            tSurface = tSurface + tDetJ * tIntegWeights( iGP );

            // add contribution to the surface from GeoInterpIG
            tSurface2 = tSurface2 + tGeoInterpIG.det_J() * tIntegWeights( iGP );
        }

        // check the integration point values
        REQUIRE( tMappedIPCheck );

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 0.25 ) < tEpsilon );
        REQUIRE( tSurfaceCheck );

        // check the surface value
        tSurfaceCheck2 = tSurfaceCheck2 && ( std::abs( tSurface2 - 0.25 ) < tEpsilon );
        REQUIRE( tSurfaceCheck2 );

        // clean up
        delete tSpaceInterpolation;
        delete tTimeInterpolation;
        delete tIntegSpaceInterpolation;
    }

    SECTION( "Interpolation mesh HEX8 - Integration mesh HEX8 " )
    {
        // define an interpolation mesh
        //------------------------------------------------------------------------------
        // define a HEX8 space element, i.e. space coordinates xHat
        Matrix< DDRMat > tXHatIP = { { 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0 },
            { 1.0, 1.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 },
            { 1.0, 0.0, 1.0 },
            { 1.0, 1.0, 1.0 },
            { 0.0, 1.0, 1.0 } };

        // define a line time element, i.e. time coordinates tHat
        Matrix< DDRMat > tTHatIP = { { 0.0 }, { 1.0 }, { 0.5 } };

        Matrix< DDRMat > tXiHatIP = { //
            { -1.0, -1.0, -1.0 },
            { 1.0, -1.0, -1.0 },
            { 1.0, 1.0, -1.0 },
            { -1.0, 1.0, -1.0 },
            { -1.0, -1.0, 1.0 },
            { 1.0, -1.0, 1.0 },
            { 1.0, 1.0, 1.0 },
            { -1.0, 1.0, 1.0 }
        };

        Matrix< DDRMat > tTauHatIP = { { -1.0 }, { 1.0 }, { 0.0 } };

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeoInterpIPRule(
                mtk::Geometry_Type::HEX,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeoInterpIP( tGeoInterpIPRule );

        // set the coefficients xHat, tHat
        tGeoInterpIP.set_space_coeff( tXHatIP );
        tGeoInterpIP.set_time_coeff( tTHatIP );

        // set the coefficients xHat, tHat
        tGeoInterpIP.set_space_param_coeff( tXiHatIP );
        tGeoInterpIP.set_time_param_coeff( tTauHatIP );

        // space and time geometry interpolations
        mtk::Interpolation_Function_Base* tSpaceInterpolation =
                tGeoInterpIPRule.create_space_interpolation_function();

        mtk::Interpolation_Function_Base* tTimeInterpolation =
                tGeoInterpIPRule.create_time_interpolation_function();

        uint tNumTimeBases     = tTimeInterpolation->get_number_of_bases();
        uint tNumTimeDim       = tTimeInterpolation->get_number_of_dimensions();
        uint tNumParamSpaceDim = tSpaceInterpolation->get_number_of_param_dimensions();

        // define an integration mesh
        //------------------------------------------------------------------------------
        // integration mesh geometry type
        mtk::Geometry_Type tGeoTypeIG = mtk::Geometry_Type::HEX;

        // define a HEX8 integration element, i.e. space param coordinates xiHat
        Matrix< DDRMat > tXiHatIG = { //
            { -1.0, -1.0, -1.0 },
            { 1.0, -1.0, -1.0 },
            { 1.0, 1.0, -1.0 },
            { -1.0, 1.0, -1.0 },
            { -1.0, -1.0, 1.0 },
            { 1.0, -1.0, 1.0 },
            { 1.0, 1.0, 1.0 },
            { -1.0, 1.0, 1.0 }
        };

        Matrix< DDRMat > tTauHatIG = { { -1.0 }, { 1.0 }, { 0.0 } };

        // the HEX8 integration element in param  coordinates xiHat, tauHat
        Matrix< DDRMat > tXHatIG{
            { 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0 },
            { 1.0, 1.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 },
            { 1.0, 0.0, 1.0 },
            { 1.0, 1.0, 1.0 },
            { 0.0, 1.0, 1.0 }
        };
        Matrix< DDRMat > tTHatIG = { { 0.0 }, { 1.0 }, { 0.5 } };

        // integration mesh interpolation rule
        mtk::Interpolation_Rule tGeoInterpIGRule(
                tGeoTypeIG,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC );

        // create a space and time geometry interpolator fot the integration element
        Geometry_Interpolator tGeoInterpIG( tGeoInterpIGRule );

        // set the coefficients xHat, tHat
        tGeoInterpIG.set_space_coeff( tXHatIG );
        tGeoInterpIG.set_time_coeff( tTHatIG );

        // set the coefficients xHat, tHat
        tGeoInterpIG.set_space_param_coeff( tXiHatIG );
        tGeoInterpIG.set_time_param_coeff( tTauHatIG );

        // space geometry interpolations for integration mesh
        mtk::Interpolation_Function_Base* tIntegSpaceInterpolation =    //
                tGeoInterpIGRule.create_space_interpolation_function();

        uint tIntegNumSpaceBases    = tIntegSpaceInterpolation->get_number_of_bases();
        uint tIntegNumParamSpaceDim = tIntegSpaceInterpolation->get_number_of_param_dimensions();

        // create a integration rule
        mtk::Integration_Rule tIntegrationRule( tGeoTypeIG,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::HEX_3x3x3,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_3 );

        // create a side integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        // get number of integration points, integration points and weights
        uint             tNumOfIntegPoints = tIntegrator.get_number_of_points();
        Matrix< DDRMat > tIntegPoints;
        tIntegrator.get_points( tIntegPoints );
        Matrix< DDRMat > tIntegWeights;
        tIntegrator.get_weights( tIntegWeights );

        // boolean for surface check
        bool tSurfaceCheck  = true;
        bool tSurfaceCheck2 = true;
        bool tMappedIPCheck = true;

        // init the surface of the integration mesh
        real tSurface  = 0;
        real tSurface2 = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // set the treated integration point location in the surface ref space for the geometry interp
            tGeoInterpIG.set_space_time( tIntegPointI );

            // get the treated integration point location in the interpolation space
            //------------------------------------------------------------------------------
            // unpack the space and time param coords of the integration point
            Matrix< DDRMat > tXi = tIntegPointI( { 0, tIntegNumParamSpaceDim - 1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, tIntegPointI( tIntegNumParamSpaceDim ) );

            Matrix< DDRMat > tNIntegSpace;
            Matrix< DDRMat > tNTime;

            // evaluate space interpolation shape functions at integration point
            tIntegSpaceInterpolation->eval_N( tXi, tNIntegSpace );

            // evaluate time interpolation shape functions at aParamPoint
            tTimeInterpolation->eval_N( tTau, tNTime );

            // build space time interpolation functions for integration mesh
            Matrix< DDRMat > tN = trans( vectorize( trans( tNIntegSpace ) * tNTime ) );

            // get the parametric coordinates of the integration mesh in interpolation mesh
            Matrix< DDRMat > tInterpParamCoords( tNumTimeDim + tNumParamSpaceDim, tNumTimeBases * tIntegNumSpaceBases );

            // get the time parametric coordinates
            Matrix< DDRMat > tTimeParamCoords;
            tTimeInterpolation->get_param_coords( tTimeParamCoords );

            // get a vector of ones
            Matrix< DDRMat > tOnes( 1, tIntegNumSpaceBases, 1.0 );

            // loop on the time bases
            for ( uint i = 0; i < tNumTimeBases; i++ )
            {
                // fill the space time parametric coordinates matrix with space coordinates
                tInterpParamCoords(
                        { 0, tIntegNumParamSpaceDim - 1 },
                        { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases - 1 } ) =
                        trans( tXiHatIG );

                // fill the space time parametric coordinates matrix with time coordinates
                tInterpParamCoords(
                        { tIntegNumParamSpaceDim, tIntegNumParamSpaceDim },
                        { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases - 1 } ) =
                        tTimeParamCoords( i ) * tOnes;
            }

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            Matrix< DDRMat > tRefIntegPointI = tInterpParamCoords * trans( tN );

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            const Matrix< DDRMat >& tRefIntegPointI2 = tGeoInterpIG.map_integration_point();

            for ( uint iCoords = 0; iCoords < 4; iCoords++ )
            {
                tMappedIPCheck = tMappedIPCheck && ( std::abs( tRefIntegPointI( iCoords ) - tRefIntegPointI2( iCoords ) ) < tEpsilon );
            }

            // evaluate detJ
            //------------------------------------------------------------------------------
            // set the global integration point for the IP geometry interp
            tGeoInterpIP.set_space_time( tRefIntegPointI );

            real tDetJ1 = tGeoInterpIP.det_J();

            // get the space jacobian
            Matrix< DDRMat > tdNSpacedXi;
            tIntegSpaceInterpolation->eval_dNdXi( tXi, tdNSpacedXi );
            Matrix< DDRMat > tSpaceJt = tdNSpacedXi * tXiHatIG;
            // get the time Jacobian
            Matrix< DDRMat > tdNTimedTau;
            tTimeInterpolation->eval_dNdXi( tTau, tdNTimedTau );
            Matrix< DDRMat > tTimeJt = tdNTimedTau * tTauHatIG;
            real             tDetJ2  = det( tSpaceJt ) * det( tTimeJt );

            real tDetJ = tDetJ1 * tDetJ2;

            // add contribution to the surface
            tSurface = tSurface + tDetJ * tIntegWeights( iGP );

            // add contribution to the surface from GeoInterpIG
            tSurface2 = tSurface2 + tGeoInterpIG.det_J() * tIntegWeights( iGP );
        }

        // check integration points values
        REQUIRE( tMappedIPCheck );

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 1.0 ) < tEpsilon );

        REQUIRE( tSurfaceCheck );

        // check the surface value
        tSurfaceCheck2 = tSurfaceCheck2 && ( std::abs( tSurface2 - 1.0 ) < tEpsilon );

        REQUIRE( tSurfaceCheck2 );

        // clean up
        delete tSpaceInterpolation;
        delete tTimeInterpolation;
        delete tIntegSpaceInterpolation;
    }

    SECTION( "Interpolation mesh HEX8 - Integration mesh TET4 " )
    {
        // define an interpolation mesh
        //------------------------------------------------------------------------------
        // define a HEX8 space element, i.e. space coordinates xHat
        Matrix< DDRMat > tXHatIP = { //
            { 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0 },
            { 1.0, 1.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 },
            { 1.0, 0.0, 1.0 },
            { 1.0, 1.0, 1.0 },
            { 0.0, 1.0, 1.0 }
        };

        // define a line time element, i.e. time coordinates tHat
        Matrix< DDRMat > tTHatIP = { { 0.0 },
            { 1.0 },
            { 0.5 } };

        // the HEX8 interpolation element in param  coordinates xiHat, tauHat
        Matrix< DDRMat > tXiHatIP = { //
            { -1.0, -1.0, -1.0 },
            { 0.0, -1.0, -1.0 },
            { 0.0, 1.0, -1.0 },
            { -1.0, 1.0, -1.0 },
            { -1.0, -1.0, 1.0 },
            { 0.0, -1.0, 1.0 },
            { 0.0, 1.0, 1.0 },
            { -1.0, 1.0, 1.0 }
        };

        Matrix< DDRMat > tTauHatIP = { { -1.0 }, { 1.0 }, { 0.0 } };

        // create a space and time geometry interpolation rule
        mtk::Interpolation_Rule tGeoInterpIPRule( mtk::Geometry_Type::HEX,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC );

        // create a space and time geometry interpolator
        Geometry_Interpolator tGeoInterpIP( tGeoInterpIPRule );

        // set the coefficients xHat, tHat
        tGeoInterpIP.set_space_coeff( tXHatIP );
        tGeoInterpIP.set_time_coeff( tTHatIP );

        // set the coefficients xiHat, tauHat
        tGeoInterpIP.set_space_param_coeff( tXiHatIP );
        tGeoInterpIP.set_time_param_coeff( tTauHatIP );

        // space and time geometry interpolations
        mtk::Interpolation_Function_Base* tSpaceInterpolation =    //
                tGeoInterpIPRule.create_space_interpolation_function();

        mtk::Interpolation_Function_Base* tTimeInterpolation =    //
                tGeoInterpIPRule.create_time_interpolation_function();

        uint tNumTimeBases     = tTimeInterpolation->get_number_of_bases();
        uint tNumTimeDim       = tTimeInterpolation->get_number_of_dimensions();
        uint tNumParamSpaceDim = tSpaceInterpolation->get_number_of_param_dimensions();

        // define an integration mesh
        //------------------------------------------------------------------------------
        // integration mesh geometry type
        mtk::Geometry_Type tGeoTypeIG = mtk::Geometry_Type::TET;

        // integration mesh interpolation rule
        mtk::Interpolation_Rule tGeoInterpIGRule( tGeoTypeIG,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::QUADRATIC );

        // define a TET4 integration element, i.e. space param coordinates xiHat
        Matrix< DDRMat > tXiHatIG = { //
            { -1.0, -1.0, -1.0 },
            { 0.0, -1.0, -1.0 },
            { -1.0, 1.0, -1.0 },
            { -1.0, -1.0, 1.0 }
        };

        Matrix< DDRMat > tTauHatIG = { { -1.0 }, { 1.0 }, { 0.0 } };

        // the TET4 integration element in space physical coordinates xHat
        Matrix< DDRMat > tXHatIG = { //
            { 0.0, 0.0, 0.0 },
            { 0.5, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { 0.0, 0.0, 1.0 }
        };

        Matrix< DDRMat > tTHatIG = { { 0.0 }, { 1.0 }, { 0.5 } };

        // create a space and time geometry interpolator fot the integration element
        Geometry_Interpolator tGeoInterpIG( tGeoInterpIGRule, tGeoInterpIPRule );

        // set the coefficients xHat, tHat
        tGeoInterpIG.set_space_coeff( tXHatIG );
        tGeoInterpIG.set_time_coeff( tTHatIG );

        // set the coefficients xiHat, tauHat
        tGeoInterpIG.set_space_param_coeff( tXiHatIG );
        tGeoInterpIG.set_time_param_coeff( tTauHatIG );

        // space geometry interpolations for integration mesh
        mtk::Interpolation_Function_Base* tIntegSpaceInterpolation =    //
                tGeoInterpIGRule.create_space_interpolation_function();

        uint tIntegNumSpaceBases    = tIntegSpaceInterpolation->get_number_of_bases();
        uint tIntegNumParamSpaceDim = tIntegSpaceInterpolation->get_number_of_param_dimensions();

        // create a integration rule
        mtk::Integration_Rule tIntegrationRule(
                tGeoTypeIG,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::TET_5,
                mtk::Integration_Type::GAUSS,
                mtk::Integration_Order::BAR_3 );

        // create a side integrator
        mtk::Integrator tIntegrator( tIntegrationRule );

        // get number of integration points, integration points and weights
        uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

        Matrix< DDRMat > tIntegPoints;
        tIntegrator.get_points( tIntegPoints );

        Matrix< DDRMat > tIntegWeights;
        tIntegrator.get_weights( tIntegWeights );

        // boolean for surface check
        bool tSurfaceCheck  = true;
        bool tSurfaceCheck2 = true;
        bool tMappedIPCheck = true;

        // init the surface of the integration mesh
        real tSurface  = 0;
        real tSurface2 = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // set the treated integration point location in the surface ref space for the geometry interp
            tGeoInterpIG.set_space_time( tIntegPointI );

            // get the treated integration point location in the interpolation space
            //------------------------------------------------------------------------------
            // unpack the space and time param coords of the integration point
            Matrix< DDRMat > tXi = tIntegPointI( { 0, tIntegNumParamSpaceDim - 1 }, { 0, 0 } );
            Matrix< DDRMat > tTau( 1, 1, tIntegPointI( tIntegNumParamSpaceDim ) );

            Matrix< DDRMat > tNIntegSpace;
            Matrix< DDRMat > tNTime;

            // evaluate space interpolation shape functions at integration point
            tIntegSpaceInterpolation->eval_N( tXi, tNIntegSpace );

            // evaluate time interpolation shape functions at aParamPoint
            tTimeInterpolation->eval_N( tTau, tNTime );

            // build space time interpolation functions for integration mesh
            Matrix< DDRMat > tN = trans( vectorize( trans( tNIntegSpace ) * tNTime ) );

            // get the parametric coordinates of the integration mesh in interpolation mesh
            Matrix< DDRMat > tInterpParamCoords( tNumTimeDim + tNumParamSpaceDim, tNumTimeBases * tIntegNumSpaceBases );

            // get the time parametric coordinates
            Matrix< DDRMat > tTimeParamCoords;
            tTimeInterpolation->get_param_coords( tTimeParamCoords );

            // get a vector of ones
            Matrix< DDRMat > tOnes( 1, tIntegNumSpaceBases, 1.0 );

            // loop on the time bases
            for ( uint i = 0; i < tNumTimeBases; i++ )
            {
                // fill the space time parametric coordinates matrix with space coordinates
                tInterpParamCoords(
                        { 0, tNumParamSpaceDim - 1 }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases - 1 } ) =    //
                        trans( tXiHatIG );

                // fill the space time parametric coordinates matrix with time coordinates
                tInterpParamCoords(
                        { tNumParamSpaceDim, tNumParamSpaceDim }, { i * tIntegNumSpaceBases, ( i + 1 ) * tIntegNumSpaceBases - 1 } ) =    //
                        tTimeParamCoords( i ) * tOnes;
            }

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            Matrix< DDRMat > tRefIntegPointI = tInterpParamCoords * trans( tN );

            // compute the parametric coordinates of the SideParamPoint in the parent reference element
            const Matrix< DDRMat >& tRefIntegPointI2 = tGeoInterpIG.map_integration_point();

            for ( uint iCoords = 0; iCoords < 4; iCoords++ )
            {
                tMappedIPCheck =
                        tMappedIPCheck &&    //
                        ( std::abs( tRefIntegPointI( iCoords ) - tRefIntegPointI2( iCoords ) ) < tEpsilon );
            }

            // evaluate detJ
            //------------------------------------------------------------------------------
            // set the global integration point for the IP geometry interp
            tGeoInterpIP.set_space_time( tRefIntegPointI );

            real tDetJ1 = tGeoInterpIP.det_J();

            // get the space jacobian
            Matrix< DDRMat > tdNSpacedXi;
            tIntegSpaceInterpolation->eval_dNdXi( tXi, tdNSpacedXi );
            Matrix< DDRMat > tSpaceJt = tdNSpacedXi * tXiHatIG;

            // get the time Jacobian
            Matrix< DDRMat > tdNTimedTau;
            tTimeInterpolation->eval_dNdXi( tTau, tdNTimedTau );
            Matrix< DDRMat > tTimeJt = tdNTimedTau * tTauHatIG;

            Matrix< DDRMat > tSpaceJt2( tIntegNumParamSpaceDim, tIntegNumParamSpaceDim, 1.0 );
            tSpaceJt2(
                    { 1, tIntegNumParamSpaceDim - 1 }, { 0, tIntegNumParamSpaceDim - 1 } ) =    //
                    trans( tSpaceJt );

            real detJSpace = det( tSpaceJt2 ) / 6.0;
            real tDetJ2    = detJSpace * det( tTimeJt );

            real tDetJ = tDetJ1 * tDetJ2;

            // add contribution to the surface
            tSurface = tSurface + tDetJ * tIntegWeights( iGP );

            // add contribution to the surface from GeoInterpIG
            tSurface2 = tSurface2 + tGeoInterpIG.det_J() * tIntegWeights( iGP );
        }

        // check integration points values
        REQUIRE( tMappedIPCheck );

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 0.083333333334 ) < tEpsilon );

        REQUIRE( tSurfaceCheck );

        // check the surface value
        tSurfaceCheck2 = tSurfaceCheck2 && ( std::abs( tSurface2 - 0.083333333334 ) < tEpsilon );

        REQUIRE( tSurfaceCheck2 );

        // clean up
        delete tSpaceInterpolation;
        delete tTimeInterpolation;
        delete tIntegSpaceInterpolation;
    }
}
