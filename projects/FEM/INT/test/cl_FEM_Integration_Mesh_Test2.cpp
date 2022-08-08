#include "catch.hpp"
// FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_MTK_Integrator.hpp"
// LINALG/src
#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "op_div.hpp"

using namespace moris;
using namespace fem;

// Testing the integration rule on the side
TEST_CASE( "Interpolation mesh QUAD4 - Integration mesh QUAD4 ",
        "[moris],[fem],[IntegMeshQUAD4]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // define an interpolation mesh
    //------------------------------------------------------------------------------
    // define a QUAD4 space element, i.e. space coordinates xHat
    Matrix< DDRMat > tXHatIP = { //
        { 0.0, 0.0 },
        { 1.0, 0.0 },
        { 1.0, 1.0 },
        { 0.0, 1.0 }
    };

    // define a line time element, i.e. time coordinates tHat
    Matrix< DDRMat > tTHatIP = { { 0.0 }, { 1.0 }, { 0.5 } };

    // create a space and time geometry interpolation rule
    mtk::Interpolation_Rule tGeoInterpIPRule(
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC );

    // create a space and time geometry interpolator
    Geometry_Interpolator tGeoInterpIP( tGeoInterpIPRule );

    // set the coefficients xHat, tHat
    tGeoInterpIP.set_space_coeff( tXHatIP );
    tGeoInterpIP.set_time_coeff( tTHatIP );

    // define an integration mesh
    //------------------------------------------------------------------------------
    // integration mesh geometry type
    mtk::Geometry_Type tSideGeoTypeIG = mtk::Geometry_Type::LINE;

    // define a QUAD4 integration element, i.e. space param coordinates xiHat
    Matrix< DDRMat > tXiHatIG = { { -1.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 1.0 }, { -1.0, 1.0 } };

    // the QUAD4 integration element in space physical coordinates xHat
    Matrix< DDRMat > tXHatIG = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };

    Matrix< DDRMat > tTHatIG = { { 0.0 }, { 1.0 } };

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule(
            tSideGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp =
            tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule for the side
    mtk::Integration_Rule tSideIntegRule(
            tSideGeoTypeIG,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_2,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_3 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );
    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // get side phys and param coords
    Cell< moris_index > tAllSideOrdinals = { 0, 1, 2, 3 };

    // boolean for surface check
    bool tSurfaceCheck    = true;
    bool tIntegPointCheck = true;

    for ( uint iSide = 0; iSide < 4; iSide++ )
    {
        // get the treated side ordinal
        uint tSideOrdinal = tAllSideOrdinals( iSide );

        // get the node ids associated to the side ordinal
        Matrix< DDSMat > tElementNodes = { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
        Matrix< DDSMat > tSideNodes    = tElementNodes.get_row( tSideOrdinal );

        // phys coords amd parm coords in IP param space for the side
        Matrix< DDRMat > tSidePhysCoords( tSideNodes.numel(), 2 );
        Matrix< DDRMat > tSideParamCoords( tSideNodes.numel(), 2 );
        for ( uint iNode = 0; iNode < tSideNodes.numel(); iNode++ )
        {
            tSidePhysCoords.get_row( iNode )  = tXHatIG.get_row( tSideNodes( iNode ) );
            tSideParamCoords.get_row( iNode ) = tXiHatIG.get_row( tSideNodes( iNode ) );
        }

        // get the side normal
        Matrix< DDRMat > tTangent = tSidePhysCoords.get_row( 1 ) - tSidePhysCoords.get_row( 0 );
        Matrix< DDRMat > tNormal  = { { tTangent( 1 ) }, { -tTangent( 0 ) } };
        tNormal                   = tNormal / norm( tNormal );

        // get the side detJ
        real tSideDetJ = norm( tTangent ) / 2.0;

        // get the time detJ
        real tTimeDetJ = ( tTHatIG( 1 ) - tTHatIG( 0 ) ) / 2.0;

        // init the surface of the integration mesh
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // get the integration point in the IP parametric space
            // via a lookup table
            Matrix< DDRMat > tRefIntegPointI( 3, 1, 0.0 );
            switch ( tSideOrdinal )
            {
                case 0:
                    tRefIntegPointI = { { tIntegPointI( 0 ) }, { -1.0 }, { tIntegPointI( 1 ) } };
                    break;
                case 1:
                    tRefIntegPointI = { { 1.0 }, { tIntegPointI( 0 ) }, { tIntegPointI( 1 ) } };
                    break;
                case 2:
                    tRefIntegPointI = { { -tIntegPointI( 0 ) }, { 1.0 }, { tIntegPointI( 1 ) } };
                    break;
                case 3:
                    tRefIntegPointI = { { -1.0 }, { -tIntegPointI( 0 ) }, { tIntegPointI( 1 ) } };
                    break;
                default:
                    MORIS_ASSERT( false, "Wrong side ordinal." );
                    break;
            }

            // get the integration point in the IP parametric space
            // via the side shape functions
            Matrix< DDRMat > tRefIntegPointI2( 3, 1, 0.0 );
            Matrix< DDRMat > tN;

            tSideSpaceInterp->eval_N( tIntegPointI, tN );

            tRefIntegPointI2( { 0, 1 }, { 0, 0 } ) = trans( tN * tSideParamCoords );
            tRefIntegPointI2( 2 )                  = tIntegPointI( 1 );

            // check the integration point in the IP parametric space
            for ( uint iCoords = 0; iCoords < 3; iCoords++ )
            {
                tIntegPointCheck = tIntegPointCheck && ( std::abs( tRefIntegPointI( iCoords ) - tRefIntegPointI2( iCoords ) ) < tEpsilon );
            }

            // evaluate detJ
            real tDetJ = tSideDetJ * tTimeDetJ;

            // add contribution to the surface
            tSurface = tSurface + tDetJ * tIntegWeights( iGP );
        }

        // check the integration points coords
        REQUIRE( tIntegPointCheck );

        // check the surface value
        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 1.0 ) < tEpsilon );
        REQUIRE( tSurfaceCheck );
    }

    // clean up
    delete tSideSpaceInterp;
}

TEST_CASE( "Interpolation mesh TRI3 - Integration mesh TRI3 ",
        "[moris],[fem],[IntegMeshTRI3]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // define an interpolation mesh
    //------------------------------------------------------------------------------
    // define a TRI3 space element, i.e. space coordinates xHat
    Matrix< DDRMat > tXHatIP = {
        { 0.0, 0.0 },
        { 1.0, 0.0 },
        { 0.0, 1.0 }
    };

    // define a line time element, i.e. time coordinates tHat
    Matrix< DDRMat > tTHatIP = { { 0.0 }, { 1.0 }, { 0.5 } };

    // create a space and time geometry interpolation rule
    mtk::Interpolation_Rule tGeoInterpIPRule(
            mtk::Geometry_Type::TRI,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC );

    // create a space and time geometry interpolator
    Geometry_Interpolator tGeoInterpIP( tGeoInterpIPRule );

    // set the coefficients xHat, tHat
    tGeoInterpIP.set_space_coeff( tXHatIP );
    tGeoInterpIP.set_time_coeff( tTHatIP );

    // define an integration mesh
    //------------------------------------------------------------------------------
    // integration mesh geometry type
    mtk::Geometry_Type tSideGeoTypeIG = mtk::Geometry_Type::LINE;

    // define a TRI3 integration element, i.e. space param coordinates xiHat
    Matrix< DDRMat > tXiHatIG = {
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 }
    };

    // the QUAD4 integration element in space physical coordinates xHat
    Matrix< DDRMat > tXHatIG = {
        { 0.0, 0.0 },
        { 1.0, 0.0 },
        { 0.0, 1.0 }
    };

    Matrix< DDRMat > tTHatIG = { { 0.0 }, { 1.0 } };

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule(
            tSideGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp = tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule for the side
    mtk::Integration_Rule tSideIntegRule(
            tSideGeoTypeIG,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_3,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_3 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint tNumOfIntegPoints = tSideIntegrator.get_number_of_points();

    Matrix< DDRMat > tIntegPoints;

    tSideIntegrator.get_points( tIntegPoints );

    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // get side phys and param coords
    Cell< moris_index > tAllSideOrdinals = { 0, 1, 2 };

    // boolean for surface check
    bool tSurfaceCheck    = true;
    bool tIntegPointCheck = true;

    for ( uint iSide = 0; iSide < 3; iSide++ )
    {
        // get the treated side ordinal
        uint tSideOrdinal = tAllSideOrdinals( iSide );

        // get the node ids associated to the side ordinal
        Matrix< DDSMat > tElementNodes = { { 0, 1 }, { 1, 2 }, { 2, 0 } };
        Matrix< DDSMat > tSideNodes    = tElementNodes.get_row( tSideOrdinal );

        // phys coords amd parm coords in IP param space for the side
        Matrix< DDRMat > tSidePhysCoords( tSideNodes.numel(), 2 );
        Matrix< DDRMat > tSideParamCoords( tSideNodes.numel(), 3 );
        for ( uint iNode = 0; iNode < tSideNodes.numel(); iNode++ )
        {
            tSidePhysCoords.get_row( iNode )  = tXHatIG.get_row( tSideNodes( iNode ) );
            tSideParamCoords.get_row( iNode ) = tXiHatIG.get_row( tSideNodes( iNode ) );
        }

        // get the side normal
        Matrix< DDRMat > tTangent = tSidePhysCoords.get_row( 1 ) - tSidePhysCoords.get_row( 0 );
        Matrix< DDRMat > tNormal  = { { tTangent( 1 ) }, { -tTangent( 0 ) } };
        tNormal                   = tNormal / norm( tNormal );

        // get the side detJ
        real tSideDetJ = norm( tTangent ) / 2.0;

        // get the time detJ
        real tTimeDetJ = ( tTHatIG( 1 ) - tTHatIG( 0 ) ) / 2.0;

        // init the surface of the integration mesh
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // get the integration point in the IP parametric space
            // via lookup tables
            Matrix< DDRMat > tRefIntegPointI( 4, 1, 0.0 );

            switch ( tSideOrdinal )
            {
                case 0:
                    tRefIntegPointI = { { 1.0 + ( -1.0 - tIntegPointI( 0 ) ) / 2.0 }, { -( -1.0 - tIntegPointI( 0 ) ) / 2.0 }, { 0.0 }, { tIntegPointI( 1 ) } };
                    break;
                case 1:
                    tRefIntegPointI = { { 0.0 }, { 1.0 + ( -1.0 - tIntegPointI( 0 ) ) / 2.0 }, { -( -1.0 - tIntegPointI( 0 ) ) / 2.0 }, { tIntegPointI( 1 ) } };
                    break;
                case 2:
                    tRefIntegPointI = { { -( -1.0 - tIntegPointI( 0 ) ) / 2.0 }, { 0.0 }, { 1.0 + ( -1.0 - tIntegPointI( 0 ) ) / 2.0 }, { tIntegPointI( 1 ) } };
                    break;
                default:
                    MORIS_ASSERT( false, "Wrong side ordinal." );
                    break;
            }

            // get the integration point in the IP parametric space
            // via the side shape functions
            Matrix< DDRMat > tRefIntegPointI2( 4, 1, 0.0 );
            Matrix< DDRMat > tN;

            tSideSpaceInterp->eval_N( tIntegPointI, tN );

            tRefIntegPointI2( { 0, 2 }, { 0, 0 } ) = trans( tN * tSideParamCoords );
            tRefIntegPointI2( 3 )                  = tIntegPointI( 1 );

            // check the integration point in the IP parametric space
            for ( uint iCoords = 0; iCoords < 3; iCoords++ )
            {
                tIntegPointCheck = tIntegPointCheck && ( std::abs( tRefIntegPointI( iCoords ) - tRefIntegPointI2( iCoords ) ) < tEpsilon );
            }

            // evaluate detJ
            real tDetJ = tSideDetJ * tTimeDetJ;

            // add contribution to the surface
            tSurface = tSurface + tDetJ * tIntegWeights( iGP );
        }

        // check the integration points coords
        REQUIRE( tIntegPointCheck );

        // check the surface value
        Matrix< DDRMat > tSurfaceExact = { { 1.0 }, { std::sqrt( 2.0 ) }, { 1.0 } };

        tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - tSurfaceExact( iSide ) ) < tEpsilon );

        REQUIRE( tSurfaceCheck );
    }

    // clean up
    delete tSideSpaceInterp;
}

TEST_CASE( "Interpolation mesh TET4 - Integration mesh TET4 ",
        "[moris],[fem],[IntegMeshTET4]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // define an interpolation mesh
    //------------------------------------------------------------------------------
    // define a HEX8 space element, i.e. space coordinates xHat
    Matrix< DDRMat > tXHatIP = {
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 }
    };

    // define a line time element, i.e. time coordinates tHat
    Matrix< DDRMat > tTHatIP = { { 0.0 },
        { 1.0 },
        { 0.5 } };

    // create a space and time geometry interpolation rule
    mtk::Interpolation_Rule tGeoInterpIPRule(
            mtk::Geometry_Type::TET,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC );

    // create a space and time geometry interpolator
    Geometry_Interpolator tGeoInterpIP( tGeoInterpIPRule );

    // set the coefficients xHat, tHat
    tGeoInterpIP.set_space_coeff( tXHatIP );
    tGeoInterpIP.set_time_coeff( tTHatIP );

    // define an integration mesh
    //------------------------------------------------------------------------------
    // integration mesh geometry type
    mtk::Geometry_Type tSideGeoTypeIG = mtk::Geometry_Type::TRI;

    // define a TET4 integration element, i.e. space param coordinates xiHat
    Matrix< DDRMat > tXiHatIG = {
        { 1.0, 0.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0, 0.0 },
        { 0.0, 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 0.0, 1.0 }
    };

    // the HEX8 integration element in space physical coordinates xHat
    Matrix< DDRMat > tXHatIG = {
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 }
    };

    Matrix< DDRMat > tTHatIG = { { 0.0 }, { 1.0 } };

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule(
            tSideGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp = tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule
    mtk::Integration_Rule tSideIntegRule(
            tSideGeoTypeIG,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TRI_3,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_3 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );
    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // get side phys and param coords
    Cell< moris_index > tAllSideOrdinals = { 0, 1, 2, 3 };

    // boolean for surface check
    bool tSurfaceCheck    = true;
    bool tIntegPointCheck = true;

    for ( uint iSide = 0; iSide < 4; iSide++ )
    {
        // get the treated side ordinal
        uint tSideOrdinal = tAllSideOrdinals( iSide );

        // get the node ids associated to the side ordinal
        Matrix< DDSMat > tElementNodes = { { 0, 1, 3 }, { 1, 2, 3 }, { 0, 3, 2 }, { 0, 2, 1 } };
        Matrix< DDSMat > tSideNodes    = tElementNodes.get_row( tSideOrdinal );

        // phys coords
        Matrix< DDRMat > tSidePhysCoords( tSideNodes.numel(), 3 );
        Matrix< DDRMat > tSideParamCoords( tSideNodes.numel(), 4 );
        for ( uint iNode = 0; iNode < tSideNodes.numel(); iNode++ )
        {
            tSidePhysCoords.get_row( iNode )  = tXHatIG.get_row( tSideNodes( iNode ) );
            tSideParamCoords.get_row( iNode ) = tXiHatIG.get_row( tSideNodes( iNode ) );
        }

        // get the side normal
        Matrix< DDRMat > tTangent1 = tSidePhysCoords.get_row( 1 ) - tSidePhysCoords.get_row( 0 );
        Matrix< DDRMat > tTangent2 = tSidePhysCoords.get_row( 2 ) - tSidePhysCoords.get_row( 0 );
        Matrix< DDRMat > tNormal   = cross( tTangent1, tTangent2 );
        tNormal                    = tNormal / norm( tNormal );

        // get the side detJ
        real tSideDetJ = ( norm( cross( tTangent1, tTangent2 ) ) / 2.0 ) / 1.0;

        // get the time detJ
        real tTimeDetJ = ( tTHatIG( 1 ) - tTHatIG( 0 ) ) / 2.0;

        // init the surface of the integration mesh
        real tSurface = 0;

        // loop over the integration points
        for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
        {
            // get the treated integration point location in integration space
            Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

            // get the integration point in the IP parametric space
            // via lookup tables
            Matrix< DDRMat > tRefIntegPointI( 5, 1, 0.0 );
            switch ( tSideOrdinal )
            {
                case 0:
                    tRefIntegPointI = { //
                        { tIntegPointI( 0 ) },
                        { tIntegPointI( 1 ) },
                        { 0.0 },
                        { 1.0 - tIntegPointI( 0 ) - tIntegPointI( 1 ) },
                        { tIntegPointI( 2 ) }
                    };
            };
            break;
            case 1:
                tRefIntegPointI = { //
                    { 0.0 },
                    { tIntegPointI( 0 ) },
                    { tIntegPointI( 1 ) },
                    { 1.0 - tIntegPointI( 0 ) - tIntegPointI( 1 ) },
                    { tIntegPointI( 2 ) }
                };
                break;
            case 2:
                tRefIntegPointI = { //
                    { tIntegPointI( 0 ) },
                    { 0.0 },
                    { 1.0 - tIntegPointI( 0 ) - tIntegPointI( 1 ) },
                    { tIntegPointI( 1 ) },
                    { tIntegPointI( 2 ) }
                };
                break;
            case 3:
                tRefIntegPointI = { //
                    { tIntegPointI( 0 ) },
                    { 1.0 - tIntegPointI( 0 ) - tIntegPointI( 1 ) },
                    { tIntegPointI( 1 ) },
                    { 0.0 },
                    { tIntegPointI( 2 ) }
                };
                break;
            default:
                MORIS_ASSERT( false, "Wrong side ordinal." );
                break;
        }

        // get the integration point in the IP parametric space
        // via the side shape functions
        Matrix< DDRMat > tRefIntegPointI2( 5, 1, 0.0 );
        Matrix< DDRMat > tN;

        tSideSpaceInterp->eval_N( tIntegPointI, tN );

        tRefIntegPointI2( { 0, 3 }, { 0, 0 } ) = trans( tN * tSideParamCoords );
        tRefIntegPointI2( 4 )                  = tIntegPointI( 2 );

        // check the integration point in the IP parametric space
        for ( uint iCoords = 0; iCoords < 5; iCoords++ )
        {
            tIntegPointCheck =
                    tIntegPointCheck &&    //
                    ( std::abs( tRefIntegPointI( iCoords ) - tRefIntegPointI2( iCoords ) ) < tEpsilon );
        }

        // evaluate detJ
        real tDetJ = tSideDetJ * tTimeDetJ;

        // add contribution to the surface
        tSurface = tSurface + tDetJ * tIntegWeights( iGP );
    }

    // check integration points coords
    REQUIRE( tIntegPointCheck );

    // check the surface value
    Matrix< DDRMat > tSurfaceExact = { { 1.0 / 2.0, std::sqrt( 3.0 ) / 2.0, 1.0 / 2.0, 1.0 / 2.0 } };

    tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - tSurfaceExact( iSide ) ) < tEpsilon );

    REQUIRE( tSurfaceCheck );
}

// clean up
delete tSideSpaceInterp;
}

TEST_CASE( "Interpolation mesh HEX8 - Integration mesh HEX8 ",
        "[moris],[fem],[IntegMeshHEX8]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // define an interpolation mesh
    //------------------------------------------------------------------------------
    // define a HEX8 space element, i.e. space coordinates xHat
    Matrix< DDRMat > tXHatIP = {
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
    Matrix< DDRMat > tTHatIP = { { 0.0 }, { 1.0 }, { 0.5 } };

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

    // define an integration mesh
    //------------------------------------------------------------------------------
    // integration mesh geometry type
    mtk::Geometry_Type tSideGeoTypeIG = mtk::Geometry_Type::QUAD;

    // define a HEX8 integration element, i.e. space param coordinates xiHat
    Matrix< DDRMat > tXiHatIG = {
        { -1.0, -1.0, -1.0 },
        { 1.0, -1.0, -1.0 },
        { 1.0, 1.0, -1.0 },
        { -1.0, 1.0, -1.0 },
        { -1.0, -1.0, 1.0 },
        { 1.0, -1.0, 1.0 },
        { 1.0, 1.0, 1.0 },
        { -1.0, 1.0, 1.0 }
    };

    // the HEX8 integration element in space physical coordinates xHat
    Matrix< DDRMat > tXHatIG = {
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 1.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
        { 1.0, 0.0, 1.0 },
        { 1.0, 1.0, 1.0 },
        { 0.0, 1.0, 1.0 }
    };

    Matrix< DDRMat > tTHatIG = { { 0.0 }, { 1.0 } };

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule(
            tSideGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp =
            tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule
    mtk::Integration_Rule tSideIntegRule(
            tSideGeoTypeIG,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint tNumOfIntegPoints = tSideIntegrator.get_number_of_points();

    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );

    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // get side phys and param coords
    Cell< moris_index > tAllSideOrdinals = { 0, 1, 2, 3, 4, 5 };

    // boolean for surface check
    bool tSurfaceCheck    = true;
    bool tIntegPointCheck = true;

    for ( uint iSide = 0; iSide < 6; iSide++ )
    {
        // get the treated side ordinal
        uint tSideOrdinal = tAllSideOrdinals( iSide );

        // get the node ids associated to the side ordinal
        Matrix< DDSMat > tElementNodes = { //
            { 0, 1, 5, 4 },
            { 1, 2, 6, 5 },
            { 2, 3, 7, 6 },
            { 0, 4, 7, 3 },
            { 0, 3, 2, 1 },
            { 4, 5, 6, 7 }
        };
    };

    Matrix< DDSMat > tSideNodes = tElementNodes.get_row( tSideOrdinal );

    // phys coords
    Matrix< DDRMat > tSidePhysCoords( tSideNodes.numel(), 3 );
    Matrix< DDRMat > tSideParamCoords( tSideNodes.numel(), 3 );

    for ( uint iNode = 0; iNode < tSideNodes.numel(); iNode++ )
    {
        tSidePhysCoords.get_row( iNode )  = tXHatIG.get_row( tSideNodes( iNode ) );
        tSideParamCoords.get_row( iNode ) = tXiHatIG.get_row( tSideNodes( iNode ) );
    }

    // get the side normal
    Matrix< DDRMat > tTangent1 = tSidePhysCoords.get_row( 1 ) - tSidePhysCoords.get_row( 0 );
    Matrix< DDRMat > tTangent2 = tSidePhysCoords.get_row( 3 ) - tSidePhysCoords.get_row( 0 );
    Matrix< DDRMat > tNormal   = cross( tTangent1, tTangent2 );
    tNormal                    = tNormal / norm( tNormal );

    // get the side detJ
    real tSideDetJ = norm( cross( tTangent1, tTangent2 ) ) / 4.0;

    // get the time detJ
    real tTimeDetJ = ( tTHatIG( 1 ) - tTHatIG( 0 ) ) / 2.0;

    // init the surface of the integration mesh
    real tSurface = 0;

    // loop over the integration points
    for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
    {
        // get the treated integration point location in integration space
        Matrix< DDRMat > tIntegPointI = tIntegPoints.get_column( iGP );

        // get the integration point in the IP parametric space
        // via lookup tables
        Matrix< DDRMat > tRefIntegPointI( 4, 1, 0.0 );

        switch ( tSideOrdinal )
        {
            case 0:
                tRefIntegPointI = { { tIntegPointI( 0 ) }, { -1.0 }, { tIntegPointI( 1 ) }, { tIntegPointI( 2 ) } };
                break;
            case 1:
                tRefIntegPointI = { { 1.0 }, { tIntegPointI( 0 ) }, { tIntegPointI( 1 ) }, { tIntegPointI( 2 ) } };
                break;
            case 2:
                tRefIntegPointI = { { -tIntegPointI( 0 ) }, { 1.0 }, { tIntegPointI( 1 ) }, { tIntegPointI( 2 ) } };
                break;
            case 3:
                tRefIntegPointI = { { -1.0 }, { tIntegPointI( 1 ) }, { tIntegPointI( 0 ) }, { tIntegPointI( 2 ) } };
                break;
            case 4:
                tRefIntegPointI = { { tIntegPointI( 1 ) }, { tIntegPointI( 0 ) }, { -1.0 }, { tIntegPointI( 2 ) } };
                break;
            case 5:
                tRefIntegPointI = { { tIntegPointI( 0 ) }, { tIntegPointI( 1 ) }, { 1.0 }, { tIntegPointI( 2 ) } };
                break;
            default:
                MORIS_ASSERT( false, "Wrong side ordinal." );
                break;
        }

        // get the integration point in the IP parametric space
        // via the side shape functions
        Matrix< DDRMat > tRefIntegPointI2( 4, 1, 0.0 );
        Matrix< DDRMat > tN;

        tSideSpaceInterp->eval_N( tIntegPointI, tN );

        tRefIntegPointI2( { 0, 2 }, { 0, 0 } ) = trans( tN * tSideParamCoords );
        tRefIntegPointI2( 3 )                  = tIntegPointI( 2 );

        // check the integration point in the IP parametric space
        for ( uint iCoords = 0; iCoords < 4; iCoords++ )
        {
            tIntegPointCheck =
                    tIntegPointCheck &&    //
                    ( std::abs( tRefIntegPointI( iCoords ) - tRefIntegPointI2( iCoords ) ) < tEpsilon );
        }

        // evaluate detJ
        real tDetJ = tSideDetJ * tTimeDetJ;

        // add contribution to the surface
        tSurface = tSurface + tDetJ * tIntegWeights( iGP );
    }

    // check integration points coords
    REQUIRE( tIntegPointCheck );

    // check the surface value
    tSurfaceCheck = tSurfaceCheck && ( std::abs( tSurface - 1.0 ) < tEpsilon );

    REQUIRE( tSurfaceCheck );
}

// clean up
delete tSideSpaceInterp;
}
