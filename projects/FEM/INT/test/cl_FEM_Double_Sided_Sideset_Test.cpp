/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Double_Sided_Sideset_Test.cpp
 *
 */

#include "catch.hpp"
#include "cl_MTK_Integrator.hpp"            //MTK/sr
#include "cl_MTK_Interpolation_Rule.hpp"    //MTK/src
#include "fn_FEM_Side_Coordinate_Map.hpp"

#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "op_div.hpp"
#include "fn_dot.hpp"
#include "fn_trans.hpp"
#include "op_times.hpp"
#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "Double sided side-set QUAD ", "[moris],[fem],[DoubleSidedQUAD]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // integration mesh geometry type
    mtk::Geometry_Type tGeoTypeIG     = mtk::Geometry_Type::QUAD;
    mtk::Geometry_Type tSideGeoTypeIG = mtk::Geometry_Type::LINE;

    // mapping of nodes per side for HEX8
    Matrix< DDSMat > tSideNodes = { //
        { 0, 1 },
        { 1, 2 },
        { 2, 3 },
        { 3, 0 }
    };

    // hex param coords
    Matrix< DDRMat > tParamCoords = { //
        { -1.0, -1.0 },
        { +1.0, -1.0 },
        { +1.0, +1.0 },
        { -1.0, +1.0 }
    };

    // interpolation rule for the volume
    mtk::Interpolation_Rule tInterpIGRule( tGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space interpolation function
    mtk::Interpolation_Function_Base* tSpaceInterp =
            tInterpIGRule.create_space_interpolation_function();

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
            mtk::Integration_Order::BAR_1 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint tNumOfIntegPoints = tSideIntegrator.get_number_of_points();

    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );

    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // loop over leader side ordinals
    for ( uint iLeaderOrd = 0; iLeaderOrd < 4; iLeaderOrd++ )
    {
        // side ordinal for leader integration mesh
        moris_index tLeaderSideOrd = iLeaderOrd;

        // define an interpolation leader element in the phys space
        Matrix< DDRMat > tXHatIP_M;
        switch ( iLeaderOrd )
        {
            case ( 0 ):
                tXHatIP_M = { //
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 },
                    { 0.0, 0.0 }
                };
                break;
            case ( 1 ):
                tXHatIP_M = { //
                    { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 }
                };
                break;
            case ( 2 ):
                tXHatIP_M = { //
                    { 0.0, 1.0 },
                    { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 }
                };
                break;
            case ( 3 ):
                tXHatIP_M = { //
                    { 1.0, 1.0 },
                    { 0.0, 1.0 },
                    { 0.0, 0.0 },
                    { 1.0, 0.0 }
                };
                break;
            default:
            {
                MORIS_ERROR( false, " wrong leader side ordinal " );
                break;
            }
        }

        // define an integration leader element in the phys space
        Matrix< DDRMat > tXHatIG_M = tXHatIP_M;

        // define the leader integration element in param space
        Matrix< DDRMat > tXiHatIG_M = tParamCoords;

        // get the node ids associated to the leader side ordinal
        Matrix< DDSMat > tLeaderSideNodes = tSideNodes.get_row( tLeaderSideOrd );

        // phys coords amd parm coords in IP param space for the side
        uint             tNumSideNodes = tLeaderSideNodes.numel();
        Matrix< DDRMat > tLeaderSidePhysCoords( tNumSideNodes, 2 );
        Matrix< DDRMat > tLeaderSideParamCoords( tNumSideNodes, 2 );
        for ( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
        {
            tLeaderSidePhysCoords.get_row( iNode ) =
                    tXHatIG_M.get_row( tLeaderSideNodes( iNode ) );

            tLeaderSideParamCoords.get_row( iNode ) =
                    tXiHatIG_M.get_row( tLeaderSideNodes( iNode ) );
        }

        // loop over follower side ordinals
        for ( uint iFollowerOrd = 0; iFollowerOrd < 4; iFollowerOrd++ )
        {
            // side ordinal for leader integration mesh
            moris_index tFollowerSideOrd = iFollowerOrd;

            // define an interpolation follower element in the phys space
            Matrix< DDRMat > tXHatIP;
            switch ( iFollowerOrd )
            {
                case ( 0 ):
                    tXHatIP = { //
                        { 1.0, 1.0 },
                        { 1.0, 0.0 },
                        { 2.0, 0.0 },
                        { 2.0, 1.0 }
                    };
                    break;
                case ( 1 ):
                    tXHatIP = { //
                        { 2.0, 1.0 },
                        { 1.0, 1.0 },
                        { 1.0, 0.0 },
                        { 2.0, 0.0 }
                    };
                    break;
                case ( 2 ):
                    tXHatIP = { //
                        { 2.0, 0.0 },
                        { 2.0, 1.0 },
                        { 1.0, 1.0 },
                        { 1.0, 0.0 }
                    };
                    break;
                case ( 3 ):
                    tXHatIP = { //
                        { 1.0, 0.0 },
                        { 2.0, 0.0 },
                        { 2.0, 1.0 },
                        { 1.0, 1.0 }
                    };
                    break;
                default:
                    MORIS_ERROR( false, " wrong follower side ordinal " );
                    break;
            }

            // for a leader follower pair
            moris_index tFollowerNode = 0;

            // define an interpolation follower element in the phys space
            Matrix< DDRMat > tXHatIP_S = tXHatIP;

            // define an integration follower element in the phys space
            Matrix< DDRMat > tXHatIG_S = tXHatIP_S;

            // define the follower integration element in param space
            Matrix< DDRMat > tXiHatIG_S = tParamCoords;

            // get the node ids associated to the follower side ordinal
            Matrix< DDSMat > tFollowerSideNodes = tSideNodes.get_row( tFollowerSideOrd );

            // phys coords and param coords in IP param space for the follower side
            Matrix< DDRMat > tFollowerSidePhysCoords( tNumSideNodes, 2 );
            Matrix< DDRMat > tFollowerSideParamCoords( tNumSideNodes, 2 );

            for ( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
            {
                tFollowerSidePhysCoords.get_row( iNode ) =
                        tXHatIG_S.get_row( tFollowerSideNodes( iNode ) );

                tFollowerSideParamCoords.get_row( iNode ) =
                        tXiHatIG_S.get_row( tFollowerSideNodes( iNode ) );
            }

            // bool for integration point check
            bool tIntegPointCheck = true;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // get the treated integration point location in integration space
                Matrix< DDRMat > tLeaderIntegPointI = tIntegPoints.get_column( iGP );

                // get the treated integration point location in integration space
                Matrix< DDRMat > tFollowerIntegPointI =
                        side_coordinate_map(
                                tSideGeoTypeIG,
                                tFollowerNode,
                                tLeaderIntegPointI );

                // get the integration point in the IP parametric space
                // via the side shape functions for the leader side
                Matrix< DDRMat > tN1;
                Matrix< DDRMat > tN2;
                Matrix< DDRMat > tN3;
                Matrix< DDRMat > tN4;

                tSideSpaceInterp->eval_N( tLeaderIntegPointI( { 0, 0 }, { 0, 0 } ), tN1 );
                tSideSpaceInterp->eval_N( tFollowerIntegPointI( { 0, 0 }, { 0, 0 } ), tN2 );

                Matrix< DDRMat > tLeaderRefIntegPointI = trans( tN1 * tLeaderSideParamCoords );

                // get the integration point in the IP parametric space
                // via the rotation matrix for the follower side
                Matrix< DDRMat > tFollowerRefIntegPointI = trans( tN2 * tFollowerSideParamCoords );

                tSpaceInterp->eval_N( tLeaderRefIntegPointI, tN3 );
                tSpaceInterp->eval_N( tFollowerRefIntegPointI, tN4 );

                // to check only
                Matrix< DDRMat > tLeaderPhysIntegPointI = trans( tN3 * tXHatIP_M );

                // to check only
                Matrix< DDRMat > tFollowerPhysIntegPointI = trans( tN4 * tXHatIP_S );

                // check the integration point in the IP parametric space
                for ( uint iCoords = 0; iCoords < 2; iCoords++ )
                {
                    tIntegPointCheck =
                            tIntegPointCheck &&    //
                            ( std::abs( tLeaderPhysIntegPointI( iCoords ) - tFollowerPhysIntegPointI( iCoords ) ) < tEpsilon );
                }
                REQUIRE( tIntegPointCheck );
            }
        }
    }

    // clean up
    delete tSpaceInterp;
    delete tSideSpaceInterp;
}

TEST_CASE( "Double sided side-set TRI ", "[moris],[fem],[DoubleSidedTRI]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // integration mesh geometry type
    mtk::Geometry_Type tGeoTypeIG     = mtk::Geometry_Type::TRI;
    mtk::Geometry_Type tSideGeoTypeIG = mtk::Geometry_Type::LINE;

    // mapping of nodes per side for HEX8
    Matrix< DDSMat > tSideNodes = { { 0, 1 }, { 1, 2 }, { 2, 0 } };

    // hex param coords
    Matrix< DDRMat > tParamCoords = { { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0 } };

    // interpolation rule for the volume
    mtk::Interpolation_Rule tInterpIGRule(
            tGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space interpolation function
    mtk::Interpolation_Function_Base* tSpaceInterp =
            tInterpIGRule.create_space_interpolation_function();

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
            mtk::Integration_Order::BAR_1 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint tNumOfIntegPoints = tSideIntegrator.get_number_of_points();

    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );

    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // loop over leader side ordinals
    for ( uint iLeaderOrd = 0; iLeaderOrd < 3; iLeaderOrd++ )
    {
        // side ordinal for leader integration mesh
        moris_index tLeaderSideOrd = iLeaderOrd;

        // define an interpolation leader element in the phys space
        Matrix< DDRMat > tXHatIP_M;
        switch ( iLeaderOrd )
        {
            case ( 0 ):
                tXHatIP_M = { { 1.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 0.0 } };
                break;
            case ( 1 ):
                tXHatIP_M = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 } };
                break;
            case ( 2 ):
                tXHatIP_M = { { 0.0, 1.0 }, { 0.0, 0.0 }, { 1.0, 0.0 } };
                break;
            default:
            {
                MORIS_ERROR( false, " wrong leader side ordinal " );
                break;
            }
        }

        // define an integration leader element in the phys space
        Matrix< DDRMat > tXHatIG_M = tXHatIP_M;

        // define the leader integration element in param space
        Matrix< DDRMat > tXiHatIG_M = tParamCoords;

        // get the node ids associated to the leader side ordinal
        Matrix< DDSMat > tLeaderSideNodes = tSideNodes.get_row( tLeaderSideOrd );

        // phys coords amd parm coords in IP param space for the side
        uint tNumSideNodes = tLeaderSideNodes.numel();

        Matrix< DDRMat > tLeaderSidePhysCoords( tNumSideNodes, 2 );
        Matrix< DDRMat > tLeaderSideParamCoords( tNumSideNodes, 3 );

        for ( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
        {
            tLeaderSidePhysCoords.get_row( iNode ) =
                    tXHatIG_M.get_row( tLeaderSideNodes( iNode ) );

            tLeaderSideParamCoords.get_row( iNode ) =
                    tXiHatIG_M.get_row( tLeaderSideNodes( iNode ) );
        }

        // loop over follower side ordinals
        for ( uint iFollowerOrd = iLeaderOrd; iFollowerOrd < 3; iFollowerOrd++ )
            for ( uint iFollowerOrd = 0; iFollowerOrd < 3; iFollowerOrd++ )
            {
                // side ordinal for leader integration mesh
                moris_index tFollowerSideOrd = iFollowerOrd;

                // define an interpolation follower element in the phys space
                Matrix< DDRMat > tXHatIP;
                switch ( iFollowerOrd )
                {
                    case ( 0 ):
                        tXHatIP = { { 0.0, 1.0 }, { 1.0, 0.0 }, { 1.0, 1.0 } };
                        break;
                    case ( 1 ):
                        tXHatIP = { { 1.0, 1.0 }, { 0.0, 1.0 }, { 1.0, 0.0 } };
                        break;
                    case ( 2 ):
                        tXHatIP = { { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
                        break;
                    default:
                        MORIS_ERROR( false, " wrong follower side ordinal " );
                        break;
                }

                // for a leader follower pair
                moris_index tFollowerNode = 0;

                // define an interpolation follower element in the phys space
                Matrix< DDRMat > tXHatIP_S = tXHatIP;

                // define an integration follower element in the phys space
                Matrix< DDRMat > tXHatIG_S = tXHatIP_S;

                // define the follower integration element in param space
                Matrix< DDRMat > tXiHatIG_S = tParamCoords;

                // get the node ids associated to the follower side ordinal
                Matrix< DDSMat > tFollowerSideNodes = tSideNodes.get_row( tFollowerSideOrd );

                // phys coords amd parm coords in IP param space for the side
                Matrix< DDRMat > tFollowerSidePhysCoords( tNumSideNodes, 2 );
                Matrix< DDRMat > tFollowerSideParamCoords( tNumSideNodes, 3 );

                for ( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
                {
                    tFollowerSidePhysCoords.get_row( iNode ) =
                            tXHatIG_S.get_row( tFollowerSideNodes( iNode ) );

                    tFollowerSideParamCoords.get_row( iNode ) =
                            tXiHatIG_S.get_row( tFollowerSideNodes( iNode ) );
                }

                // bool for integration point check
                bool tIntegPointCheck = true;

                // loop over the integration points
                for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tLeaderIntegPointI = tIntegPoints.get_column( iGP );

                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tFollowerIntegPointI =
                            side_coordinate_map(
                                    tSideGeoTypeIG,
                                    tFollowerNode,
                                    tLeaderIntegPointI );

                    // get the integration point in the IP parametric space
                    // via the side shape functions for the leader side
                    Matrix< DDRMat > tLeaderRefIntegPointI( 4, 1, 0.0 );
                    Matrix< DDRMat > tN1;

                    tSideSpaceInterp->eval_N( tLeaderIntegPointI( { 0, 0 }, { 0, 0 } ), tN1 );

                    tLeaderRefIntegPointI( { 0, 2 }, { 0, 0 } ) = trans( tN1 * tLeaderSideParamCoords );
                    tLeaderRefIntegPointI( 3 )                  = tLeaderIntegPointI( 1 );

                    // get the integration point in the IP parametric space
                    // via the rotation matrix for the follower side
                    Matrix< DDRMat > tFollowerRefIntegPointI( 4, 1, 0.0 );
                    Matrix< DDRMat > tN2;

                    tSideSpaceInterp->eval_N( tFollowerIntegPointI( { 0, 0 }, { 0, 0 } ), tN2 );

                    tFollowerRefIntegPointI( { 0, 2 }, { 0, 0 } ) = trans( tN2 * tFollowerSideParamCoords );
                    tFollowerRefIntegPointI( 3 )                  = tFollowerIntegPointI( 1 );

                    // to check only
                    Matrix< DDRMat > tN3;
                    tSpaceInterp->eval_N( tLeaderRefIntegPointI, tN3 );
                    Matrix< DDRMat > tLeaderPhysIntegPointI = trans( tN3 * tXHatIP_M );

                    // to check only
                    Matrix< DDRMat > tN4;
                    tSpaceInterp->eval_N( tFollowerRefIntegPointI, tN4 );
                    Matrix< DDRMat > tFollowerPhysIntegPointI = trans( tN4 * tXHatIP_S );

                    // check the integration point in the IP parametric space
                    for ( uint iCoords = 0; iCoords < 2; iCoords++ )
                    {
                        tIntegPointCheck =
                                tIntegPointCheck &&    //
                                ( std::abs( tLeaderPhysIntegPointI( iCoords ) - tFollowerPhysIntegPointI( iCoords ) ) < tEpsilon );
                    }
                    REQUIRE( tIntegPointCheck );
                }
            }
    }

    // clean up
    delete tSpaceInterp;
    delete tSideSpaceInterp;
}

TEST_CASE( "Double sided side-set TET ", "[moris],[fem],[DoubleSidedTET]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // integration mesh geometry type
    mtk::Geometry_Type tGeoTypeIG     = mtk::Geometry_Type::TET;
    mtk::Geometry_Type tSideGeoTypeIG = mtk::Geometry_Type::TRI;

    // mapping of nodes per side for HEX8
    Matrix< DDSMat > tSideNodes = { //
        { 0, 1, 3 },
        { 1, 2, 3 },
        { 0, 3, 2 },
        { 0, 2, 1 }
    };
    Matrix< DDSMat > tFacingSideNodes = { //
        { 2 },
        { 0 },
        { 1 },
        { 3 }
    };

    // possible permutation cases
    Matrix< IndexMat > tPermuteCases = { //
        { 0, 1, 2 },
        { 1, 2, 0 },
        { 2, 0, 1 }
    };

    // tet param coords
    Matrix< DDRMat > tParamCoords = {
        { 1.0, 0.0, 0.0 },    //
        { 0.0, 1.0, 0.0 },    //
        { 0.0, 0.0, 0.0 },    //
        { 0.0, 0.0, 1.0 }     //
    };

    // interpolation rule for the volume
    mtk::Interpolation_Rule tInterpIGRule(
            tGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space interpolation function
    mtk::Interpolation_Function_Base* tSpaceInterp =    //
            tInterpIGRule.create_space_interpolation_function();

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule(
            tSideGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp =    //
            tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule for the side
    mtk::Integration_Rule tSideIntegRule(
            tSideGeoTypeIG,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TRI_3,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points (space & time) and weights on side
    uint tNumOfIntegPoints = tSideIntegrator.get_number_of_points();

    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );

    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // loop over leader side ordinals
    for ( uint iLeaderOrd = 0; iLeaderOrd < 4; iLeaderOrd++ )
    {
        // side ordinal for leader integration mesh
        moris_index tLeaderSideOrd = iLeaderOrd;

        // define an interpolation leader element in the phys space
        Matrix< DDRMat > tXHatIP_M;

        switch ( tLeaderSideOrd )
        {
            case ( 0 ):
            {
                tXHatIP_M = { //
                    { 1.0, 1.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 1.0 }
                };
                break;
            }
            case ( 1 ):
            {
                tXHatIP_M = { //
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 0.0, 1.0 }
                };
                break;
            }
            case ( 2 ):
            {
                tXHatIP_M = { //
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 0.0, 0.0 }
                };
                break;
            }
            case ( 3 ):
            {
                tXHatIP_M = { //
                    { 1.0, 1.0, 0.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0 }
                };
                break;
            }
            default:
            {
                MORIS_ERROR( false, " wrong leader side ordinal " );
                break;
            }
        }

        // define an integration leader element in the phys space
        Matrix< DDRMat > tXHatIG_M = tXHatIP_M;

        // define the leader integration element in param space
        Matrix< DDRMat > tXiHatIG_M = tParamCoords;

        // get the node ids associated to the leader side ordinal
        Matrix< DDSMat > tLeaderSideNodes = tSideNodes.get_row( tLeaderSideOrd );

        // phys coords amd parm coords in IP param space for the side
        uint tNumSideNodes = tLeaderSideNodes.numel();

        Matrix< DDRMat > tLeaderSidePhysCoords( tNumSideNodes, 3 );
        Matrix< DDRMat > tLeaderSideParamCoords( tNumSideNodes, 3 );

        for ( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
        {
            tLeaderSidePhysCoords.get_row( iNode ) =    //
                    tXHatIG_M.get_row( tLeaderSideNodes( iNode ) );

            tLeaderSideParamCoords.get_row( iNode ) =    //
                    tXiHatIG_M.get_row( tLeaderSideNodes( iNode ) );
        }

        // loop over follower side ordinals
        for ( uint iFollowerOrd = 0; iFollowerOrd < 4; iFollowerOrd++ )
        {
            // side ordinal for leader integration mesh
            moris_index tFollowerSideOrd = iFollowerOrd;

            // sort the nodes on side or not
            Matrix< DDSMat > tNodesOnFollowerSide    = tSideNodes.get_row( iFollowerOrd );
            Matrix< DDSMat > tNodesNotOnFollowerSide = tFacingSideNodes.get_row( iFollowerOrd );

            // define an interpolation follower element in the phys space
            Matrix< DDRMat > tXHatIP;

            switch ( tFollowerSideOrd )
            {
                case ( 0 ):
                {
                    tXHatIP = { //
                        { 1.0, 1.0, 0.0 },
                        { 1.0, 0.0, 1.0 },
                        { 2.0, 0.0, 0.0 },
                        { 1.0, 0.0, 0.0 }
                    };
                    break;
                }
                case ( 1 ):
                {
                    tXHatIP = { //
                        { 2.0, 0.0, 0.0 },
                        { 1.0, 1.0, 0.0 },
                        { 1.0, 0.0, 1.0 },
                        { 1.0, 0.0, 0.0 }
                    };
                    break;
                }
                case ( 2 ):
                {
                    tXHatIP = { //
                        { 1.0, 1.0, 0.0 },
                        { 2.0, 0.0, 0.0 },
                        { 1.0, 0.0, 0.0 },
                        { 1.0, 0.0, 1.0 }
                    };
                    break;
                }
                case ( 3 ):
                {
                    tXHatIP = { //
                        { 1.0, 1.0, 0.0 },
                        { 1.0, 0.0, 0.0 },
                        { 1.0, 0.0, 1.0 },
                        { 2.0, 0.0, 0.0 }
                    };
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, " wrong follower side ordinal " );
                    break;
                }
            }

            // node permutation for $ rotation in plane
            for ( uint iPermute = 0; iPermute < tPermuteCases.n_rows(); iPermute++ )
            {
                // get the first corresponding node on the follower side
                moris_index tFollowerNode = iPermute;

                // for a leader follower pair

                // define an interpolation follower element in the phys space
                Matrix< DDRMat > tXHatIP_S( 4, 3, 0.0 );

                // fill the phys coords
                for ( uint i = 0; i < 3; i++ )
                {
                    moris_index tNodeSide    = tNodesOnFollowerSide( tPermuteCases( iPermute, i ) );
                    moris_index tNodeNotSide = tNodesNotOnFollowerSide( 0 );

                    // place the interface nodes
                    tXHatIP_S.get_row( tNodeSide )    = tXHatIP.get_row( tNodesOnFollowerSide( i ) );
                    tXHatIP_S.get_row( tNodeNotSide ) = tXHatIP.get_row( tNodesNotOnFollowerSide( 0 ) );
                }

                // define an integration follower element in the phys space
                Matrix< DDRMat > tXHatIG_S = tXHatIP_S;

                // define the follower integration element in param space
                Matrix< DDRMat > tXiHatIG_S = tParamCoords;

                // get the node ids associated to the leader side ordinal
                Matrix< DDSMat > tFollowerSideNodes = tSideNodes.get_row( tFollowerSideOrd );

                // phys coords amd parm coords in IP param space for the side
                Matrix< DDRMat > tFollowerSidePhysCoords( tNumSideNodes, 3 );
                Matrix< DDRMat > tFollowerSideParamCoords( tNumSideNodes, 3 );

                for ( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
                {
                    tFollowerSidePhysCoords.get_row( iNode )  = tXHatIG_S.get_row( tFollowerSideNodes( iNode ) );
                    tFollowerSideParamCoords.get_row( iNode ) = tXiHatIG_S.get_row( tFollowerSideNodes( iNode ) );
                }

                // bool for integration point check
                bool tIntegPointCheck = true;

                // loop over the integration points
                for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tLeaderIntegPointI = tIntegPoints.get_column( iGP );

                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tFollowerIntegPointI =
                            side_coordinate_map(
                                    tSideGeoTypeIG,
                                    tFollowerNode,
                                    tLeaderIntegPointI );

                    // get the integration point in the IP parametric space
                    // via the side shape functions for the leader side
                    Matrix< DDRMat > tN1;
                    tSideSpaceInterp->eval_N( tLeaderIntegPointI( { 0, 1 }, { 0, 0 } ), tN1 );
                    Matrix< DDRMat > tLeaderRefIntegPointI = trans( tN1 * tLeaderSideParamCoords );

                    // get the integration point in the IP parametric space
                    // via the rotation matrix for the follower side
                    Matrix< DDRMat > tN2;
                    tSideSpaceInterp->eval_N( tFollowerIntegPointI( { 0, 1 }, { 0, 0 } ), tN2 );
                    Matrix< DDRMat > tFollowerRefIntegPointI = trans( tN2 * tFollowerSideParamCoords );

                    // to check only
                    Matrix< DDRMat > tN3;
                    tSpaceInterp->eval_N( tLeaderRefIntegPointI, tN3 );
                    Matrix< DDRMat > tLeaderPhysIntegPointI = trans( tN3 * tXHatIP_M );

                    // to check only
                    Matrix< DDRMat > tN4;
                    tSpaceInterp->eval_N( tFollowerRefIntegPointI, tN4 );
                    Matrix< DDRMat > tFollowerPhysIntegPointI = trans( tN4 * tXHatIP_S );

                    // check the integration point in the IP parametric space
                    for ( uint iCoords = 0; iCoords < 3; iCoords++ )
                    {
                        tIntegPointCheck =
                                tIntegPointCheck &&    //
                                ( std::abs( tLeaderPhysIntegPointI( iCoords ) - tFollowerPhysIntegPointI( iCoords ) ) < tEpsilon );
                    }
                    REQUIRE( tIntegPointCheck );
                }
            }
        }
    }
    // clean up
    delete tSpaceInterp;
    delete tSideSpaceInterp;
}

TEST_CASE( "Double sided side-set HEX ", "[moris],[fem],[DoubleSidedHEX]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-12;

    // integration mesh geometry type
    mtk::Geometry_Type tGeoTypeIG     = mtk::Geometry_Type::HEX;
    mtk::Geometry_Type tSideGeoTypeIG = mtk::Geometry_Type::QUAD;

    // mapping of nodes per side for HEX8
    Matrix< DDSMat > tSideNodes = { //
        { 0, 1, 5, 4 },
        { 1, 2, 6, 5 },
        { 2, 3, 7, 6 },
        { 0, 4, 7, 3 },
        { 0, 3, 2, 1 },
        { 4, 5, 6, 7 }
    };
    Matrix< DDSMat > tFacingSideNodes = { //
        { 3, 2, 6, 7 },
        { 0, 3, 7, 4 },
        { 1, 0, 4, 5 },
        { 1, 5, 6, 2 },
        { 4, 7, 6, 5 },
        { 0, 1, 2, 3 }
    };

    // possible permutation cases
    Matrix< IndexMat > tPermuteCases = { //
        { 0, 1, 2, 3 },
        { 1, 2, 3, 0 },
        { 2, 3, 0, 1 },
        { 3, 0, 1, 2 }
    };

    // hex param coords
    Matrix< DDRMat > tHexParamCoords = { //
        { -1.0, -1.0, -1.0 },
        { +1.0, -1.0, -1.0 },
        { +1.0, +1.0, -1.0 },
        { -1.0, +1.0, -1.0 },
        { -1.0, -1.0, +1.0 },
        { +1.0, -1.0, +1.0 },
        { +1.0, +1.0, +1.0 },
        { -1.0, +1.0, +1.0 }
    };

    // interpolation rule for the volume
    mtk::Interpolation_Rule tInterpIGRule(
            tGeoTypeIG,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space interpolation function
    mtk::Interpolation_Function_Base* tSpaceInterp =    //
            tInterpIGRule.create_space_interpolation_function();

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
            mtk::Integration_Order::QUAD_5x5,
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

    // loop over leader side ordinals
    for ( uint iLeaderOrd = 0; iLeaderOrd < 6; iLeaderOrd++ )
    {
        // side ordinal for leader integration mesh
        moris_index tLeaderSideOrd = iLeaderOrd;

        // define an interpolation leader element in the phys space
        Matrix< DDRMat > tXHatIP_M;
        switch ( iLeaderOrd )
        {
            case ( 0 ):
                tXHatIP_M = { //
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 },
                    { 0.0, 0.0, 1.0 }
                };
                break;
            case ( 1 ):
                tXHatIP_M = { //
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 }
                };
                break;
            case ( 2 ):
                tXHatIP_M = { //
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 1.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 }
                };
                break;
            case ( 3 ):
                tXHatIP_M = { //
                    { 1.0, 0.0, 0.0 },
                    { 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 1.0, 1.0 },
                    { 1.0, 1.0, 1.0 }
                };
                break;
            case ( 4 ):
                tXHatIP_M = { //
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 0.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 0.0, 1.0, 1.0 },
                    { 0.0, 1.0, 0.0 }
                };
                break;
            case ( 5 ):
                tXHatIP_M = { //
                    { 0.0, 0.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 1.0, 1.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 1.0, 1.0, 1.0 },
                    { 1.0, 0.0, 1.0 }
                };
                break;
            default:
            {
                MORIS_ERROR( false, " wrong leader side ordinal " );
                break;
            }
        }

        // define an integration leader element in the phys space
        Matrix< DDRMat > tXHatIG_M = tXHatIP_M;

        // define the leader integration element in param space
        Matrix< DDRMat > tXiHatIG_M = tHexParamCoords;

        // get the node ids associated to the leader side ordinal
        Matrix< DDSMat > tLeaderSideNodes = tSideNodes.get_row( tLeaderSideOrd );

        // phys coords amd parm coords in IP param space for the side
        uint tNumSideNodes = tLeaderSideNodes.numel();

        Matrix< DDRMat > tLeaderSidePhysCoords( tNumSideNodes, 3 );
        Matrix< DDRMat > tLeaderSideParamCoords( tNumSideNodes, 3 );

        for ( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
        {
            tLeaderSidePhysCoords.get_row( iNode )  = tXHatIG_M.get_row( tLeaderSideNodes( iNode ) );
            tLeaderSideParamCoords.get_row( iNode ) = tXiHatIG_M.get_row( tLeaderSideNodes( iNode ) );
        }

        // loop over follower side ordinals
        // for ( uint iFollowerOrd = iLeaderOrd; iFollowerOrd < 6; iFollowerOrd++ )
        for ( uint iFollowerOrd = 0; iFollowerOrd < 6; iFollowerOrd++ )
        {
            // side ordinal for leader integration mesh
            moris_index tFollowerSideOrd = iFollowerOrd;

            // sort the nodes on side or not
            Matrix< DDSMat > tNodesOnFollowerSide   = tSideNodes.get_row( iFollowerOrd );
            Matrix< DDSMat > tNodeNotOnFollowerSide = tFacingSideNodes.get_row( iFollowerOrd );

            // define an interpolation follower element in the phys space
            Matrix< DDRMat > tXHatIP;
            switch ( iFollowerOrd )
            {
                case ( 0 ):
                    tXHatIP = { //
                        { 1.0, 0.0, 0.0 },
                        { 1.0, 0.0, 1.0 },
                        { 2.0, 0.0, 1.0 },
                        { 2.0, 0.0, 0.0 },
                        { 1.0, 1.0, 0.0 },
                        { 1.0, 1.0, 1.0 },
                        { 2.0, 1.0, 1.0 },
                        { 2.0, 1.0, 0.0 }
                    };
                    break;
                case ( 1 ):
                    tXHatIP = { //
                        { 2.0, 0.0, 0.0 },
                        { 1.0, 0.0, 0.0 },
                        { 1.0, 0.0, 1.0 },
                        { 2.0, 0.0, 1.0 },
                        { 2.0, 1.0, 0.0 },
                        { 1.0, 1.0, 0.0 },
                        { 1.0, 1.0, 1.0 },
                        { 2.0, 1.0, 1.0 }
                    };
                    break;
                case ( 2 ):
                    tXHatIP = { //
                        { 2.0, 0.0, 1.0 },
                        { 2.0, 0.0, 0.0 },
                        { 1.0, 0.0, 0.0 },
                        { 1.0, 0.0, 1.0 },
                        { 2.0, 1.0, 1.0 },
                        { 2.0, 1.0, 0.0 },
                        { 1.0, 1.0, 0.0 },
                        { 1.0, 1.0, 1.0 }
                    };
                    break;
                case ( 3 ):
                    tXHatIP = { //
                        { 1.0, 0.0, 0.0 },
                        { 2.0, 0.0, 0.0 },
                        { 2.0, 1.0, 0.0 },
                        { 1.0, 1.0, 0.0 },
                        { 1.0, 0.0, 1.0 },
                        { 2.0, 0.0, 1.0 },
                        { 2.0, 1.0, 1.0 },
                        { 1.0, 1.0, 1.0 }
                    };
                    break;
                case ( 4 ):
                    tXHatIP = { //
                        { 1.0, 0.0, 0.0 },
                        { 1.0, 1.0, 0.0 },
                        { 1.0, 1.0, 1.0 },
                        { 1.0, 0.0, 1.0 },
                        { 2.0, 0.0, 0.0 },
                        { 2.0, 1.0, 0.0 },
                        { 2.0, 1.0, 1.0 },
                        { 2.0, 0.0, 1.0 }
                    };
                    break;
                case ( 5 ):
                    tXHatIP = { //
                        { 2.0, 0.0, 0.0 },
                        { 2.0, 0.0, 1.0 },
                        { 2.0, 1.0, 1.0 },
                        { 2.0, 1.0, 0.0 },
                        { 1.0, 0.0, 0.0 },
                        { 1.0, 0.0, 1.0 },
                        { 1.0, 1.0, 1.0 },
                        { 1.0, 1.0, 0.0 }
                    };
                    break;
                default:
                    MORIS_ERROR( false, " wrong follower side ordinal " );
                    break;
            }

            // node permutation for $ rotation in plane
            for ( uint iPermute = 0; iPermute < tPermuteCases.n_rows(); iPermute++ )
            {
                // get the first corresponding node on the follower side
                moris_index tFollowerNode = iPermute;

                // for a leader follower pair

                // define an interpolation follower element in the phys space
                Matrix< DDRMat > tXHatIP_S( 8, 3, 0.0 );

                // fill the phys coords
                for ( uint i = 0; i < 4; i++ )
                {
                    moris_index tNodeSide    = tNodesOnFollowerSide( tPermuteCases( iPermute, i ) );
                    moris_index tNodeNotSide = tNodeNotOnFollowerSide( tPermuteCases( iPermute, i ) );

                    // place the interface nodes
                    tXHatIP_S.get_row( tNodeSide )    = tXHatIP.get_row( tNodesOnFollowerSide( i ) );
                    tXHatIP_S.get_row( tNodeNotSide ) = tXHatIP.get_row( tNodeNotOnFollowerSide( i ) );
                }

                // define an integration follower element in the phys space
                Matrix< DDRMat > tXHatIG_S = tXHatIP_S;

                // define the follower integration element in param space
                Matrix< DDRMat > tXiHatIG_S = tHexParamCoords;

                // get the node ids associated to the leader side ordinal
                Matrix< DDSMat > tFollowerSideNodes = tSideNodes.get_row( tFollowerSideOrd );

                // phys coords and parm coords in IP param space for the side
                Matrix< DDRMat > tFollowerSidePhysCoords( tNumSideNodes, 3 );
                Matrix< DDRMat > tFollowerSideParamCoords( tNumSideNodes, 3 );

                for ( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
                {
                    tFollowerSidePhysCoords.get_row( iNode )  = tXHatIG_S.get_row( tFollowerSideNodes( iNode ) );
                    tFollowerSideParamCoords.get_row( iNode ) = tXiHatIG_S.get_row( tFollowerSideNodes( iNode ) );
                }

                // bool for integration point check
                bool tIntegPointCheck = true;

                // loop over the integration points
                for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tLeaderIntegPointI = tIntegPoints.get_column( iGP );

                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tFollowerIntegPointI =
                            side_coordinate_map(
                                    tSideGeoTypeIG,
                                    tFollowerNode,
                                    tLeaderIntegPointI );

                    // get the integration point in the IP parametric space
                    // via the side shape functions for the leader side
                    Matrix< DDRMat > tN1;
                    tSideSpaceInterp->eval_N( tLeaderIntegPointI( { 0, 1 }, { 0, 0 } ), tN1 );
                    Matrix< DDRMat > tLeaderRefIntegPointI = trans( tN1 * tLeaderSideParamCoords );

                    // get the integration point in the IP parametric space
                    // via the rotation matrix for the follower side
                    Matrix< DDRMat > tN2;
                    tSideSpaceInterp->eval_N( tFollowerIntegPointI( { 0, 1 }, { 0, 0 } ), tN2 );
                    Matrix< DDRMat > tFollowerRefIntegPointI = trans( tN2 * tFollowerSideParamCoords );

                    // to check only
                    Matrix< DDRMat > tN3;
                    tSpaceInterp->eval_N( tLeaderRefIntegPointI, tN3 );
                    Matrix< DDRMat > tLeaderPhysIntegPointI = trans( tN3 * tXHatIP_M );

                    // to check only
                    Matrix< DDRMat > tN4;
                    tSpaceInterp->eval_N( tFollowerRefIntegPointI, tN4 );
                    Matrix< DDRMat > tFollowerPhysIntegPointI = trans( tN4 * tXHatIP_S );

                    // check the integration point in the IP parametric space
                    for ( uint iCoords = 0; iCoords < 3; iCoords++ )
                    {
                        tIntegPointCheck = tIntegPointCheck && ( std::abs( tLeaderPhysIntegPointI( iCoords ) - tFollowerPhysIntegPointI( iCoords ) ) < tEpsilon );
                    }
                    REQUIRE( tIntegPointCheck );
                }
            }
        }
    }
    // clean up
    delete tSpaceInterp;
    delete tSideSpaceInterp;
}
