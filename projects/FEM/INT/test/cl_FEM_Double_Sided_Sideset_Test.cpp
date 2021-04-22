#include "catch.hpp"
#include "cl_MTK_Integrator.hpp" //MTK/sr
#include "cl_MTK_Interpolation_Rule.hpp" //MTK/src
#include "fn_FEM_Rotation_Matrix.hpp"

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
    Matrix< DDSMat > tSideNodes = {{ 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }};

    // hex param coords
    Matrix< DDRMat > tParamCoords = {{ -1.0, -1.0 }, {  1.0, -1.0 },
                                     {  1.0,  1.0 }, { -1.0,  1.0 }};

    // interpolation rule for the volume
    mtk::Interpolation_Rule tInterpIGRule( tGeoTypeIG,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR );

    // create a space interpolation function
    mtk::Interpolation_Function_Base* tSpaceInterp = tInterpIGRule.create_space_interpolation_function();

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule( tSideGeoTypeIG,
                                          mtk::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR,
                                          mtk::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp = tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule for the side
    mtk::Integration_Rule tSideIntegRule( tSideGeoTypeIG,
                                     mtk::Integration_Type::GAUSS,
                                     mtk::Integration_Order::BAR_2,
                                     mtk::Integration_Type::GAUSS,
                                     mtk::Integration_Order::BAR_1 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );
    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // loop over master side ordinals
    for ( uint iMasterOrd = 0; iMasterOrd < 4; iMasterOrd++ )
    {
        // side ordinal for master integration mesh
        moris_index tMasterSideOrd = iMasterOrd;

        // define an interpolation master element in the phys space
        Matrix< DDRMat > tXHatIP_M;
        switch ( iMasterOrd )
        {
            case( 0 ):
                tXHatIP_M = {{ 1.0, 0.0 }, { 1.0, 1.0 },
                             { 0.0, 1.0 }, { 0.0, 0.0 }};
                break;
            case( 1 ):
                tXHatIP_M = {{ 0.0, 0.0 }, { 1.0, 0.0 },
                             { 1.0, 1.0 }, { 0.0, 1.0 }};
                break;
            case( 2 ):
                tXHatIP_M = {{ 0.0, 1.0 }, { 0.0, 0.0 },
                             { 1.0, 0.0 }, { 1.0, 1.0 }};
                break;
            case( 3 ):
                tXHatIP_M = {{ 1.0, 1.0 }, { 0.0, 1.0 },
                             { 0.0, 0.0 }, { 1.0, 0.0 }};
                break;
            default:
            {
                MORIS_ERROR( false, " wrong master side ordinal ");
                break;
            }
        }

        // define an integration master element in the phys space
        Matrix< DDRMat > tXHatIG_M = tXHatIP_M;

        // define the master integration element in param space
        Matrix< DDRMat > tXiHatIG_M = tParamCoords;

        // get the node ids associated to the master side ordinal
        Matrix< DDSMat > tMasterSideNodes = tSideNodes.get_row( tMasterSideOrd );

        // phys coords amd parm coords in IP param space for the side
        uint tNumSideNodes = tMasterSideNodes.numel();
        Matrix< DDRMat > tMasterSidePhysCoords( tNumSideNodes, 2 );
        Matrix< DDRMat > tMasterSideParamCoords( tNumSideNodes, 2 );
        for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
        {
            tMasterSidePhysCoords.get_row( iNode )  = tXHatIG_M.get_row( tMasterSideNodes( iNode ) );
            tMasterSideParamCoords.get_row( iNode ) = tXiHatIG_M.get_row( tMasterSideNodes( iNode ) );
        }

        // loop over slave side ordinals
        for ( uint iSlaveOrd = 0; iSlaveOrd < 4; iSlaveOrd++ )
        {
            // side ordinal for master integration mesh
            moris_index tSlaveSideOrd = iSlaveOrd;

            // define an interpolation slave element in the phys space
            Matrix< DDRMat > tXHatIP;
            switch ( iSlaveOrd )
            {
                case( 0 ):
                    tXHatIP = {{ 1.0, 1.0 }, { 1.0, 0.0 },
                               { 2.0, 0.0 }, { 2.0, 1.0 }};
                    break;
                case( 1 ):
                    tXHatIP = {{ 2.0, 1.0 }, { 1.0, 1.0 },
                               { 1.0, 0.0 }, { 2.0, 0.0 }};
                    break;
                case( 2 ):
                    tXHatIP = {{ 2.0, 0.0 }, { 2.0, 1.0 },
                               { 1.0, 1.0 }, { 1.0, 0.0 }};
                    break;
                case( 3 ):
                    tXHatIP = {{ 1.0, 0.0 }, { 2.0, 0.0 },
                               { 2.0, 1.0 }, { 1.0, 1.0 }};
                    break;
                default:
                    MORIS_ERROR( false, " wrong slave side ordinal ");
                    break;
            }

            // for a master slave pair
            moris_index tSlaveNode = 0;
            //std::cout<<"A master slave pair + slave node"<<" "<<tMasterSideOrd<<" - "<<tSlaveSideOrd<<" - "<<tSlaveNode<<"-----------------------"<<std::endl;

            // get the rotation matrix
            Matrix< DDRMat > tR;
            rotation_matrix( tSideGeoTypeIG, tSlaveNode, tR );

            // define an interpolation slave element in the phys space
            Matrix< DDRMat > tXHatIP_S = tXHatIP;

            // define an integration slave element in the phys space
            Matrix< DDRMat > tXHatIG_S = tXHatIP_S;

            // define the slave integration element in param space
            Matrix< DDRMat > tXiHatIG_S = tParamCoords;

            // get the node ids associated to the slave side ordinal
            Matrix< DDSMat > tSlaveSideNodes = tSideNodes.get_row( tSlaveSideOrd );

            // phys coords and param coords in IP param space for the slave side
            Matrix< DDRMat > tSlaveSidePhysCoords( tNumSideNodes, 2 );
            Matrix< DDRMat > tSlaveSideParamCoords( tNumSideNodes, 2 );
            for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
            {
                tSlaveSidePhysCoords.get_row( iNode )  = tXHatIG_S.get_row( tSlaveSideNodes( iNode ) );
                tSlaveSideParamCoords.get_row( iNode ) = tXiHatIG_S.get_row( tSlaveSideNodes( iNode ) );
            }

            // bool for integration point check
            bool tIntegPointCheck = true;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // get the treated integration point location in integration space
                Matrix< DDRMat > tMasterIntegPointI = tIntegPoints.get_column( iGP );

                // get the treated integration point location in integration space
                Matrix< DDRMat > tSlaveIntegPointI = tMasterIntegPointI;
                tSlaveIntegPointI({0,0}, {0,0}) = tR * tMasterIntegPointI({0,0}, {0,0});

                // get the integration point in the IP parametric space
                // via the side shape functions for the master side
                Matrix< DDRMat > tN1;
                Matrix< DDRMat > tN2;
                Matrix< DDRMat > tN3;
                Matrix< DDRMat > tN4;
                tSideSpaceInterp->eval_N( tMasterIntegPointI({0,0}, {0,0}), tN1 );
                tSideSpaceInterp->eval_N( tSlaveIntegPointI({0,0}, {0,0}), tN2 );

                Matrix< DDRMat > tMasterRefIntegPointI = trans( tN1  * tMasterSideParamCoords );

                // get the integration point in the IP parametric space
                // via the rotation matrix for the slave side
                Matrix< DDRMat > tSlaveRefIntegPointI = trans( tN2 * tSlaveSideParamCoords );

                tSpaceInterp->eval_N( tMasterRefIntegPointI, tN3 );
                tSpaceInterp->eval_N( tSlaveRefIntegPointI, tN4 );
                // to check only
                Matrix< DDRMat > tMasterPhysIntegPointI = trans( tN3 * tXHatIP_M );

                // to check only
                Matrix< DDRMat > tSlavePhysIntegPointI = trans( tN4 * tXHatIP_S );

                // check the integration point in the IP parametric space
                for( uint iCoords = 0; iCoords < 2; iCoords++ )
                {
                    tIntegPointCheck = tIntegPointCheck && ( std::abs( tMasterPhysIntegPointI( iCoords ) - tSlavePhysIntegPointI( iCoords ) ) < tEpsilon );
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
    Matrix< DDSMat > tSideNodes = {{ 0, 1 }, { 1, 2 }, { 2, 0 }};

    // hex param coords
    Matrix< DDRMat > tParamCoords = {{ 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, {  0.0, 0.0, 1.0 }};

    // interpolation rule for the volume
    mtk::Interpolation_Rule tInterpIGRule( tGeoTypeIG,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR );

    // create a space interpolation function
    mtk::Interpolation_Function_Base* tSpaceInterp = tInterpIGRule.create_space_interpolation_function();

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule( tSideGeoTypeIG,
                                          mtk::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR,
                                          mtk::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp = tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule for the side
    mtk::Integration_Rule tSideIntegRule( tSideGeoTypeIG,
                                     mtk::Integration_Type::GAUSS,
                                     mtk::Integration_Order::BAR_2,
                                     mtk::Integration_Type::GAUSS,
                                     mtk::Integration_Order::BAR_1 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );
    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // loop over master side ordinals
    for ( uint iMasterOrd = 0; iMasterOrd < 3; iMasterOrd++ )
    {
        // side ordinal for master integration mesh
        moris_index tMasterSideOrd = iMasterOrd;

        // define an interpolation master element in the phys space
        Matrix< DDRMat > tXHatIP_M;
        switch ( iMasterOrd )
        {
            case( 0 ):
                tXHatIP_M = {{ 1.0, 0.0 }, { 0.0, 1.0 }, { 0.0, 0.0 }};
                break;
            case( 1 ):
                tXHatIP_M = {{ 0.0,  0.0 }, { 1.0, 0.0 }, { 0.0, 1.0 }};
                break;
            case( 2 ):
                tXHatIP_M = {{ 0.0, 1.0 }, { 0.0, 0.0 }, { 1.0, 0.0 }};
                break;
            default:
            {
                MORIS_ERROR( false, " wrong master side ordinal ");
                break;
            }
        }

        // define an integration master element in the phys space
        Matrix< DDRMat > tXHatIG_M = tXHatIP_M;

        // define the master integration element in param space
        Matrix< DDRMat > tXiHatIG_M = tParamCoords;

        // get the node ids associated to the master side ordinal
        Matrix< DDSMat > tMasterSideNodes = tSideNodes.get_row( tMasterSideOrd );

        // phys coords amd parm coords in IP param space for the side
        uint tNumSideNodes = tMasterSideNodes.numel();
        Matrix< DDRMat > tMasterSidePhysCoords( tNumSideNodes, 2 );
        Matrix< DDRMat > tMasterSideParamCoords( tNumSideNodes, 3 );
        for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
        {
            tMasterSidePhysCoords.get_row( iNode )  = tXHatIG_M.get_row( tMasterSideNodes( iNode ) );
            tMasterSideParamCoords.get_row( iNode ) = tXiHatIG_M.get_row( tMasterSideNodes( iNode ) );
        }

        // loop over slave side ordinals
        for ( uint iSlaveOrd = iMasterOrd; iSlaveOrd < 3; iSlaveOrd++ )
        for ( uint iSlaveOrd = 0; iSlaveOrd < 3; iSlaveOrd++ )
        {
            // side ordinal for master integration mesh
            moris_index tSlaveSideOrd = iSlaveOrd;

            // define an interpolation slave element in the phys space
            Matrix< DDRMat > tXHatIP;
            switch ( iSlaveOrd )
            {
                case( 0 ):
                    tXHatIP = {{ 0.0, 1.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }};
                    break;
                case( 1 ):
                    tXHatIP = {{ 1.0, 1.0 }, { 0.0, 1.0 }, { 1.0, 0.0 }};
                    break;
                case( 2 ):
                    tXHatIP = {{ 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 }};
                    break;
                default:
                    MORIS_ERROR( false, " wrong slave side ordinal ");
                    break;
            }

            // for a master slave pair
            moris_index tSlaveNode = 0;
            //std::cout<<"A master slave pair + slave node"<<" "<<tMasterSideOrd<<" - "<<tSlaveSideOrd<<" - "<<tSlaveNode<<"-----------------------"<<std::endl;

            // get the rotation matrix
            Matrix< DDRMat > tR;
            rotation_matrix( tSideGeoTypeIG, tSlaveNode, tR );

            // define an interpolation slave element in the phys space
            Matrix< DDRMat > tXHatIP_S = tXHatIP;

            // define an integration slave element in the phys space
            Matrix< DDRMat > tXHatIG_S = tXHatIP_S;

            // define the slave integration element in param space
            Matrix< DDRMat > tXiHatIG_S = tParamCoords;

            // get the node ids associated to the slave side ordinal
            Matrix< DDSMat > tSlaveSideNodes = tSideNodes.get_row( tSlaveSideOrd );

            // phys coords amd parm coords in IP param space for the side
            Matrix< DDRMat > tSlaveSidePhysCoords( tNumSideNodes, 2 );
            Matrix< DDRMat > tSlaveSideParamCoords( tNumSideNodes, 3 );
            for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
            {
                tSlaveSidePhysCoords.get_row( iNode )  = tXHatIG_S.get_row( tSlaveSideNodes( iNode ) );
                tSlaveSideParamCoords.get_row( iNode ) = tXiHatIG_S.get_row( tSlaveSideNodes( iNode ) );
            }

            // bool for integration point check
            bool tIntegPointCheck = true;

            // loop over the integration points
            for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
            {
                // get the treated integration point location in integration space
                Matrix< DDRMat > tMasterIntegPointI = tIntegPoints.get_column( iGP );

                // get the treated integration point location in integration space
                Matrix< DDRMat > tSlaveIntegPointI = tMasterIntegPointI;
                tSlaveIntegPointI({0,0}, {0,0}) = tR * tMasterIntegPointI({0,0}, {0,0});

                // get the integration point in the IP parametric space
                // via the side shape functions for the master side
                Matrix< DDRMat > tMasterRefIntegPointI( 4, 1, 0.0 );
                Matrix< DDRMat > tN1;
                tSideSpaceInterp->eval_N( tMasterIntegPointI({0,0},{0,0}),tN1 );
                tMasterRefIntegPointI( { 0, 2 }, { 0, 0 } ) = trans( tN1 * tMasterSideParamCoords );
                tMasterRefIntegPointI( 3 ) = tMasterIntegPointI( 1 );

                // get the integration point in the IP parametric space
                // via the rotation matrix for the slave side
                Matrix< DDRMat > tSlaveRefIntegPointI( 4, 1, 0.0 );
                Matrix< DDRMat > tN2;
                tSideSpaceInterp->eval_N( tSlaveIntegPointI({0,0},{0,0}),tN2);
                tSlaveRefIntegPointI( { 0, 2 }, { 0, 0 } ) = trans( tN2 * tSlaveSideParamCoords );
                tSlaveRefIntegPointI( 3 ) = tSlaveIntegPointI( 1 );

                // to check only
                Matrix< DDRMat > tN3;
                tSpaceInterp->eval_N( tMasterRefIntegPointI,tN3 );
                Matrix< DDRMat > tMasterPhysIntegPointI = trans( tN3 * tXHatIP_M );

                // to check only
                Matrix< DDRMat > tN4;
                tSpaceInterp->eval_N( tSlaveRefIntegPointI,tN4 );
                Matrix< DDRMat > tSlavePhysIntegPointI = trans( tN4 * tXHatIP_S );

                // check the integration point in the IP parametric space
                for( uint iCoords = 0; iCoords < 2; iCoords++ )
                {
                    tIntegPointCheck = tIntegPointCheck && ( std::abs( tMasterPhysIntegPointI( iCoords ) - tSlavePhysIntegPointI( iCoords ) ) < tEpsilon );
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
    Matrix< DDSMat > tSideNodes = { { 0, 1, 3 }, { 1, 2, 3 }, { 0, 3, 2 }, { 0, 2, 1 } };
    Matrix< DDSMat > tFacingSideNodes = { { 2 }, { 0 }, { 1 }, { 3 } };

    // possible permutation cases
    Matrix< IndexMat > tPermuteCases = { { 0, 1, 2 }, { 1, 2, 0 }, { 2, 0, 1 } };

    // hex param coords
    Matrix< DDRMat > tParamCoords = {{ 1.0, 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 },
                                     { 0.0, 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0, 1.0 }};

    // interpolation rule for the volume
    mtk::Interpolation_Rule tInterpIGRule( tGeoTypeIG,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR );

    // create a space interpolation function
    mtk::Interpolation_Function_Base* tSpaceInterp = tInterpIGRule.create_space_interpolation_function();

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule( tSideGeoTypeIG,
                                          mtk::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR,
                                          mtk::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp = tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule for the side
    mtk::Integration_Rule tSideIntegRule( tSideGeoTypeIG,
                                     mtk::Integration_Type::GAUSS,
                                     mtk::Integration_Order::TRI_3,
                                     mtk::Integration_Type::GAUSS,
                                     mtk::Integration_Order::BAR_1 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );
    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // loop over master side ordinals
    for ( uint iMasterOrd = 0; iMasterOrd < 4; iMasterOrd++ )
    {
        // side ordinal for master integration mesh
        moris_index tMasterSideOrd = iMasterOrd;

        // define an interpolation master element in the phys space
        Matrix< DDRMat > tXHatIP_M;
        switch( tMasterSideOrd )
        {
            case( 0 ):
            {
                tXHatIP_M = {{ 1.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 }}; break;
            }
            case( 1 ):
            {
                tXHatIP_M = {{ 0.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 }}; break;
            }
            case( 2 ):
            {
                tXHatIP_M = {{ 1.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 0.0, 0.0 }}; break;
            }
            case( 3 ):
            {
                tXHatIP_M = {{ 1.0, 1.0, 0.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }}; break;
            }
            default:
            {
                MORIS_ERROR( false, " wrong master side ordinal "); break;
            }
        }

        // define an integration master element in the phys space
        Matrix< DDRMat > tXHatIG_M = tXHatIP_M;

        // define the master integration element in param space
        Matrix< DDRMat > tXiHatIG_M = tParamCoords;

        // get the node ids associated to the master side ordinal
        Matrix< DDSMat > tMasterSideNodes = tSideNodes.get_row( tMasterSideOrd );

        // phys coords amd parm coords in IP param space for the side
        uint tNumSideNodes = tMasterSideNodes.numel();
        Matrix< DDRMat > tMasterSidePhysCoords( tNumSideNodes, 3 );
        Matrix< DDRMat > tMasterSideParamCoords( tNumSideNodes, 4 );
        for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
        {
            tMasterSidePhysCoords.get_row( iNode )  = tXHatIG_M.get_row( tMasterSideNodes( iNode ) );
            tMasterSideParamCoords.get_row( iNode ) = tXiHatIG_M.get_row( tMasterSideNodes( iNode ) );
        }

        // loop over slave side ordinals
        for ( uint iSlaveOrd = 0; iSlaveOrd < 4; iSlaveOrd++ )
        {
            // side ordinal for master integration mesh
            moris_index tSlaveSideOrd = iSlaveOrd;

            // sort the nodes on side or not
            Matrix< DDSMat > tNodesOnSlaveSide    = tSideNodes.get_row( iSlaveOrd );
            Matrix< DDSMat > tNodesNotOnSlaveSide = tFacingSideNodes.get_row( iSlaveOrd );

            // define an interpolation slave element in the phys space
            Matrix< DDRMat > tXHatIP;
            switch( tSlaveSideOrd )
            {
                case( 0 ):
                {
                    tXHatIP = {{ 1.0, 1.0, 0.0 }, { 1.0, 0.0, 1.0 }, { 2.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }}; break;
                }
                case( 1 ):
                {
                    tXHatIP = {{ 2.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 0.0, 0.0 }}; break;
                }
                case( 2 ):
                {
                    tXHatIP = {{ 1.0, 1.0, 0.0 }, { 2.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 }}; break;
                }
                case( 3 ):
                {
                    tXHatIP = {{ 1.0, 1.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 }, { 2.0, 0.0, 0.0 }}; break;
                }
                default:
                {
                    MORIS_ERROR( false, " wrong slave side ordinal ");
                    break;
                }
            }

            // node permutation for $ rotation in plane
            for( uint  iPermute = 0; iPermute < tPermuteCases.n_rows(); iPermute++ )
            {
                // get the first corresponding node on the slave side
                moris_index tSlaveNode = iPermute;

                // for a master slave pair
                //std::cout<<"A master slave pair + slave node"<<" "<<tMasterSideOrd<<" - "<<tSlaveSideOrd<<" - "<<tSlaveNode<<"-----------------------"<<std::endl;

                // get the rotation matrix
                Matrix< DDRMat > tR;
                rotation_matrix( tSideGeoTypeIG, tSlaveNode, tR );

                // define an interpolation slave element in the phys space
                Matrix< DDRMat > tXHatIP_S( 4, 3, 0.0 );

                // fill the phys coords
                for(uint  i = 0; i < 3 ; i++)
                {
                    moris_index tNodeSide    = tNodesOnSlaveSide( tPermuteCases( iPermute, i ) );
                    moris_index tNodeNotSide = tNodesNotOnSlaveSide( 0 );
                    // place the interface nodes
                    tXHatIP_S.get_row( tNodeSide )    = tXHatIP.get_row( tNodesOnSlaveSide( i ) );
                    tXHatIP_S.get_row( tNodeNotSide ) = tXHatIP.get_row( tNodesNotOnSlaveSide( 0 ) );
                }

                // define an integration slave element in the phys space
                Matrix< DDRMat > tXHatIG_S = tXHatIP_S;

                // define the slave integration element in param space
                Matrix< DDRMat > tXiHatIG_S = tParamCoords;

                // get the node ids associated to the master side ordinal
                Matrix< DDSMat > tSlaveSideNodes = tSideNodes.get_row( tSlaveSideOrd );

                // phys coords amd parm coords in IP param space for the side
                Matrix< DDRMat > tSlaveSidePhysCoords( tNumSideNodes, 3 );
                Matrix< DDRMat > tSlaveSideParamCoords( tNumSideNodes, 4 );
                for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
                {
                    tSlaveSidePhysCoords.get_row( iNode )  = tXHatIG_S.get_row( tSlaveSideNodes( iNode ) );
                    tSlaveSideParamCoords.get_row( iNode ) = tXiHatIG_S.get_row( tSlaveSideNodes( iNode ) );
                }

                // bool for integration point check
                bool tIntegPointCheck = true;

                // loop over the integration points
                for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tMasterIntegPointI = tIntegPoints.get_column( iGP );

                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tSlaveIntegPointI = tMasterIntegPointI;
                    tSlaveIntegPointI({0,2}, {0,0}) = tR * tMasterIntegPointI({0,2}, {0,0});

                    // get the integration point in the IP parametric space
                    // via the side shape functions for the master side
                    Matrix< DDRMat > tN1;
                    tSideSpaceInterp->eval_N( tMasterIntegPointI({0,2}, {0,0}),tN1 );
                    Matrix< DDRMat > tMasterRefIntegPointI = trans( tN1 * tMasterSideParamCoords );

                    // get the integration point in the IP parametric space
                    // via the rotation matrix for the slave side
                    Matrix< DDRMat > tN2;
                    tSideSpaceInterp->eval_N( tSlaveIntegPointI({0,2}, {0,0}), tN2 );
                    Matrix< DDRMat > tSlaveRefIntegPointI = trans( tN2 * tSlaveSideParamCoords );

                    // to check only
                    Matrix< DDRMat > tN3;
                    tSpaceInterp->eval_N( tMasterRefIntegPointI, tN3 );
                    Matrix< DDRMat > tMasterPhysIntegPointI = trans( tN3 * tXHatIP_M );

                    // to check only
                    Matrix< DDRMat > tN4;
                    tSpaceInterp->eval_N( tSlaveRefIntegPointI, tN4 );
                    Matrix< DDRMat > tSlavePhysIntegPointI = trans( tN4 * tXHatIP_S );

                    // check the integration point in the IP parametric space
                    for( uint iCoords = 0; iCoords < 3; iCoords++ )
                    {
                        tIntegPointCheck = tIntegPointCheck && ( std::abs( tMasterPhysIntegPointI( iCoords ) - tSlavePhysIntegPointI( iCoords ) ) < tEpsilon );
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
    Matrix< DDSMat > tSideNodes = { { 0, 1, 5, 4 }, { 1, 2, 6, 5 }, { 2, 3, 7, 6 },
                                    { 0, 4, 7, 3 }, { 0, 3, 2, 1 }, { 4, 5, 6, 7 } };
    Matrix< DDSMat > tFacingSideNodes = { { 3, 2, 6, 7 }, { 0, 3, 7, 4 }, { 1, 0, 4, 5 },
                                          { 1, 5, 6, 2 }, { 4, 7, 6, 5 }, { 0, 1, 2, 3 } };

    // possible permutation cases
    Matrix< IndexMat > tPermuteCases = { { 0, 1, 2, 3 }, { 1, 2, 3, 0 }, { 2, 3, 0, 1 }, { 3, 0, 1, 2 } };

    // hex param coords
    Matrix< DDRMat > tHexParamCoords = {{ -1.0, -1.0, -1.0 }, {  1.0, -1.0, -1.0 },
                                        {  1.0,  1.0, -1.0 }, { -1.0,  1.0, -1.0 },
                                        { -1.0, -1.0,  1.0 }, {  1.0, -1.0,  1.0 },
                                        {  1.0,  1.0,  1.0 }, { -1.0,  1.0,  1.0 }};

    // interpolation rule for the volume
    mtk::Interpolation_Rule tInterpIGRule( tGeoTypeIG,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR );

    // create a space interpolation function
    mtk::Interpolation_Function_Base* tSpaceInterp = tInterpIGRule.create_space_interpolation_function();

    // interpolation rule for the side
    mtk::Interpolation_Rule tSideInterpIGRule( tSideGeoTypeIG,
                                          mtk::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR,
                                          mtk::Interpolation_Type::LAGRANGE,
                                          mtk::Interpolation_Order::LINEAR );

    // create a side space interpolation function
    mtk::Interpolation_Function_Base* tSideSpaceInterp = tSideInterpIGRule.create_space_interpolation_function();

    // create a integration rule for the side
    mtk::Integration_Rule tSideIntegRule( tSideGeoTypeIG,
                                     mtk::Integration_Type::GAUSS,
                                     mtk::Integration_Order::QUAD_5x5,
                                     mtk::Integration_Type::GAUSS,
                                     mtk::Integration_Order::BAR_1 );

    // create a side integrator
    mtk::Integrator tSideIntegrator( tSideIntegRule );

    // get number of integration points, integration points and weights
    uint             tNumOfIntegPoints = tSideIntegrator.get_number_of_points();
    Matrix< DDRMat > tIntegPoints;
    tSideIntegrator.get_points( tIntegPoints );
    Matrix< DDRMat > tIntegWeights;
    tSideIntegrator.get_weights( tIntegWeights );

    // loop over master side ordinals
    for ( uint iMasterOrd = 0; iMasterOrd < 6; iMasterOrd++ )
    {
        // side ordinal for master integration mesh
        moris_index tMasterSideOrd = iMasterOrd;

        // define an interpolation master element in the phys space
        Matrix< DDRMat > tXHatIP_M;
        switch ( iMasterOrd )
        {
            case( 0 ):
                tXHatIP_M = {{ 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 },
                             { 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 },
                             { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 },
                             { 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 }};
                break;
            case( 1 ):
                tXHatIP_M = {{ 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 },
                             { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
                             { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 },
                             { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 }};
                break;
            case( 2 ):
                tXHatIP_M = {{ 0.0, 1.0, 0.0 }, { 0.0, 0.0, 0.0 },
                             { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 },
                             { 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 },
                             { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }};
                break;
            case( 3 ):
                tXHatIP_M = {{ 1.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 },
                             { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 },
                             { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
                             { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }};
                break;
            case( 4 ):
                tXHatIP_M = {{ 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 },
                             { 1.0, 1.0, 1.0 }, { 1.0, 1.0, 0.0 },
                             { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 1.0 },
                             { 0.0, 1.0, 1.0 }, { 0.0, 1.0, 0.0 }};
                break;
            case( 5 ):
                tXHatIP_M = {{ 0.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 },
                             { 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 },
                             { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 },
                             { 1.0, 1.0, 1.0 }, { 1.0, 0.0, 1.0 }};
                break;
            default:
            {
                MORIS_ERROR( false, " wrong master side ordinal ");
                break;
            }
        }

        // define an integration master element in the phys space
        Matrix< DDRMat > tXHatIG_M = tXHatIP_M;

        // define the master integration element in param space
        Matrix< DDRMat > tXiHatIG_M = tHexParamCoords;

        // get the node ids associated to the master side ordinal
        Matrix< DDSMat > tMasterSideNodes = tSideNodes.get_row( tMasterSideOrd );

        // phys coords amd parm coords in IP param space for the side
        uint tNumSideNodes = tMasterSideNodes.numel();
        Matrix< DDRMat > tMasterSidePhysCoords( tNumSideNodes, 3 );
        Matrix< DDRMat > tMasterSideParamCoords( tNumSideNodes, 3 );
        for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
        {
            tMasterSidePhysCoords.get_row( iNode )  = tXHatIG_M.get_row( tMasterSideNodes( iNode ) );
            tMasterSideParamCoords.get_row( iNode ) = tXiHatIG_M.get_row( tMasterSideNodes( iNode ) );
        }

        // loop over slave side ordinals
        //for ( uint iSlaveOrd = iMasterOrd; iSlaveOrd < 6; iSlaveOrd++ )
        for ( uint iSlaveOrd = 0; iSlaveOrd < 6; iSlaveOrd++ )
        {
            // side ordinal for master integration mesh
            moris_index tSlaveSideOrd = iSlaveOrd;

            // sort the nodes on side or not
            Matrix< DDSMat > tNodesOnSlaveSide   = tSideNodes.get_row( iSlaveOrd );
            Matrix< DDSMat > tNodeNotOnSlaveSide = tFacingSideNodes.get_row( iSlaveOrd );

            // define an interpolation slave element in the phys space
            Matrix< DDRMat > tXHatIP;
            switch ( iSlaveOrd )
            {
                case( 0 ):
                    tXHatIP = {{ 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 },
                               { 2.0, 0.0, 1.0 }, { 2.0, 0.0, 0.0 },
                               { 1.0, 1.0, 0.0 }, { 1.0, 1.0, 1.0 },
                               { 2.0, 1.0, 1.0 }, { 2.0, 1.0, 0.0 }};
                    break;
                case( 1 ):
                    tXHatIP = {{ 2.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 },
                               { 1.0, 0.0, 1.0 }, { 2.0, 0.0, 1.0 },
                               { 2.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 },
                               { 1.0, 1.0, 1.0 }, { 2.0, 1.0, 1.0 }};
                    break;
                case( 2 ):
                    tXHatIP = {{ 2.0, 0.0, 1.0 }, { 2.0, 0.0, 0.0 },
                               { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 },
                               { 2.0, 1.0, 1.0 }, { 2.0, 1.0, 0.0 },
                               { 1.0, 1.0, 0.0 }, { 1.0, 1.0, 1.0 }};
                    break;
                case( 3 ):
                    tXHatIP = {{ 1.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0 },
                               { 2.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 },
                               { 1.0, 0.0, 1.0 }, { 2.0, 0.0, 1.0 },
                               { 2.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }};
                    break;
                case( 4 ):
                    tXHatIP = {{ 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 },
                               { 1.0, 1.0, 1.0 }, { 1.0, 0.0, 1.0 },
                               { 2.0, 0.0, 0.0 }, { 2.0, 1.0, 0.0 },
                               { 2.0, 1.0, 1.0 }, { 2.0, 0.0, 1.0 }};
                    break;
                case( 5 ):
                    tXHatIP = {{ 2.0, 0.0, 0.0 }, { 2.0, 0.0, 1.0 },
                               { 2.0, 1.0, 1.0 }, { 2.0, 1.0, 0.0 },
                               { 1.0, 0.0, 0.0 }, { 1.0, 0.0, 1.0 },
                               { 1.0, 1.0, 1.0 }, { 1.0, 1.0, 0.0 }};
                    break;
                default:
                    MORIS_ERROR( false, " wrong slave side ordinal ");
                    break;
            }


            // node permutation for $ rotation in plane
            for( uint  iPermute = 0; iPermute < tPermuteCases.n_rows(); iPermute++ )
            {
                // get the first corresponding node on the slave side
                moris_index tSlaveNode = iPermute;

                // for a master slave pair
                //std::cout<<"A master slave pair + slave node"<<" "<<tMasterSideOrd<<" - "<<tSlaveSideOrd<<" - "<<tSlaveNode<<"-----------------------"<<std::endl;

                // get the rotation matrix
                Matrix< DDRMat > tR;
                rotation_matrix( tSideGeoTypeIG, tSlaveNode, tR );

                // define an interpolation slave element in the phys space
                Matrix< DDRMat > tXHatIP_S( 8, 3, 0.0 );

                // fill the phys coords
                for(uint  i = 0; i < 4 ; i++)
                {
                    moris_index tNodeSide    = tNodesOnSlaveSide( tPermuteCases( iPermute, i ) );
                    moris_index tNodeNotSide = tNodeNotOnSlaveSide( tPermuteCases( iPermute, i ) );

                    // place the interface nodes
                    tXHatIP_S.get_row( tNodeSide )    = tXHatIP.get_row( tNodesOnSlaveSide( i ) );
                    tXHatIP_S.get_row( tNodeNotSide ) = tXHatIP.get_row( tNodeNotOnSlaveSide( i ) );
                }
                //print(tXHatIP_S,"tXHatIP_S");

                // define an integration slave element in the phys space
                Matrix< DDRMat > tXHatIG_S = tXHatIP_S;

                // define the slave integration element in param space
                Matrix< DDRMat > tXiHatIG_S = tHexParamCoords;

                // get the node ids associated to the master side ordinal
                Matrix< DDSMat > tSlaveSideNodes = tSideNodes.get_row( tSlaveSideOrd );

                // phys coords and parm coords in IP param space for the side
                Matrix< DDRMat > tSlaveSidePhysCoords( tNumSideNodes, 3 );
                Matrix< DDRMat > tSlaveSideParamCoords( tNumSideNodes, 3 );
                for( uint iNode = 0; iNode < tNumSideNodes; iNode++ )
                {
                    tSlaveSidePhysCoords.get_row( iNode )  = tXHatIG_S.get_row( tSlaveSideNodes( iNode ) );
                    tSlaveSideParamCoords.get_row( iNode ) = tXiHatIG_S.get_row( tSlaveSideNodes( iNode ) );
                }

                // bool for integration point check
                bool tIntegPointCheck = true;

                // loop over the integration points
                for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
                {
                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tMasterIntegPointI = tIntegPoints.get_column( iGP );

                    // get the treated integration point location in integration space
                    Matrix< DDRMat > tSlaveIntegPointI = tMasterIntegPointI;
                    tSlaveIntegPointI({0,1}, {0,0}) = tR * tMasterIntegPointI({0,1}, {0,0});

                    // get the integration point in the IP parametric space
                    // via the side shape functions for the master side
                    Matrix< DDRMat > tN1;
                    tSideSpaceInterp->eval_N( tMasterIntegPointI({0,1}, {0,0}),tN1 );
                    Matrix< DDRMat > tMasterRefIntegPointI = trans( tN1 * tMasterSideParamCoords );

                    // get the integration point in the IP parametric space
                    // via the rotation matrix for the slave side
                    Matrix< DDRMat > tN2;
                    tSideSpaceInterp->eval_N( tSlaveIntegPointI({0,1}, {0,0}), tN2 );
                    Matrix< DDRMat > tSlaveRefIntegPointI = trans(tN2  * tSlaveSideParamCoords );

                    // to check only
                    Matrix< DDRMat > tN3;
                    tSpaceInterp->eval_N( tMasterRefIntegPointI,tN3 );
                    Matrix< DDRMat > tMasterPhysIntegPointI = trans( tN3 * tXHatIP_M );

                    // to check only
                    Matrix< DDRMat > tN4;
                    tSpaceInterp->eval_N( tSlaveRefIntegPointI,tN4 );
                    Matrix< DDRMat > tSlavePhysIntegPointI = trans( tN4 * tXHatIP_S );

                    // check the integration point in the IP parametric space
                    for( uint iCoords = 0; iCoords < 3; iCoords++ )
                    {
                        tIntegPointCheck = tIntegPointCheck && ( std::abs( tMasterPhysIntegPointI( iCoords ) - tSlavePhysIntegPointI( iCoords ) ) < tEpsilon );
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
