/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Struc_Linear_Contact_Nitsche.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private public
// FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
// MTK/src
#include "cl_MTK_Enums.hpp"
// FEM//INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
// LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

void Test_IWG_Struc_Contact_Remap(
        const IWG_Type aIWGType,
        bool           aUseDeformedConfig,
        bool           aUseConsistentDeformedConfig,
        const real&    aDisplacementScaling )
{
    // define an epsilon environment
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType     = mtk::Geometry_Type::UNDEFINED;
    mtk::Geometry_Type tSideGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    Vector< mtk::Interpolation_Order > tInterpolationOrders = {
        mtk::Interpolation_Order::LINEAR,
        mtk::Interpolation_Order::QUADRATIC,
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };

    // dof type list
    Vector< Vector< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( { { { 10.0 } } } );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( { { { 0.3 } } } );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropFollowerEMod = std::make_shared< fem::Property >();
    tPropFollowerEMod->set_parameters( { { { 20.0 } } } );
    tPropFollowerEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropFollowerNu = std::make_shared< fem::Property >();
    tPropFollowerNu->set_parameters( { { { 0.3 } } } );
    tPropFollowerNu->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMFollowerStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMFollowerStrucLinIso->set_property( tPropFollowerEMod, "YoungsModulus" );
    tCMFollowerStrucLinIso->set_property( tPropFollowerNu, "PoissonRatio" );
    tCMFollowerStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { { { 100.0 } } } );
    tSPNitscheInterface->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );
    tSPNitscheInterface->set_property( tPropFollowerEMod, "Material", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitscheInterface->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( aIWGType );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerStrucLinIso, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set gap computing parameters
    tIWG->mUseDeformedGeometryForGap           = aUseDeformedConfig;
    tIWG->mUseConsistentDeformedGeometryForGap = aUseConsistentDeformedConfig;

    // set error flag
    bool tPassedCheck = true;

    // loop on the space dimension - Note: this test is only for 2D
    for ( uint iSpaceDim = 2; iSpaceDim < 3; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set space dim for CM
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );
        tCMFollowerStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMFollowerStrucLinIso->set_space_dim( iSpaceDim );
        tSPNitscheInterface->set_space_dim( iSpaceDim );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type of background and side
                tGeometryType     = mtk::Geometry_Type::QUAD;
                tSideGeometryType = mtk::Geometry_Type::LINE;
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType     = mtk::Geometry_Type::HEX;
                tSideGeometryType = mtk::Geometry_Type::TRI;
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // loop on the field interpolation order
        // only up to quadratic checked for now
        for ( uint iInterpOrder = 1; iInterpOrder < 3; iInterpOrder++ )
        {
            // create an interpolation order
            mtk::Interpolation_Order tGIInterpolationOrder   = tInterpolationOrders( iInterpOrder - 1 );
            mtk::Interpolation_Order tSideInterpolationOrder = aUseConsistentDeformedConfig ? tGIInterpolationOrder : tInterpolationOrders( 0 );

            // space and time geometry interpolators
            //------------------------------------------------------------------------------
            // create a space geometry interpolation rule (same for leader and follower)
            mtk::Interpolation_Rule tGIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tGIInterpolationOrder,
                    mtk::Interpolation_Type::CONSTANT,
                    mtk::Interpolation_Order::CONSTANT );

            mtk::Interpolation_Rule tSideGIRule(
                    tSideGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tSideInterpolationOrder,
                    mtk::Interpolation_Type::CONSTANT,
                    mtk::Interpolation_Order::CONSTANT );

            // create a space time geometry interpolator
            Geometry_Interpolator tLeaderGI   = Geometry_Interpolator( tGIRule );
            Geometry_Interpolator tFollowerGI = Geometry_Interpolator( tGIRule );

            Geometry_Interpolator tLeaderSideGI   = Geometry_Interpolator( tSideGIRule, mtk::CellShape::GENERAL, true );
            Geometry_Interpolator tFollowerSideGI = Geometry_Interpolator( tSideGIRule, mtk::CellShape::GENERAL, true );

            // create time coeff tHat
            Matrix< DDRMat > tTHat = { { 0.0 } };
            Matrix< DDRMat > tLeaderXHat;
            Matrix< DDRMat > tFollowerXHat;

            Matrix< DDRMat > tLeaderSideNodesParam;
            Matrix< DDRMat > tFollowerSideNodesParam;

            Matrix< DDRMat > tLeaderDOFHatDisp;
            Matrix< DDRMat > tFollowerDOFHatDisp;

            Matrix< DDRMat > tIntegPoints;

            switch ( iSpaceDim )
            {
                case 2:
                {
                    tIntegPoints = { { 0.123 }, { 0.0 } };

                    switch ( iInterpOrder )
                    {
                        case 1:
                        {
                            tLeaderXHat = { { 0.0, 0.0 },
                                { 1.0, 0.0 },
                                { 1.0, 1.0 },
                                { 0.0, 1.0 } };

                            tFollowerXHat = tLeaderXHat;

                            tLeaderSideNodesParam   = { { -0.9, -0.3 }, { 0.8, -0.1 } };
                            tFollowerSideNodesParam = { { -0.8, 0.8 }, { 0.9, 0.5 } };

                            tLeaderDOFHatDisp   = { { -0.2, 0.1 }, { -0.1, 0.1 }, { -0.1, -0.9 }, { -0.2, 0.1 } };
                            tFollowerDOFHatDisp = { { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.5 }, { 0.0, 0.0 } };
                            break;
                        }
                        case 2:
                        {
                            tLeaderXHat = {
                                { 0.0, 0.0 },
                                { 1.0, 0.0 },
                                { 1.0, 1.0 },
                                { 0.0, 1.0 },
                                { 0.5, 0.0 },
                                { 1.0, 0.5 },
                                { 0.5, 1.0 },
                                { 0.0, 0.5 },
                                { 0.5, 0.5 }
                            };

                            tFollowerXHat = tLeaderXHat;

                            if ( aUseConsistentDeformedConfig )
                            {
                                tLeaderSideNodesParam   = { { -0.9, -0.3 }, { 0.8, -0.1 }, { 0.3, -0.1 } };
                                tFollowerSideNodesParam = { { -0.8, 0.8 }, { 0.9, 0.5 }, { 0.3, 0.7 } };
                            }
                            else
                            {
                                tLeaderSideNodesParam   = { { -0.9, -0.3 }, { 0.8, -0.1 } };
                                tFollowerSideNodesParam = { { -0.8, 0.8 }, { 0.9, 0.5 } };
                            }

                            tLeaderDOFHatDisp = {
                                { -0.2, 0.1 },
                                { -0.1, 0.1 },
                                { -0.1, -0.9 },
                                { -0.2, 0.1 },
                                { -0.15, 0.1 },
                                { -0.1, -0.4 },
                                { -0.15, -0.4 },
                                { -0.175, 0.1 },
                                { -0.15, -0.1 }
                            };
                            tFollowerDOFHatDisp = {
                                { 0.0, 0.0 },
                                { 0.0, 0.0 },
                                { 0.0, 0.5 },
                                { 0.0, 0.0 },
                                { 0.0, 0.0 },
                                { 0.0, 0.25 },
                                { 0.0, 0.25 },
                                { 0.0, 0.0 },
                                { 0.0, 0.1 }
                            };
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "Interpolation order not implemented." );
                        }
                    }
                    break;
                }
                default:
                    MORIS_ERROR( false, "Space dimension not implemented." );
            }

            // set the coefficients xHat, tHat
            tLeaderGI.set_coeff( tLeaderXHat, tTHat );
            tFollowerGI.set_coeff( tFollowerXHat, tTHat );

            // set the coefficients for the leader side geometry interpolators from parametric coordinates
            Matrix< DDRMat > tLeaderSideXHat( tLeaderSideNodesParam.n_rows(), iSpaceDim );

            for ( uint iNode = 0; iNode < tLeaderSideNodesParam.n_rows(); iNode++ )
            {
                tLeaderGI.set_space( trans( tLeaderSideNodesParam.get_row( iNode ) ) );
                tLeaderSideXHat.get_row( iNode ) = tLeaderGI.valx().matrix_data();
            }

            tLeaderSideGI.set_coeff( tLeaderSideXHat, tTHat );
            tLeaderSideGI.set_space_param_coeff( tLeaderSideNodesParam );
            tLeaderSideGI.set_time_param_coeff( tTHat );

            // set the coefficients for the follower side geometry interpolators from parametric coordinates
            Matrix< DDRMat > tFollowerSideXHat( tFollowerSideNodesParam.n_rows(), iSpaceDim );

            for ( uint iNode = 0; iNode < tFollowerSideNodesParam.n_rows(); iNode++ )
            {
                tFollowerGI.set_space( trans( tFollowerSideNodesParam.get_row( iNode ) ) );
                tFollowerSideXHat.get_row( iNode ) = tFollowerGI.valx().matrix_data();
            }

            tFollowerSideGI.set_coeff( tFollowerSideXHat, tTHat );
            tFollowerSideGI.set_space_param_coeff( tFollowerSideNodesParam );
            tFollowerSideGI.set_time_param_coeff( tTHat );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // number of dof for interpolation order
            uint tNumCoeff = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // get number of dof per type
            int tNumDofDisp = tNumCoeff * iSpaceDim;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::CONSTANT,
                    mtk::Interpolation_Order::CONSTANT );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tLeaderGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderDOFHatDisp );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tFollowerGI, tDispDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerDOFHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDispDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, 2 * tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofDisp - 1 },
                { tNumDofDisp, 2 * tNumDofDisp - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofDisp, 2 * tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes   = tDispDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes = tDispDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >   tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;

            Field_Interpolator_Manager tLeaderFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tFollowerFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderFIManager.mFI                       = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator   = &tLeaderGI;
            tLeaderFIManager.mIGGeometryInterpolator   = &tLeaderSideGI;
            tFollowerFIManager.mFI                     = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator = &tFollowerGI;
            tFollowerFIManager.mIGGeometryInterpolator = &tFollowerSideGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager   = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager = &tFollowerFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tLeaderFIManager, mtk::Leader_Follower::LEADER );
            tIWG->set_field_interpolator_manager( &tFollowerFIManager, mtk::Leader_Follower::FOLLOWER );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPointLeader = tIntegPoints.get_column( iGP );

                // set integration point on follower side
                tIWG->mSet->mLeaderFIManager->set_space_time_from_local_IG_point( tParamPointLeader );

                Matrix< DDRMat > tParamPointfollower = tIWG->remap_nonconformal_rays(
                        tIWG->mUseDeformedGeometryForGap,
                        tIWG->mUseConsistentDeformedGeometryForGap,
                        tDispDofTypes( 0 ),
                        tIWG->mSet->mLeaderFIManager,
                        tIWG->mSet->mFollowerFIManager,
                        tIWG->mGapData );

                // get reference to current gap data
                const std::unique_ptr< GapData >& tCurrentGapData = tIWG->get_gap_data();

                // save mapping results
                std::unique_ptr< GapData > tNominalGapData = std::make_unique< GapData >();
                tNominalGapData->copy( tCurrentGapData );

                // FD wrt displacements on leader side
                std::unique_ptr< GapData > tPertGapData = std::make_unique< GapData >();

                uint tCounter = 0;
                uint tNumDofs = tLeaderDOFHatDisp.numel();

                for ( uint idir = 0; idir < tLeaderDOFHatDisp.n_cols(); idir++ )
                {
                    for ( uint in = 0; in < tLeaderDOFHatDisp.n_rows(); in++ )
                    {
                        tLeaderDOFHatDisp( in, idir ) += tPerturbation;
                        tLeaderFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tIWG->mUseConsistentDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mGapData );
                        tPertGapData->copy( tCurrentGapData );

                        tLeaderDOFHatDisp( in, idir ) -= 2.0 * tPerturbation;
                        tLeaderFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tIWG->mUseConsistentDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mGapData );

                        real             tFdGapdu   = ( tPertGapData->mGap - tCurrentGapData->mGap ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdGap2du2 = ( tPertGapData->mdGapdu - tCurrentGapData->mdGapdu ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdGap2duv = ( tPertGapData->mdGapdv - tCurrentGapData->mdGapdv ) / 2 / tPerturbation;

                        real tErrordGapdu   = 100.0 * std::abs( tFdGapdu - tNominalGapData->mdGapdu( tCounter ) ) / std::abs( tFdGapdu );
                        real tErrordGap2du2 = 100.0 * norm( tFdGap2du2 - tNominalGapData->mdGap2du2.get_row( tCounter ) ) / norm( tFdGap2du2 );
                        real tErrordGap2duv = 100.0 * norm( tFdGap2duv - tNominalGapData->mdGap2duv.get_row( tCounter ) ) / norm( tFdGap2duv );

                        if ( tErrordGapdu + tErrordGap2du2 + tErrordGap2duv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error dGapdu:  %e\n", tErrordGapdu );
                            fprintf( stdout, "percent error dGap2du2 %e\n", tErrordGap2du2 );
                            fprintf( stdout, "percent error dGap2duv %e\n", tErrordGap2duv );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFdEtadu   = ( tPertGapData->mEta - tCurrentGapData->mEta ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdEta2du2 = ( tPertGapData->mdEtadu - tCurrentGapData->mdEtadu ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdEta2duv = ( tPertGapData->mdEtadv - tCurrentGapData->mdEtadv ) / 2 / tPerturbation;

                        real tErrordEtadu   = 100.0 * norm( tFdEtadu - tNominalGapData->mdEtadu( tCounter ) ) / norm( tFdEtadu );
                        real tErrordEta2du2 = 100.0 * norm( tFdEta2du2 - tNominalGapData->mdEta2du2.get_row( tCounter ) ) / norm( tFdEta2du2 );
                        real tErrordEta2duv = 100.0 * norm( tFdEta2duv - tNominalGapData->mdEta2duv.get_row( tCounter ) ) / norm( tFdEta2duv );

                        if ( tErrordEtadu + tErrordEta2du2 + tErrordEta2duv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error dEtadu:   %e\n", tErrordEtadu );
                            fprintf( stdout, "percent error dEta2du2: %e\n", tErrordEta2du2 );
                            fprintf( stdout, "percent error dEta2duv: %e\n", tErrordEta2duv );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFdLeaderNormaldu    = ( tPertGapData->mLeaderNormal - tCurrentGapData->mLeaderNormal ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdLeaderdNormal2du2 = ( tPertGapData->mLeaderdNormaldu - tCurrentGapData->mLeaderdNormaldu ) / 2 / tPerturbation;

                        std::pair< moris::size_t, moris::size_t > tRowBlk = { 0, iSpaceDim - 1 };
                        std::pair< moris::size_t, moris::size_t > tColBlk = { tCounter * tNumDofs, ( tCounter + 1 ) * tNumDofs - 1 };

                        real tErrordLeaderNormaldu    = 100.0 * norm( tFdLeaderNormaldu - tNominalGapData->mLeaderdNormaldu.get_column( tCounter ) ) / norm( tFdLeaderNormaldu );
                        real tErrordLeaderdNormal2du2 = 100.0 * norm( tFdLeaderdNormal2du2 - tNominalGapData->mLeaderdNormal2du2( tRowBlk, tColBlk ) ) / norm( tFdLeaderdNormal2du2 );

                        if ( tErrordLeaderNormaldu + tErrordLeaderdNormal2du2 > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error dLeaderNormaldu:   %e\n", tErrordLeaderNormaldu );
                            fprintf( stdout, "percent error dLeaderdNormal2du2: %e\n", tErrordLeaderdNormal2du2 );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFdLeaderRefNormaldu = ( tPertGapData->mLeaderRefNormal - tCurrentGapData->mLeaderRefNormal ) / 2 / tPerturbation;

                        if ( norm( tFdLeaderRefNormaldu ) > MORIS_REAL_EPS )
                        {
                            fprintf( stdout, "Error - Reference normal changes with leader displacement\n" );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFdGapVecdu   = ( tPertGapData->mGapVec - tCurrentGapData->mGapVec ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdGapvec2du2 = ( tPertGapData->mdGapvecdu - tCurrentGapData->mdGapvecdu ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdGapvec2duv = ( tPertGapData->mdGapvecdv - tCurrentGapData->mdGapvecdv ) / 2 / tPerturbation;

                        real tErrordGapVectordu = 100.0 * norm( tFdGapVecdu - tNominalGapData->mdGapvecdu.get_column( tCounter ) ) / norm( tFdGapVecdu );
                        real tErrordGapvec2du2  = 100.0 * norm( tFdGapvec2du2 - tNominalGapData->mdGapvec2du2( tRowBlk, tColBlk ) ) / norm( tFdGapvec2du2 );
                        real tErrordGapvec2duv  = 100.0 * norm( tFdGapvec2duv - tNominalGapData->mdGapvec2duv( tRowBlk, tColBlk ) ) / norm( tFdGapvec2duv );

                        if ( tErrordGapVectordu + tErrordGapvec2du2 + tErrordGapvec2duv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error dGapVecdu:   %e\n", tErrordGapVectordu );
                            fprintf( stdout, "percent error dGapvec2du2: %e\n", tErrordGapvec2du2 );
                            fprintf( stdout, "percent error dGapvec2duv: %e\n", tErrordGapvec2duv );
                            tPassedCheck = false;
                        }

                        tLeaderDOFHatDisp( in, idir ) += tPerturbation;
                        tCounter++;
                    }
                }
                // reset the leader dof hat displacement
                tLeaderFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderDOFHatDisp );

                // FD wrt displacements on leader side
                tCounter       = 0;
                uint tNumNodes = tFollowerDOFHatDisp.n_rows();

                for ( uint idir = 0; idir < tFollowerDOFHatDisp.n_cols(); idir++ )
                {
                    for ( uint in = 0; in < tFollowerDOFHatDisp.n_rows(); in++ )
                    {
                        tFollowerDOFHatDisp( in, idir ) += tPerturbation;
                        tFollowerFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tIWG->mUseConsistentDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mGapData );
                        tPertGapData->copy( tCurrentGapData );

                        tFollowerDOFHatDisp( in, idir ) -= 2.0 * tPerturbation;
                        tFollowerFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tIWG->mUseConsistentDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mGapData );

                        real             tFdGapdv   = ( tPertGapData->mGap - tCurrentGapData->mGap ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdGap2dv2 = ( tPertGapData->mdGapdv - tCurrentGapData->mdGapdv ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdGap2duv = ( tPertGapData->mdGapdu - tCurrentGapData->mdGapdu ) / 2 / tPerturbation;

                        real tErrordGapdv   = 100.0 * std::abs( tFdGapdv - tNominalGapData->mdGapdv( tCounter ) ) / std::abs( tFdGapdv );
                        real tErrordGap2dv2 = 100.0 * norm( tFdGap2dv2 - tNominalGapData->mdGap2dv2.get_row( tCounter ) ) / norm( tFdGap2dv2 );
                        real tErrordGap2duv = 100.0 * norm( tFdGap2duv - trans( tNominalGapData->mdGap2duv.get_column( tCounter ) ) ) / norm( tFdGap2duv );

                        if ( tErrordGapdv + tErrordGap2dv2 + tErrordGap2duv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error dGapdv:  %e\n", tErrordGapdv );
                            fprintf( stdout, "percent error dGap2dv2 %e\n", tErrordGap2dv2 );
                            fprintf( stdout, "percent error dGap2duv %e\n", tErrordGap2duv );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFdEtadv   = ( tPertGapData->mEta - tCurrentGapData->mEta ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdEta2dv2 = ( tPertGapData->mdEtadv - tCurrentGapData->mdEtadv ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdEta2duv = ( tPertGapData->mdEtadu - tCurrentGapData->mdEtadu ) / 2 / tPerturbation;

                        real tErrordEtadv   = 100.0 * norm( tFdEtadv - tNominalGapData->mdEtadv( tCounter ) ) / norm( tFdEtadv );
                        real tErrordEta2dv2 = 100.0 * norm( tFdEta2dv2 - tNominalGapData->mdEta2dv2.get_row( tCounter ) ) / norm( tFdEta2dv2 );
                        real tErrordEta2duv = 100.0 * norm( tFdEta2duv - trans( tNominalGapData->mdEta2duv.get_column( tCounter ) ) ) / norm( tFdEta2duv );

                        if ( tErrordEtadv + tErrordEta2dv2 + tErrordEta2duv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error dEtadv:   %e\n", tErrordEtadv );
                            fprintf( stdout, "percent error dEta2dv2: %e\n", tErrordEta2dv2 );
                            fprintf( stdout, "percent error dEta2duv: %e\n", tErrordEta2duv );
                            tPassedCheck = false;
                        }
                        Matrix< DDRMat > tFdLeaderNormaldu    = ( tPertGapData->mLeaderNormal - tCurrentGapData->mLeaderNormal ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdLeaderdNormal2du2 = ( tPertGapData->mLeaderdNormaldu - tCurrentGapData->mLeaderdNormaldu ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdLeaderRefNormaldu = ( tPertGapData->mLeaderRefNormal - tCurrentGapData->mLeaderRefNormal ) / 2 / tPerturbation;

                        if ( norm( tFdLeaderNormaldu ) + norm( tFdLeaderdNormal2du2 ) + norm( tFdLeaderRefNormaldu ) > MORIS_REAL_EPS )
                        {
                            fprintf( stdout, "Error - Deformed normal, its derivative, and/or reference normal changes with follower displacement\n" );
                            tPassedCheck = false;
                        }

                        std::pair< moris::size_t, moris::size_t > tRowBlk = { 0, iSpaceDim - 1 };
                        std::pair< moris::size_t, moris::size_t > tColBlk = { tCounter * tNumDofs, ( tCounter + 1 ) * tNumDofs - 1 };

                        Matrix< DDRMat > tFdGapVecdv   = ( tPertGapData->mGapVec - tCurrentGapData->mGapVec ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdGapvec2dv2 = ( tPertGapData->mdGapvecdv - tCurrentGapData->mdGapvecdv ) / 2 / tPerturbation;
                        Matrix< DDRMat > tFdGapvec2duv = ( tPertGapData->mdGapvecdu - tCurrentGapData->mdGapvecdu ) / 2 / tPerturbation;

                        real tErrordGapVectordv = 100.0 * norm( tFdGapVecdv - tNominalGapData->mdGapvecdv.get_column( tCounter ) ) / norm( tFdGapVecdv );
                        real tErrordGapvec2dv2  = 100.0 * norm( tFdGapvec2dv2 - tNominalGapData->mdGapvec2dv2( tRowBlk, tColBlk ) ) / norm( tFdGapvec2dv2 );

                        Matrix< DDRMat > tdGapvec2duvExtraction( iSpaceDim, tNumDofs );
                        for ( uint id = 0; id < tNumDofs; id++ )
                        {
                            tdGapvec2duvExtraction.get_column( id ) = tNominalGapData->mdGapvec2duv.get_column( in + idir * tNumNodes + id * tNumDofs );
                        }

                        real tErrordGapvec2duv = 100.0 * norm( tFdGapvec2duv - tdGapvec2duvExtraction ) / norm( tFdGapvec2duv );

                        if ( tErrordGapVectordv + tErrordGapvec2dv2 + tErrordGapvec2duv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error dGapVecdv:   %e\n", tErrordGapVectordv );
                            fprintf( stdout, "percent error dGapvec2dv2: %e\n", tErrordGapvec2dv2 );
                            fprintf( stdout, "percent error dGapvec2duv: %e\n", tErrordGapvec2duv );
                            tPassedCheck = false;
                        }

                        tFollowerDOFHatDisp( in, idir ) += tPerturbation;
                        tCounter++;
                    }
                }

                // clean up
                tLeaderFIs.clear();
                tFollowerFIs.clear();

                // check if the test passed
                REQUIRE( tPassedCheck );
            }
        }
    }
}

TEST_CASE( "IWG_Struc_Contact_Remap_Deformed_Linear", "[moris],[fem],[IWG_Struc_Contact_Remap_Deformed_Linear]" )
{
    // check with contact
    Test_IWG_Struc_Contact_Remap(
            fem::IWG_Type::STRUC_LINEAR_CONTACT_SYMMETRIC_NITSCHE,
            true,     // use deformed geometry
            false,    // use linear displacement assumption
            1.0 );
}

TEST_CASE( "IWG_Struc_Contact_Remap_Deformed_Nonlinear", "[moris],[fem],[IWG_Struc_Contact_Remap_Deformed_Nonlinear]" )
{
    // check with contact
    Test_IWG_Struc_Contact_Remap(
            fem::IWG_Type::STRUC_LINEAR_CONTACT_SYMMETRIC_NITSCHE,
            true,    // use deformed geometry
            true,    // use consistent deformed geometry
            1.0 );
}

/* END_TEST_CASE */
