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

void Test_IWG_Struc_Nonlinear_Contact_Nitsche(
        const IWG_Type               aIWGType,
        const fem::Constitutive_Type aConstitutiveModel,
        const real&                  aNitscheParameter,
        bool                         aUseDeformedGeometryForGap,
        bool                         aUseConsistentDeformedGeometryForGap )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

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

    // Nitsche parameters: used to trigger in and out of contact condition
    Vector< real > tNitscheParameter = { 1.0, 100 };

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

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStruct =
            tCMFactory.create_CM( aConstitutiveModel );
    tCMLeaderStruct->set_dof_type_list( tDispDofTypes );
    tCMLeaderStruct->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStruct->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStruct->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerStruct =
            tCMFactory.create_CM( aConstitutiveModel );
    tCMFollowerStruct->set_dof_type_list( tDispDofTypes );
    tCMFollowerStruct->set_property( tPropFollowerEMod, "YoungsModulus" );
    tCMFollowerStruct->set_property( tPropFollowerNu, "PoissonRatio" );
    tCMFollowerStruct->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { { { aNitscheParameter } } } );
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
    tIWG->set_constitutive_model( tCMLeaderStruct, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerStruct, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::NONCONFORMAL_SIDESET );
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

    // set parameter for gap computation
    tIWG->mUseDeformedGeometryForGap           = aUseDeformedGeometryForGap;
    tIWG->mUseConsistentDeformedGeometryForGap = aUseConsistentDeformedGeometryForGap;

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
        tCMLeaderStruct->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMLeaderStruct->set_space_dim( iSpaceDim );

        tCMFollowerStruct->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMFollowerStruct->set_space_dim( iSpaceDim );

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
            mtk::Interpolation_Order tGIInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // set side interpolation order to linear if inconsistent deformed geometry is used
            mtk::Interpolation_Order tSideInterpolationOrder =
                    aUseConsistentDeformedGeometryForGap ? tGIInterpolationOrder : tInterpolationOrders( 0 );

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

                            if ( aUseConsistentDeformedGeometryForGap )
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
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDisp );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tFollowerGI, tDispDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatDisp );

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

                Matrix< DDRMat > tParamPointfollower = tIWG->remap_nonconformal_rays( tLeaderFIs( 0 ), tFollowerFIs( 0 ) );

                tIWG->mSet->mFollowerFIManager->set_space_time_from_local_IG_point( tParamPointfollower );

                tIWG->reset_eval_flags();

                // check evaluation of the residual for IWG
                //------------------------------------------------------------------------------
                // reset residual
                tIWG->mSet->mResidual( 0 )
                        .fill( 0.0 );

                // compute residual
                tIWG->compute_residual( 1.0 );

                // check evaluation of the jacobian by FD
                //------------------------------------------------------------------------------
                // reset jacobian
                tIWG->mSet->mJacobian.fill( 0.0 );

                // compute jacobian
                tIWG->compute_jacobian( 1.0 );

                // save analytic jacobian
                Matrix< DDRMat > tJacobian = tIWG->mSet->mJacobian;

                // init the jacobian for IWG and FD evaluation
                Matrix< DDRMat > tJacobianFD( 2 * tNumDofDisp, 2 * tNumDofDisp, 0.0 );
                tIWG->mSet->mJacobian.fill( 0.0 );

                uint tCounter = 0;
                for ( uint idir = 0; idir < tLeaderDOFHatDisp.n_cols(); idir++ )
                {
                    for ( uint in = 0; in < tLeaderDOFHatDisp.n_rows(); in++ )
                    {
                        tLeaderDOFHatDisp( in, idir ) += tPerturbation;
                        tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDisp );

                        tIWG->reset_eval_flags();
                        tIWG->mSet->mResidual( 0 ).fill( 0.0 );
                        tIWG->compute_residual( 1.0 );

                        tJacobianFD.get_column( tCounter ) = 0.5 / tPerturbation * tIWG->mSet->mResidual( 0 );

                        tLeaderDOFHatDisp( in, idir ) -= 2.0 * tPerturbation;
                        tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDisp );

                        tIWG->reset_eval_flags();
                        tIWG->mSet->mResidual( 0 ).fill( 0.0 );
                        tIWG->compute_residual( 1.0 );

                        tJacobianFD.get_column( tCounter ) -= 0.5 / tPerturbation * tIWG->mSet->mResidual( 0 );

                        tLeaderDOFHatDisp( in, idir ) += tPerturbation;
                        tCounter++;
                    }
                }

                // reset the leader dof hat displacement
                tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDisp );

                // FD wrt displacements on leader side
                for ( uint idir = 0; idir < tFollowerDOFHatDisp.n_cols(); idir++ )
                {
                    for ( uint in = 0; in < tFollowerDOFHatDisp.n_rows(); in++ )
                    {
                        tFollowerDOFHatDisp( in, idir ) += tPerturbation;
                        tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatDisp );

                        tIWG->reset_eval_flags();
                        tIWG->mSet->mResidual( 0 ).fill( 0.0 );
                        tIWG->compute_residual( 1.0 );

                        tJacobianFD.get_column( tCounter ) = 0.5 / tPerturbation * tIWG->mSet->mResidual( 0 );

                        tFollowerDOFHatDisp( in, idir ) -= 2.0 * tPerturbation;
                        tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatDisp );

                        tIWG->reset_eval_flags();
                        tIWG->mSet->mResidual( 0 ).fill( 0.0 );
                        tIWG->compute_residual( 1.0 );

                        tJacobianFD.get_column( tCounter ) -= 0.5 / tPerturbation * tIWG->mSet->mResidual( 0 );

                        tFollowerDOFHatDisp( in, idir ) += tPerturbation;
                        tCounter++;
                    }
                }

                // reset the follower dof hat displacement
                tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatDisp );
                tIWG->reset_eval_flags();

                real errorNorm = norm( tJacobian - tJacobianFD ) / norm( tJacobianFD );

                if ( errorNorm > tEpsilon )
                {
                    tPassedCheck = false;
                    fprintf( stdout, "Jacobian relative error norm = %e\n", errorNorm );
                    print( tJacobian, "AnalyticalJacobian" );
                    print( tJacobianFD, "FDJacobian" );
                    print( tJacobian - tJacobianFD, "Jacobian error" );
                }

                // check jacobian by FD with standard method

                bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                        tEpsilon,
                        1.0,
                        tJacobian,
                        tJacobianFD,
                        true,
                        true );

                if ( !tCheckJacobian )
                {
                    errorNorm = norm( tJacobian - tJacobianFD ) / norm( tJacobianFD );
                    fprintf( stdout, "Jacobian relative error norm = %e\n", errorNorm );
                    print( tJacobian, "AnalyticalJacobian" );
                    print( tJacobianFD, "FDJacobian" );
                    print( tJacobian - tJacobianFD, "Jacobian error" );
                }

                // clean up
                tLeaderFIs.clear();
                tFollowerFIs.clear();

                // check if the test passed
                REQUIRE( tPassedCheck );
                REQUIRE( tCheckJacobian );
            }
        }
    }
}

//----------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_Linear_Contact_Nitsche_Unbiased_SYM", "[moris],[fem],[IWG_Struc_Linear_Contact_Nitsche_Unbiased_SYM]" )
{
    // check with contact
    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_SYMMETRIC_NITSCHE_UNBIASED,
            // fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            fem::Constitutive_Type::STRUC_LIN_ISO,
            10.0,
            false,     // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation

    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_SYMMETRIC_NITSCHE_UNBIASED,
            // fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            fem::Constitutive_Type::STRUC_LIN_ISO,
            0.01,
            false,     // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation
}

//----------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_Nonlinear_Contact_Nitsche_NEOHOOK_MLIKA_SYM", "[moris],[fem],[IWG_Struc_Nonlinear_Contact_Nitsche_NEOHOOK_MLIKA_SYM]" )
{
    // check with contact
    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_SYMMETRIC,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            10.0,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation

    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_SYMMETRIC,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            0.01,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation
}

//----------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_Nonlinear_Contact_Nitsche_NEOHOOK_MLIKA_UNSYM", "[moris],[fem],[IWG_Struc_Nonlinear_Contact_Nitsche_NEOHOOK_MLIKA_UNSYM]" )
{
    // check with contact
    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_UNSYMMETRIC,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            10.0,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation

    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_UNSYMMETRIC,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            0.01,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation
}

//----------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_Nonlinear_Contact_Nitsche_NEOHOOK_MLIKA_NEUTRAL", "[moris],[fem],[IWG_Struc_Nonlinear_Contact_Nitsche_NEOHOOK_MLIKA_NEUTRAL]" )
{
    // check with contact
    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_NEUTRAL,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            10.0,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation

    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_NEUTRAL,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            0.01,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation
}

//----------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_Nonlinear_Contact_Nitsche_SVK_MLIKA_SYM", "[moris],[fem],[IWG_Struc_Nonlinear_Contact_Nitsche_SVK_MLIKA_SYM]" )
{
    // check with contact
    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_SYMMETRIC,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF,
            10.0,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation

    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_SYMMETRIC,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF,
            0.01,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation
}

//----------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_Nonlinear_Contact_Nitsche_Linearelastic_MLIKA_SYM", "[moris],[fem],[IWG_Struc_Nonlinear_Contact_Nitsche_Linearelastic_MLIKA_SYM]" )
{
    // check with contact
    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_LINEAR_UNBIASED_SYMMETRIC,
            fem::Constitutive_Type::STRUC_LIN_ISO,
            10.0,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation

    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_LINEAR_UNBIASED_SYMMETRIC,
            fem::Constitutive_Type::STRUC_LIN_ISO,
            0.01,
            true,      // compute gap on undeformed configuration
            true );    // compute gap consistent with displacement interpolation
}

//----------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_Nonlinear_Contact_Nitsche_NEOHOOK_MLIKA_SYM_LinDisp", "[moris],[fem],[IWG_Struc_Nonlinear_Contact_Nitsche_NEOHOOK_MLIKA_SYM_LinDisp]" )
{

    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_SYMMETRIC,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            0.01,
            true,       // compute gap on undeformed configuration
            false );    // compute gap using linear displacement along contact interface

    // check with contact
    Test_IWG_Struc_Nonlinear_Contact_Nitsche(
            IWG_Type::STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_SYMMETRIC,
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
            10.0,
            true,       // compute gap on undeformed configuration
            false );    // compute gap using linear displacement along contact interface
}

//----------------------------------------------------------------------------------------------------------------------------

/* END_TEST_CASE */
