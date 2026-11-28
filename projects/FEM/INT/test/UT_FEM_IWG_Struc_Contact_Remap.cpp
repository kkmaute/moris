/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Struc_Linear_Contact_Nitsche.cpp
 *
 */

#include "cl_FEM_Contact_Gap_Data.hpp"
#include "cl_FEM_Contact_Geometry.hpp"
#include <string>
#include <functional>
#include <catch.hpp>
#include <memory>
#include <utility>
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
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
#include "fn_isfinite.hpp"
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
    tIWG->mUseDeformedGeometryForGap = aUseDeformedConfig;

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
            mtk::Interpolation_Order tSideInterpolationOrder = tGIInterpolationOrder;

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

            Geometry_Interpolator tLeaderPreviousGI   = Geometry_Interpolator( tGIRule );
            Geometry_Interpolator tFollowerPreviousGI = Geometry_Interpolator( tGIRule );

            Geometry_Interpolator tLeaderPreviousSideGI   = Geometry_Interpolator( tSideGIRule, mtk::CellShape::GENERAL, true );
            Geometry_Interpolator tFollowerPreviousSideGI = Geometry_Interpolator( tSideGIRule, mtk::CellShape::GENERAL, true );

            // create time coeff tHat
            Matrix< DDRMat > tTHat = { { 0.0 } };
            Matrix< DDRMat > tLeaderXHat;
            Matrix< DDRMat > tFollowerXHat;

            Matrix< DDRMat > tLeaderSideNodesParam;
            Matrix< DDRMat > tFollowerSideNodesParam;

            Matrix< DDRMat > tLeaderDOFHatDisp;
            Matrix< DDRMat > tFollowerDOFHatDisp;

            Matrix< DDRMat > tLeaderPreviousDOFHatDisp;
            Matrix< DDRMat > tFollowerPreviousDOFHatDisp;

            Matrix< DDRMat > tIntegPoints;
            Matrix< DDRMat > tSomeSpaceVector;

            switch ( iSpaceDim )
            {
                case 2:
                {
                    tIntegPoints     = { { 0.123 }, { 0.0 } };
                    tSomeSpaceVector = { { 0.75 }, { -0.15 } };

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

                            tLeaderSideNodesParam   = { { -0.9, -0.3 }, { 0.8, -0.1 }, { 0.3, -0.1 } };
                            tFollowerSideNodesParam = { { -0.8, 0.8 }, { 0.9, 0.5 }, { 0.3, 0.7 } };

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

            // initialize the previous time step dofs
            real dux = 0.15;
            real duy = 0.075;

            tLeaderPreviousDOFHatDisp   = tLeaderDOFHatDisp;
            tFollowerPreviousDOFHatDisp = tFollowerDOFHatDisp;

            tLeaderPreviousDOFHatDisp( 2, 0 ) -= dux;
            tLeaderPreviousDOFHatDisp( 0, 1 ) += duy;
            for ( uint iRow = 0; iRow < tFollowerPreviousDOFHatDisp.n_rows(); ++iRow )
            {
                tFollowerPreviousDOFHatDisp( iRow, 1 ) += duy;
            }

            // set the coefficients xHat, tHat
            tLeaderGI.set_coeff( tLeaderXHat, tTHat );
            tFollowerGI.set_coeff( tFollowerXHat, tTHat );

            tLeaderPreviousGI.set_coeff( tLeaderXHat, tTHat );
            tFollowerPreviousGI.set_coeff( tFollowerXHat, tTHat );

            // set the coefficients for the leader side geometry interpolators from parametric coordinates
            Matrix< DDRMat > tLeaderSideXHat( tLeaderSideNodesParam.n_rows(), iSpaceDim );

            for ( uint iNode = 0; iNode < tLeaderSideNodesParam.n_rows(); iNode++ )
            {
                tLeaderGI.set_space( trans( tLeaderSideNodesParam.get_row( iNode ) ) );
                tLeaderPreviousGI.set_space( trans( tLeaderSideNodesParam.get_row( iNode ) ) );

                tLeaderSideXHat.get_row( iNode ) = tLeaderGI.valx().matrix_data();
            }

            tLeaderSideGI.set_coeff( tLeaderSideXHat, tTHat );
            tLeaderSideGI.set_space_param_coeff( tLeaderSideNodesParam );
            tLeaderSideGI.set_time_param_coeff( tTHat );

            tLeaderPreviousSideGI.set_coeff( tLeaderSideXHat, tTHat );
            tLeaderPreviousSideGI.set_space_param_coeff( tLeaderSideNodesParam );
            tLeaderPreviousSideGI.set_time_param_coeff( tTHat );

            // set the coefficients for the follower side geometry interpolators from parametric coordinates
            Matrix< DDRMat > tFollowerSideXHat( tFollowerSideNodesParam.n_rows(), iSpaceDim );

            for ( uint iNode = 0; iNode < tFollowerSideNodesParam.n_rows(); iNode++ )
            {
                tFollowerGI.set_space( trans( tFollowerSideNodesParam.get_row( iNode ) ) );
                tFollowerPreviousGI.set_space( trans( tFollowerSideNodesParam.get_row( iNode ) ) );

                tFollowerSideXHat.get_row( iNode ) = tFollowerGI.valx().matrix_data();
            }

            tFollowerSideGI.set_coeff( tFollowerSideXHat, tTHat );
            tFollowerSideGI.set_space_param_coeff( tFollowerSideNodesParam );
            tFollowerSideGI.set_time_param_coeff( tTHat );

            tFollowerPreviousSideGI.set_coeff( tFollowerSideXHat, tTHat );
            tFollowerPreviousSideGI.set_space_param_coeff( tFollowerSideNodesParam );
            tFollowerPreviousSideGI.set_time_param_coeff( tTHat );

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

            // create a cell of field interpolators for IWG - Leader
            Vector< Field_Interpolator* > tLeaderFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tLeaderGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderDOFHatDisp );

            // create a cell of field interpolators for IWG - Follower
            Vector< Field_Interpolator* > tFollowerFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tFollowerGI, tDispDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerDOFHatDisp );

            // create a cell of field interpolators for IWG - Leader previous time step
            Vector< Field_Interpolator* > tLeaderPreviousFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tLeaderPreviousFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tLeaderGI, tDispDofTypes( 0 ) );
            tLeaderPreviousFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderPreviousDOFHatDisp );

            // create a cell of field interpolators for IWG - Follower previous time step
            Vector< Field_Interpolator* > tFollowerPreviousFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tFollowerPreviousFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tFollowerGI, tDispDofTypes( 0 ) );
            tFollowerPreviousFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerPreviousDOFHatDisp );

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

            // create a field interpolator manager for previous time step
            Field_Interpolator_Manager tLeaderPreviousFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tFollowerPreviousFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderPreviousFIManager.mFI                       = tLeaderPreviousFIs;
            tLeaderPreviousFIManager.mIPGeometryInterpolator   = &tLeaderPreviousGI;
            tLeaderPreviousFIManager.mIGGeometryInterpolator   = &tLeaderPreviousSideGI;
            tFollowerPreviousFIManager.mFI                     = tFollowerPreviousFIs;
            tFollowerPreviousFIManager.mIPGeometryInterpolator = &tFollowerPreviousGI;
            tFollowerPreviousFIManager.mIGGeometryInterpolator = &tFollowerPreviousSideGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager   = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager = &tFollowerFIManager;

            tIWG->mLeaderFIManager         = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager = &tFollowerFIManager;

            tIWG->mSet->mLeaderPreviousFIManager   = &tLeaderPreviousFIManager;
            tIWG->mSet->mFollowerPreviousFIManager = &tFollowerPreviousFIManager;

            tIWG->mLeaderPreviousFIManager         = &tLeaderPreviousFIManager;
            tIWG->mSet->mFollowerPreviousFIManager = &tFollowerPreviousFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tLeaderFIManager, mtk::Leader_Follower::LEADER );
            tIWG->set_field_interpolator_manager( &tFollowerFIManager, mtk::Leader_Follower::FOLLOWER );

            tIWG->set_field_interpolator_manager_previous_time( &tLeaderPreviousFIManager, mtk::Leader_Follower::LEADER );
            tIWG->set_field_interpolator_manager_previous_time( &tFollowerPreviousFIManager, mtk::Leader_Follower::FOLLOWER );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPointLeader = tIntegPoints.get_column( iGP );

                // set integration point on leader side
                tIWG->mSet->mLeaderFIManager->set_space_time_from_local_IG_point( tParamPointLeader );
                tIWG->mSet->mLeaderPreviousFIManager->set_space_time_from_local_IG_point( tParamPointLeader );

                Matrix< DDRMat > tParamPointfollower = tIWG->remap_nonconformal_rays(
                        tIWG->mUseDeformedGeometryForGap,
                        tDispDofTypes( 0 ),
                        tIWG->mSet->mLeaderFIManager,
                        tIWG->mSet->mFollowerFIManager,
                        tIWG->mSet->mLeaderPreviousFIManager,
                        tIWG->mSet->mFollowerPreviousFIManager,
                        tIWG->mGapData );

                // get reference to current gap data
                const std::unique_ptr< GapData >& tCurrentGapData = tIWG->get_gap_data();

                // save mapping results
                std::unique_ptr< GapData > tNominalGapData = std::make_unique< GapData >();
                tNominalGapData->copy( tCurrentGapData );

                Matrix< DDRMat > tNominaldTangentVectordu = tNominalGapData->multiply_leader_dtangentprojectiondu( tSomeSpaceVector );

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
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mSet->mLeaderPreviousFIManager,
                                tIWG->mSet->mFollowerPreviousFIManager,
                                tIWG->mGapData );
                        tPertGapData->copy( tCurrentGapData );

                        tLeaderDOFHatDisp( in, idir ) -= 2.0 * tPerturbation;
                        tLeaderFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mSet->mLeaderPreviousFIManager,
                                tIWG->mSet->mFollowerPreviousFIManager,
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

                        Matrix< DDRMat > tFollowerdDisplacementsdu = ( tPertGapData->mFollowerDisplacements - tCurrentGapData->mFollowerDisplacements ) / 2 / tPerturbation;

                        real tErrorFollowerdDisplacementsdu = 100.0 * norm( tFollowerdDisplacementsdu - tNominalGapData->mFollowerdDisplacementsdu.get_column( tCounter ) ) / norm( tFollowerdDisplacementsdu );

                        if ( tErrorFollowerdDisplacementsdu > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error FollowerdDisplacementsdu:   %e\n", tErrorFollowerdDisplacementsdu );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFollowerPreviousdDisplacementsdu = ( tPertGapData->mFollowerPreviousDisplacements - tCurrentGapData->mFollowerPreviousDisplacements ) / 2 / tPerturbation;

                        real tErrorFollowerPreviousdDisplacementsdu = 100.0 * norm( tFollowerPreviousdDisplacementsdu - tNominalGapData->mFollowerPreviousdDisplacementsdu.get_column( tCounter ) ) / norm( tFollowerPreviousdDisplacementsdu );

                        if ( tErrorFollowerPreviousdDisplacementsdu > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error FollowerPreviousdDisplacementsdu:   %e\n", tErrorFollowerPreviousdDisplacementsdu );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFDTangentVectordu = 0.5 / tPerturbation * ( tPertGapData->mLeaderTangentProjection - tCurrentGapData->mLeaderTangentProjection ) * tSomeSpaceVector;

                        real tErrorTangentVectordu = 100.0 * norm( tFDTangentVectordu - tNominaldTangentVectordu.get_column( tCounter ) ) / norm( tFDTangentVectordu );

                        if ( tErrorTangentVectordu > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error TangentVectordu:   %e\n", tErrorTangentVectordu );
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
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mSet->mLeaderPreviousFIManager,
                                tIWG->mSet->mFollowerPreviousFIManager,
                                tIWG->mGapData );
                        tPertGapData->copy( tCurrentGapData );

                        tFollowerDOFHatDisp( in, idir ) -= 2.0 * tPerturbation;
                        tFollowerFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mSet->mLeaderPreviousFIManager,
                                tIWG->mSet->mFollowerPreviousFIManager,
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

                        Matrix< DDRMat > tLeaderFdTangentVelocitydv = ( tPertGapData->mLeaderTangentVelocity - tCurrentGapData->mLeaderTangentVelocity ) / 2 / tPerturbation;

                        real tErrorLeaderdTangentVelocitydv = 100.0 * norm( tLeaderFdTangentVelocitydv - tNominalGapData->mLeaderdTangentVelocitydv.get_column( tCounter ) ) / norm( tLeaderFdTangentVelocitydv );

                        if ( tErrorLeaderdTangentVelocitydv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error LeaderdTangentVelocitydv:   %e\n", tErrorLeaderdTangentVelocitydv );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFollowerdDisplacementsdv = ( tPertGapData->mFollowerDisplacements - tCurrentGapData->mFollowerDisplacements ) / 2 / tPerturbation;

                        real tErrorFollowerdDisplacementsdv = 100.0 * norm( tFollowerdDisplacementsdv - tNominalGapData->mFollowerdDisplacementsdv.get_column( tCounter ) ) / norm( tFollowerdDisplacementsdv );

                        if ( tErrorFollowerdDisplacementsdv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error FollowerdDisplacementsdv:   %e\n", tErrorFollowerdDisplacementsdv );
                            tPassedCheck = false;
                        }

                        Matrix< DDRMat > tFollowerPreviousdDisplacementsdv = ( tPertGapData->mFollowerPreviousDisplacements - tCurrentGapData->mFollowerPreviousDisplacements ) / 2 / tPerturbation;

                        real tErrorFollowerPreviousdDisplacementsdv = 100.0 * norm( tFollowerPreviousdDisplacementsdv - tNominalGapData->mFollowerPreviousdDisplacementsdv.get_column( tCounter ) ) / norm( tFollowerPreviousdDisplacementsdv );

                        if ( tErrorFollowerPreviousdDisplacementsdv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error FollowerPreviousdDisplacementsdv:   %e\n", tErrorFollowerPreviousdDisplacementsdv );
                            tPassedCheck = false;
                        }

                        tFollowerDOFHatDisp( in, idir ) += tPerturbation;
                        tCounter++;
                    }
                }
                // reset the leader dof hat displacement
                tFollowerFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerDOFHatDisp );

                // FD wrt displacements of previous time step on leader side
                tCounter  = 0;
                tNumNodes = tLeaderPreviousDOFHatDisp.n_rows();

                for ( uint idir = 0; idir < tLeaderPreviousDOFHatDisp.n_cols(); idir++ )
                {
                    for ( uint in = 0; in < tLeaderPreviousDOFHatDisp.n_rows(); in++ )
                    {
                        tLeaderPreviousDOFHatDisp( in, idir ) += tPerturbation;
                        tLeaderPreviousFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderPreviousDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mSet->mLeaderPreviousFIManager,
                                tIWG->mSet->mFollowerPreviousFIManager,
                                tIWG->mGapData );
                        tPertGapData->copy( tCurrentGapData );

                        tLeaderPreviousDOFHatDisp( in, idir ) -= 2.0 * tPerturbation;
                        tLeaderPreviousFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderPreviousDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mSet->mLeaderPreviousFIManager,
                                tIWG->mSet->mFollowerPreviousFIManager,
                                tIWG->mGapData );

                        Matrix< DDRMat > tLeaderFdTangentVelocitydu = ( tPertGapData->mLeaderTangentVelocity - tCurrentGapData->mLeaderTangentVelocity ) / 2 / tPerturbation;

                        real tErrorLeaderdTangentVelocitydu = 100.0 * norm( tLeaderFdTangentVelocitydu - tNominalGapData->mLeaderdTangentVelocityduprev.get_column( tCounter ) ) / norm( tLeaderFdTangentVelocitydu );

                        if ( tErrorLeaderdTangentVelocitydu > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error LeaderdTangentVelocityduprev:   %e\n", tErrorLeaderdTangentVelocitydu );
                            tPassedCheck = false;
                        }

                        tLeaderPreviousDOFHatDisp( in, idir ) += tPerturbation;
                        tCounter++;
                    }
                }
                // reset the leader dof hat displacement
                tLeaderPreviousFIs( 0 )->set_coeff( aDisplacementScaling * tLeaderPreviousDOFHatDisp );

                // FD wrt displacements on leader side
                tCounter  = 0;
                tNumNodes = tFollowerPreviousDOFHatDisp.n_rows();

                for ( uint idir = 0; idir < tFollowerPreviousDOFHatDisp.n_cols(); idir++ )
                {
                    for ( uint in = 0; in < tFollowerPreviousDOFHatDisp.n_rows(); in++ )
                    {
                        tFollowerPreviousDOFHatDisp( in, idir ) += tPerturbation;
                        tFollowerPreviousFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerPreviousDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mSet->mLeaderPreviousFIManager,
                                tIWG->mSet->mFollowerPreviousFIManager,
                                tIWG->mGapData );
                        tPertGapData->copy( tCurrentGapData );

                        tFollowerPreviousDOFHatDisp( in, idir ) -= 2.0 * tPerturbation;
                        tFollowerPreviousFIs( 0 )->set_coeff( aDisplacementScaling * tFollowerPreviousDOFHatDisp );
                        tIWG->reset_eval_flags();
                        tIWG->remap_nonconformal_rays(
                                tIWG->mUseDeformedGeometryForGap,
                                tDispDofTypes( 0 ),
                                tIWG->mSet->mLeaderFIManager,
                                tIWG->mSet->mFollowerFIManager,
                                tIWG->mSet->mLeaderPreviousFIManager,
                                tIWG->mSet->mFollowerPreviousFIManager,
                                tIWG->mGapData );

                        Matrix< DDRMat > tLeaderFdTangentVelocitydv = ( tPertGapData->mLeaderTangentVelocity - tCurrentGapData->mLeaderTangentVelocity ) / 2 / tPerturbation;

                        real tErrorLeaderdTangentVelocitydv = 100.0 * norm( tLeaderFdTangentVelocitydv - tNominalGapData->mLeaderdTangentVelocitydvprev.get_column( tCounter ) ) / norm( tLeaderFdTangentVelocitydv );

                        if ( tErrorLeaderdTangentVelocitydv > 100.0 * tEpsilon )
                        {
                            fprintf( stdout, "percent error LeaderdTangentVelocitydvprev:   %e\n", tErrorLeaderdTangentVelocitydv );
                            tPassedCheck = false;
                        }

                        tFollowerPreviousDOFHatDisp( in, idir ) += tPerturbation;
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

TEST_CASE( "IWG_Struc_Contact_Remap_Deformed_Nonlinear", "[moris],[fem],[IWG_Struc_Contact_Remap_Deformed_Nonlinear]" )
{
    // check with contact
    Test_IWG_Struc_Contact_Remap(
            fem::IWG_Type::STRUC_LINEAR_CONTACT_SYMMETRIC_NITSCHE,
            true,    // use deformed geometry
            1.0 );
}


namespace
{
    enum class Ray_Face_Type
    {
        COINCIDENT,
        PARALLEL,
        COLLAPSED
    };

    struct Ray_Mapping_Result
    {
        Matrix< DDRMat >           point;
        std::unique_ptr< GapData > gap;
    };

    Ray_Mapping_Result map_straight_contact_ray(
            real                                                                                                          aGeometryScale,
            real                                                                                                          aStretch,
            real                                                                                                          aLeaderCoordinate,
            Ray_Face_Type                                                                                                 aFaceType                  = Ray_Face_Type::COINCIDENT,
            bool                                                                                                          aEvaluateSecondDerivatives = true,
            const std::function< void( Field_Interpolator_Manager&, Field_Interpolator_Manager&, Ray_Mapping_Result& ) >& aCheck                     = {} )
    {
        Set tSet;
        tSet.set_set_type( Element_Type::DOUBLE_SIDESET );
        tSet.mUniqueDofTypeList = { MSI::Dof_Type::UX };
        tSet.mLeaderDofTypeMap.set_size( static_cast< uint >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tSet.mLeaderDofTypeMap( static_cast< uint >( MSI::Dof_Type::UX ) ) = 0;
        tSet.mFollowerDofTypeMap                                           = tSet.mLeaderDofTypeMap;

        const mtk::Interpolation_Rule tCellRule( mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::CONSTANT,
                mtk::Interpolation_Order::CONSTANT );
        const mtk::Interpolation_Rule tSideRule( mtk::Geometry_Type::LINE,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::CONSTANT,
                mtk::Interpolation_Order::CONSTANT );
        Geometry_Interpolator         tLeaderGI( tCellRule );
        Geometry_Interpolator         tFollowerGI( tCellRule );
        Geometry_Interpolator         tLeaderSideGI( tSideRule, mtk::CellShape::GENERAL, true );
        Geometry_Interpolator         tFollowerSideGI( tSideRule, mtk::CellShape::GENERAL, true );

        Matrix< DDRMat > tLeaderCoordinates      = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
        Matrix< DDRMat > tFollowerCoordinates    = { { 0.0, -1.0 }, { 1.0, -1.0 }, { 1.0, 0.0 }, { 0.0, 0.0 } };
        Matrix< DDRMat > tFollowerSideParameters = { { 1.0, 1.0 }, { -1.0, 1.0 } };
        if ( aFaceType == Ray_Face_Type::PARALLEL )
        {
            // A vertical follower face passes through the leader midpoint.
            // Its tangent is parallel to the ray, making the ray Jacobian
            // singular even though the initial residual is exactly zero.
            tFollowerCoordinates    = { { 0.0, -0.5 }, { 1.0, -0.5 }, { 1.0, 0.5 }, { 0.0, 0.5 } };
            tFollowerSideParameters = { { 0.0, -1.0 }, { 0.0, 1.0 } };
        }
        tLeaderCoordinates           = aGeometryScale * tLeaderCoordinates;
        tFollowerCoordinates         = aGeometryScale * tFollowerCoordinates;
        const Matrix< DDRMat > tTime = { { 0.0 } };
        tLeaderGI.set_coeff( tLeaderCoordinates, tTime );
        tFollowerGI.set_coeff( tFollowerCoordinates, tTime );

        auto tSetSide = [ & ]( Geometry_Interpolator& aCell, Geometry_Interpolator& aSide, const Matrix< DDRMat >& aParameters ) {
            Matrix< DDRMat > tCoordinates( aParameters.n_rows(), 2 );
            for ( uint i = 0; i < aParameters.n_rows(); ++i )
            {
                aCell.set_space( trans( aParameters.get_row( i ) ) );
                tCoordinates.get_row( i ) = aCell.valx().matrix_data();
            }
            aSide.set_coeff( tCoordinates, tTime );
            aSide.set_space_param_coeff( aParameters );
            aSide.set_time_param_coeff( tTime );
        };
        tSetSide( tLeaderGI, tLeaderSideGI, Matrix< DDRMat >{ { -1.0, -1.0 }, { 1.0, -1.0 } } );
        tSetSide( tFollowerGI, tFollowerSideGI, tFollowerSideParameters );

        // Managers own the heap field interpolators and borrow the stack
        // geometries. The set is not given ownership of these stack managers.
        const Vector< Vector< MSI::Dof_Type > > tDofTypes = { { MSI::Dof_Type::UX } };
        Field_Interpolator_Manager              tLeaderManager( tDofTypes, &tSet );
        Field_Interpolator_Manager              tFollowerManager( tDofTypes, &tSet, mtk::Leader_Follower::FOLLOWER );
        auto                                    tLeaderFI   = std::make_unique< Field_Interpolator >( 2, tCellRule, &tLeaderGI, tDofTypes( 0 ) );
        auto                                    tFollowerFI = std::make_unique< Field_Interpolator >( 2, tCellRule, &tFollowerGI, tDofTypes( 0 ) );
        Matrix< DDRMat >                        tLeaderDisplacements( 4, 2, 0.0 );
        Matrix< DDRMat >                        tFollowerDisplacements( 4, 2, 0.0 );
        for ( uint i = 0; i < 4; ++i )
        {
            tLeaderDisplacements( i, 0 )   = ( aStretch - 1.0 ) * tLeaderCoordinates( i, 0 );
            tFollowerDisplacements( i, 0 ) = ( aStretch - 1.0 ) * tFollowerCoordinates( i, 0 );
        }
        if ( aFaceType == Ray_Face_Type::COLLAPSED )
        {
            // Collapse the deformed follower face onto the leader midpoint.
            // Its reference geometry stays valid, so rejection must come from
            // the current ray Jacobian rather than fixture initialization.
            for ( uint i = 0; i < 4; ++i )
            {
                tFollowerDisplacements( i, 0 ) = 0.5 * aStretch * aGeometryScale - tFollowerCoordinates( i, 0 );
            }
        }
        tLeaderFI->set_coeff( tLeaderDisplacements );
        tFollowerFI->set_coeff( tFollowerDisplacements );
        tLeaderManager.mFI                       = { tLeaderFI.release() };
        tFollowerManager.mFI                     = { tFollowerFI.release() };
        tLeaderManager.mIPGeometryInterpolator   = &tLeaderGI;
        tLeaderManager.mIGGeometryInterpolator   = &tLeaderSideGI;
        tFollowerManager.mIPGeometryInterpolator = &tFollowerGI;
        tFollowerManager.mIGGeometryInterpolator = &tFollowerSideGI;
        tLeaderManager.set_space_time_from_local_IG_point( Matrix< DDRMat >{ { aLeaderCoordinate }, { 0.0 } } );
        tFollowerManager.set_space_time_from_local_IG_point( Matrix< DDRMat >{ { 0.0 }, { 0.0 } } );

        // Exercise mapping independently of an IWG; the caller owns this cache.
        Ray_Mapping_Result tResult;
        tResult.point = Contact_Geometry::remap_nonconformal_rays( true, tDofTypes( 0 ), &tLeaderManager, &tFollowerManager, nullptr, nullptr, tResult.gap, aEvaluateSecondDerivatives );
        if ( aCheck )
        {
            aCheck( tLeaderManager, tFollowerManager, tResult );
        }
        return tResult;
    }
}    // namespace

// Coincident faces with opposite normals have an exact initial intersection.
// Check that it is accepted and that the inverse needed by gap derivatives is
// available, independently of geometric scale and an applied uniform stretch.
TEST_CASE( "Contact ray accepts exact initial intersections", "[fem][ContactRayConvergence]" )
{
    for ( real tScale : { 1.0e-6, 1.0, 1.0e6 } )
    {
        for ( real tStretch : { 1.0, 32.0 } )
        {
            INFO( "Geometry scale=" << tScale << " stretch=" << tStretch );
            auto tResult = map_straight_contact_ray( tScale, tStretch, 0.0 );
            CHECK( tResult.point( 0 ) == 0.0 );
            CHECK( tResult.gap->mGap == Approx( 0.0 ).margin( 1.0e-12 * tScale * tStretch ) );
            REQUIRE( tResult.gap->mdGapdu.numel() > 0 );
            CHECK( isfinite( tResult.gap->mdGapdu ) );
        }
    }
}

// A tiny nonzero residual is accepted using the deformed face length. The
// 32-fold stretch makes it exceed the reference-length tolerance, so this also
// checks that the length scale includes displacement. An off-center root still
// requires Newton updates.
TEST_CASE( "Contact ray scales its absolute residual tolerance", "[fem][ContactRayConvergence]" )
{
    for ( real tScale : { 1.0e-6, 1.0, 1.0e6 } )
    {
        for ( real tStretch : { 1.0, 32.0 } )
        {
            INFO( "Geometry scale=" << tScale << " stretch=" << tStretch );
            auto tNear = map_straight_contact_ray( tScale, tStretch, 1.0e-13 );
            CHECK( tNear.point( 0 ) == 0.0 );
            auto tOffCenter = map_straight_contact_ray( tScale, tStretch, 0.123 );
            CHECK( tOffCenter.point( 0 ) == Approx( -0.123 ).margin( 1.0e-12 ) );
        }
    }
}

// A zero residual alone must not accept an underdetermined intersection.
// The Jacobian check must still reject parallel rays and collapsed faces.
TEST_CASE( "Contact ray rejects singular exact intersections", "[fem][ContactRayConvergence]" )
{
    for ( real tScale : { 1.0e-6, 1.0, 1.0e6 } )
    {
        INFO( "Geometry scale=" << tScale );
        for ( auto tFaceType : { Ray_Face_Type::PARALLEL, Ray_Face_Type::COLLAPSED } )
        {
            INFO( "Face type=" << static_cast< uint >( tFaceType ) );
            auto tResult = map_straight_contact_ray( tScale, 1.0, 0.0, tFaceType );
            CHECK( tResult.point( 0 ) == -2.0 );
            CHECK( tResult.gap->mEta( 0 ) == -2.0 );
        }
    }
}

// Residual requests must avoid every Hessian buffer. Upgrading that same
// point for a Jacobian must agree with a fresh full evaluation, and edits or
// failed intersections must never expose Hessians from the previous geometry.
TEST_CASE( "Contact geometry evaluates Hessians on request", "[fem][ContactHessianEvaluation]" )
{
    using Gap_Matrix                                      = Matrix< DDRMat >       GapData::*;
    const std::initializer_list< Gap_Matrix > tFirstOrder = {
        &GapData::mEta, &GapData::mdGapdu, &GapData::mdGapdv, &GapData::mdEtadu, &GapData::mdEtadv, &GapData::mLeaderNormal, &GapData::mLeaderRefNormal, &GapData::mLeaderdNormaldu, &GapData::mGapVec, &GapData::mdGapvecdu, &GapData::mdGapvecdv, &GapData::mFollowerDisplacements, &GapData::mFollowerdDisplacementsdu, &GapData::mFollowerdDisplacementsdv, &GapData::mFollowerShapeFunctiondu, &GapData::mFollowerShapeFunctiondv
    };
    const std::initializer_list< Gap_Matrix > tHessians = {
        &GapData::mdGap2du2, &GapData::mdGap2dv2, &GapData::mdGap2duv, &GapData::mdEta2du2, &GapData::mdEta2dv2, &GapData::mdEta2duv, &GapData::mLeaderdNormal2du2, &GapData::mdGapvec2du2, &GapData::mdGapvec2dv2, &GapData::mdGapvec2duv
    };

    auto tCheckEmpty = [ & ]( const GapData& aGap ) {
        REQUIRE_FALSE( aGap.mHasSecondOrderDerivatives );
        for ( auto tMember : tHessians )
        {
            REQUIRE( ( aGap.*tMember ).numel() == 0 );
        }
    };
    auto tCompare = [ & ]( const GapData& aActual, const GapData& aReference, bool aSecondOrder ) {
        CHECK( aActual.mGap == Approx( aReference.mGap ).margin( 1.0e-12 ) );
        auto tCompareMatrices = [ & ]( const std::initializer_list< Gap_Matrix >& aMembers ) {
            for ( auto tMember : aMembers )
            {
                const auto& tActual    = aActual.*tMember;
                const auto& tReference = aReference.*tMember;
                REQUIRE( tActual.n_rows() == tReference.n_rows() );
                REQUIRE( tActual.n_cols() == tReference.n_cols() );
                CHECK( norm( tActual - tReference ) <= 1.0e-12 * ( 1.0 + norm( tReference ) ) );
            }
        };
        tCompareMatrices( tFirstOrder );
        if ( aSecondOrder )
        {
            REQUIRE( aActual.mHasSecondOrderDerivatives );
            tCompareMatrices( tHessians );
        }
    };

    map_straight_contact_ray( 1.0, 1.0, 0.123, Ray_Face_Type::COINCIDENT, false, [ & ]( Field_Interpolator_Manager& aLeader, Field_Interpolator_Manager& aFollower, Ray_Mapping_Result& aResult ) {
        auto tMap = [ & ]( std::unique_ptr< GapData >& aGap, bool aSecondOrder ) {
            return Contact_Geometry::remap_nonconformal_rays( true, { MSI::Dof_Type::UX }, &aLeader, &aFollower, nullptr, nullptr, aGap, aSecondOrder );
        };

        tCheckEmpty( *aResult.gap );
        REQUIRE_FALSE( aResult.gap->mEval );
        std::unique_ptr< GapData > tFull;
        tMap( tFull, true );
        REQUIRE( norm( tFull->mLeaderdNormal2du2 ) > 1.0e-6 );
        tCompare( *aResult.gap, *tFull, false );

        // A residual-to-Jacobian upgrade preserves locked contact state.
        aResult.gap->mIsLocked  = true;
        aResult.gap->mIsContact = true;
        tMap( aResult.gap, true );
        tCompare( *aResult.gap, *tFull, true );
        CHECK( aResult.gap->mIsLocked );
        CHECK( aResult.gap->mIsContact );
        tMap( aResult.gap, false );
        tCompare( *aResult.gap, *tFull, true );

        // Edit both the leader cell and its extracted side. The owning
        // caller invalidates the contact cache after changing geometry.
        for ( auto tGeometry : { aLeader.get_IP_geometry_interpolator(), aLeader.get_IG_geometry_interpolator() } )
        {
            Matrix< DDRMat > tCoordinates = tGeometry->get_space_coeff();
            for ( uint i = 0; i < tCoordinates.n_rows(); ++i )
            {
                tCoordinates( i, 1 ) -= 0.05;
            }
            tGeometry->set_space_coeff( tCoordinates );
        }
        aLeader.set_space_time_from_local_IG_point( Matrix< DDRMat >{ { 0.123 }, { 0.0 } } );
        aResult.gap->mEval = true;
        tMap( aResult.gap, false );
        tCheckEmpty( *aResult.gap );
        CHECK( std::abs( aResult.gap->mGap ) > 0.01 );
        tFull.reset();
        tMap( tFull, true );
        tCompare( *aResult.gap, *tFull, false );
        tMap( aResult.gap, true );
        tCompare( *aResult.gap, *tFull, true );

        // Collapse the follower onto a point. A failed mapping clears
        // old Hessians, and requesting them again retains the sentinel.
        auto*            tFollowerFI    = aFollower.get_field_interpolators_for_type( MSI::Dof_Type::UX );
        Matrix< DDRMat > tDisplacements = tFollowerFI->get_coeff();
        const auto&      tCoordinates   = aFollower.get_IP_geometry_interpolator()->get_space_coeff();
        for ( uint i = 0; i < tCoordinates.n_rows(); ++i )
        {
            tDisplacements( i, 0 ) = 0.5 - tCoordinates( i, 0 );
        }
        tFollowerFI->set_coeff( tDisplacements );
        aResult.gap->mEval = true;
        CHECK( tMap( aResult.gap, true )( 0 ) == -2.0 );
        tCheckEmpty( *aResult.gap );
        CHECK( tMap( aResult.gap, true )( 0 ) == -2.0 );
        tCheckEmpty( *aResult.gap );
    } );
}

/* END_TEST_CASE */
