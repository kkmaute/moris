/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Struc_Linear_Interface.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM//INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
//LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Struc_Linear_Interface", "[moris],[fem],[IWG_Struc_Linear_Interface]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    Vector< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< Vector< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property > ();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property > ();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropFollowerEMod = std::make_shared< fem::Property > ();
    tPropFollowerEMod->set_parameters( {{{ 20.0 }}} );
    tPropFollowerEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropFollowerNu = std::make_shared< fem::Property > ();
    tPropFollowerNu->set_parameters( {{{ 0.3 }}} );
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
    tSPNitscheInterface->set_parameters( { {{ 100.0 }} } );
    tSPNitscheInterface->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );
    tSPNitscheInterface->set_property( tPropFollowerEMod, "Material", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitscheInterface->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerStrucLinIso, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
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
    tIWG->mSet->mFollowerDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
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
        switch( iSpaceDim )
        {
            case 2 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;
                break;
            }
            case 3 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // create an interpolation order
            mtk::Interpolation_Order tGIInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // space and time geometry interpolators
            //------------------------------------------------------------------------------
            // create a space geometry interpolation rule
            mtk::Interpolation_Rule tGIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tGIInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create a space time geometry interpolator
            Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

            // create time coeff tHat
            Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

            Matrix< DDRMat > tXHat;
            fill_xhat_Elast( tXHat, iSpaceDim, iInterpOrder );

            // set the coefficients xHat, tHat
            tGI.set_coeff( tXHat, tTHat );

            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder = tIntegrationOrders( iSpaceDim - 2 );

            // create an integration rule
            mtk::Integration_Rule tIntegrationRule(
                    tGeometryType,
                    mtk::Integration_Type::GAUSS,
                    tIntegrationOrder,
                    mtk::Geometry_Type::LINE,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );

            // create an integrator
            mtk::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // number of dof for interpolation order
            uint tNumCoeff = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // get number of dof per type
            int tNumDofDisp  = tNumCoeff * iSpaceDim;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatDisp;
            fill_uhat_Elast( tLeaderDOFHatDisp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDisp );

            // fill random coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatDisp;
            fill_uhat_Elast( tFollowerDOFHatDisp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDispDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, 2 * tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofDisp - 1 },
                    { tNumDofDisp, 2 * tNumDofDisp - 1 } };
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
            tIWG->mRequestedLeaderGlobalDofTypes = tDispDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDispDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > > tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tLeaderFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tFollowerFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderFIManager.mFI = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator = &tGI;
            tLeaderFIManager.mIGGeometryInterpolator = &tGI;
            tFollowerFIManager.mFI = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator = &tGI;
            tFollowerFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager = &tFollowerFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tLeaderFIManager, mtk::Leader_Follower::LEADER );
            tIWG->set_field_interpolator_manager( &tFollowerFIManager, mtk::Leader_Follower::FOLLOWER );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mFollowerFIManager->set_space_time( tParamPoint );

                // check evaluation of the residual for IWG
                //------------------------------------------------------------------------------
                // reset residual
                tIWG->mSet->mResidual( 0 ).fill( 0.0 );

                // compute residual
                tIWG->compute_residual( 1.0 );

                // check evaluation of the jacobian by FD
                //------------------------------------------------------------------------------
                // reset jacobian
                tIWG->mSet->mJacobian.fill( 0.0 );

                // init the jacobian for IWG and FD evaluation
                Matrix< DDRMat > tJacobian;
                Matrix< DDRMat > tJacobianFD;

                // check jacobian by FD
                bool tCheckJacobian = tIWG->check_jacobian(
                        tPerturbation,
                        tEpsilon,
                        1.0,
                        tJacobian,
                        tJacobianFD,
                        true );

                // print for debug
                if( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << '\n';
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Linear_Elasticity_Damage_Interface",
        "[moris],[fem],[IWG_Struc_Linear_Elasticity_Damage_Interface]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 5E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-5;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    Vector< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
        Vector< MSI::Dof_Type > tDisplDofTypes = { MSI::Dof_Type::UX };

        Vector< Vector< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
        Vector< Vector< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
        Vector< Vector< MSI::Dof_Type > > tDofTypes           = { tDisplDofTypes, tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // Damage setup - 0,2,0
    // Local equivalent strain - 0 - energy release rate
    // Damage law - 2 - smooth exponential
    // Smoothing law - 0 - no smoothing
    Vector< Matrix< DDRMat > > tDamageParameters = { //
            { { 0.0 } },                                      //
            { { 2.0, 1.0e-3, 10.0 } },                        //
            { { 0.0 } }
    };

    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property > ();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property > ();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropFollowerEMod = std::make_shared< fem::Property > ();
    tPropFollowerEMod->set_parameters( {{{ 20.0 }}} );
    tPropFollowerEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropFollowerNu = std::make_shared< fem::Property > ();
    tPropFollowerNu->set_parameters( {{{ 0.3 }}} );
    tPropFollowerNu->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
                tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeaderStrucLinIso->set_dof_type_list( tDofTypes );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();
    tCMLeaderStrucLinIso->set_parameters( tDamageParameters );

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerStrucLinIso =
                    tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMFollowerStrucLinIso->set_dof_type_list( tDofTypes );
    tCMFollowerStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMFollowerStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMFollowerStrucLinIso->set_local_properties();
    tCMFollowerStrucLinIso->set_parameters( tDamageParameters );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { {{ 100.0 }} } );
    tSPNitscheInterface->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );
    tSPNitscheInterface->set_property( tPropFollowerEMod, "Material", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitscheInterface->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( { tDisplDofTypes } );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerStrucLinIso, "ElastLinIso", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
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
        switch( iSpaceDim )
        {
            case 2 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;
                break;
            }
            case 3 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // create an interpolation order
            mtk::Interpolation_Order tGIInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // space and time geometry interpolators
            //------------------------------------------------------------------------------
            // create a space geometry interpolation rule
            mtk::Interpolation_Rule tGIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tGIInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create a space time geometry interpolator
            Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

            // create time coeff tHat
            Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

            Matrix< DDRMat > tXHat;
            fill_xhat_Elast( tXHat, iSpaceDim, iInterpOrder );

            // set the coefficients xHat, tHat
            tGI.set_coeff( tXHat, tTHat );

            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder = tIntegrationOrders( iSpaceDim - 2 );

            // create an integration rule
            mtk::Integration_Rule tIntegrationRule(
                    tGeometryType,
                    mtk::Integration_Type::GAUSS,
                    tIntegrationOrder,
                    mtk::Geometry_Type::LINE,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );

            // create an integrator
            mtk::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // number of dof for interpolation order
            uint tNumCoeff = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // get number of dof per type
            int tNumDofDispl      = tNumCoeff * iSpaceDim;
            int tNumDofNlEqStrain = tNumCoeff;
            int tNumDofHistory    = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FIs
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 1e-3*tLeaderDOFHatNlEqStrain;
            Matrix< DDRMat > tLeaderDOFHatHistory;
            fill_phat_Elast( tLeaderDOFHatHistory, iSpaceDim, iInterpOrder );
            tLeaderDOFHatHistory = 1e-3*tLeaderDOFHatHistory;

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDisplDofTypes );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatHistory );

            // fill coefficients for leader FIs
            Matrix< DDRMat > tFollowerDOFHatDisp;
            fill_uhat_Elast( tFollowerDOFHatDisp, iSpaceDim, iInterpOrder );
            tFollowerDOFHatDisp = 0.88 * tFollowerDOFHatDisp;
            Matrix< DDRMat > tFollowerDOFHatNlEqStrain;
            fill_phat_Elast( tFollowerDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tFollowerDOFHatNlEqStrain = 0.77 * 1e-3* tFollowerDOFHatNlEqStrain;
            Matrix< DDRMat > tFollowerDOFHatHistory;
            fill_phat_Elast( tFollowerDOFHatHistory, iSpaceDim, iInterpOrder );
            tFollowerDOFHatHistory = 0.66 *1e-3* tFollowerDOFHatHistory;

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator displacement
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDisplDofTypes );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatDisp );

            // create the field interpolator nonlocal equivalent strain
            tFollowerFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            tFollowerFIs( 1 )->set_coeff( tFollowerDOFHatNlEqStrain );

            // create the field interpolator nonlocal equivalent strain
            tFollowerFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            tFollowerFIs( 2 )->set_coeff( tFollowerDOFHatHistory );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDispl - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDispl, tNumDofNlEqStrain + tNumDofDispl - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofDispl + tNumDofNlEqStrain, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, 2 * tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, 2 * ( tNumDofDispl + tNumDofNlEqStrain ) + tNumDofHistory - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofDispl + tNumDofNlEqStrain ) + tNumDofHistory, 2 * ( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofDispl-1 },
                    { tNumDofDispl, tNumDofDispl + tNumDofNlEqStrain - 1 },
                    { tNumDofDispl + tNumDofNlEqStrain, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 },
                    { tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, 2 * tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 },
                    { 2 * tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, 2 * ( tNumDofDispl + tNumDofNlEqStrain ) + tNumDofHistory - 1 },
                    { 2 * ( tNumDofDispl + tNumDofNlEqStrain ) + tNumDofHistory, 2 * ( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory ) - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 4 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 5 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    2 * ( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory ),
                    2 * ( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory ),
                    0.0 );
            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > > tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tLeaderFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tFollowerFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderFIManager.mFI = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator = &tGI;
            tLeaderFIManager.mIGGeometryInterpolator = &tGI;
            tFollowerFIManager.mFI = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator = &tGI;
            tFollowerFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager = &tFollowerFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tLeaderFIManager, mtk::Leader_Follower::LEADER );
            tIWG->set_field_interpolator_manager( &tFollowerFIManager, mtk::Leader_Follower::FOLLOWER );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mFollowerFIManager->set_space_time( tParamPoint );

                // check evaluation of the residual for IWG
                //------------------------------------------------------------------------------
                // reset residual
                tIWG->mSet->mResidual( 0 ).fill( 0.0 );

                // compute residual
                tIWG->compute_residual( 1.0 );

                // check evaluation of the jacobian by FD
                //------------------------------------------------------------------------------
                // reset jacobian
                tIWG->mSet->mJacobian.fill( 0.0 );

                // init the jacobian for IWG and FD evaluation
                Matrix< DDRMat > tJacobian;
                Matrix< DDRMat > tJacobianFD;

                // check jacobian by FD
                bool tCheckJacobian = tIWG->check_jacobian(
                        tPerturbation,
                        tEpsilon,
                        1.0,
                        tJacobian,
                        tJacobianFD,
                        true );

                // print for debug
                if( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << '\n';
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
}/* END_TEST_CASE */

