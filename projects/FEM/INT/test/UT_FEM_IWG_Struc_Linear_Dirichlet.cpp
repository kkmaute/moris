/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Struc_Linear_Dirichlet.cpp
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
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
//LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Struc_Dirichlet_Symmetric_Nitsche",
        "[moris],[fem],[IWG_Struc_Dirichlet_Symmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-5;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = tDispDofTypes;

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderSelect = std::make_shared< fem::Property >();
    tPropLeaderSelect->set_val_function( tMValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderStrucLinIso->set_dof_type_list( { tDispDofTypes } );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 100.0 }} } );
    tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPDirichletNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderSelect, "Select" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
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

        // set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set selection matrix parameters
        tPropLeaderSelect->set_parameters( { {{ static_cast< real >( iSpaceDim ) }} } );

        // Dirichlet BC
        Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.5 );
        tPropLeaderDirichlet->set_parameters( { tImposedDisp } );

        // set space dimension to CM, SP
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );

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
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp-1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp, tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Dirichlet_Unsymmetric_Nitsche",
        "[moris],[fem],[IWG_Struc_Dirichlet_Unsymmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-5;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = tDispDofTypes;

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderSelect = std::make_shared< fem::Property >();
    tPropLeaderSelect->set_val_function( tMValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderStrucLinIso->set_dof_type_list( { tDispDofTypes } );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 100.0 }} } );
    tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPDirichletNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderSelect, "Select" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
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

        // set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set selection matrix parameters
        tPropLeaderSelect->set_parameters( { {{ static_cast< real >( iSpaceDim ) }} } );

        // Dirichlet BC
        Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.5 );
        tPropLeaderDirichlet->set_parameters( { tImposedDisp } );

        // set space dimension to CM, SP
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );

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
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp-1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp, tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Dirichlet_Mixed_Displacement_Symmetric_Nitsche",
        "[IWG_Struc_Dirichlet_Mixed_Displacement_Symmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tPDofTypes = { MSI::Dof_Type::P };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tDispDofTypes( 0 ), tPDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderSelect = std::make_shared< fem::Property >();
    tPropLeaderSelect->set_val_function( tMValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderStrucLinIso->set_dof_type_list( tDofTypes, { "Displacement", "Pressure" } );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPDirichletNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDofTypes );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderSelect, "Select" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set selection matrix parameters
        tPropLeaderSelect->set_parameters( { {{ static_cast< real >( iSpaceDim ) }} } );

        // Dirichlet BC
        Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.5 );
        tPropLeaderDirichlet->set_parameters( { tImposedDisp } );

        // set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

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

        // set space dimension to CM, SP
        tCMLeaderStrucLinIso->set_model_type(fem::Model_Type::PLANE_STRESS);
        tCMLeaderStrucLinIso->set_model_type(fem::Model_Type::DEVIATORIC);
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );

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
            int tNumDofP     = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDofHatP;
            fill_phat_Elast( tLeaderDofHatP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDofHatP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, tNumDofDisp + tNumDofP - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofDisp - 1 },
                    { tNumDofDisp, tNumDofDisp + tNumDofP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp + tNumDofP, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp + tNumDofP, tNumDofDisp + tNumDofP, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Dirichlet_Mixed_Displacement_Unsymmetric_Nitsche",
        "[IWG_Struc_Dirichlet_Mixed_Displacement_Unsymmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tPDofTypes = { MSI::Dof_Type::P };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tDispDofTypes( 0 ), tPDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderSelect = std::make_shared< fem::Property >();
    tPropLeaderSelect->set_val_function( tMValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderStrucLinIso->set_dof_type_list( tDofTypes, { "Displacement", "Pressure" } );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPDirichletNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDofTypes );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderSelect, "Select" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set selection matrix parameters
        tPropLeaderSelect->set_parameters( { {{ static_cast< real >( iSpaceDim ) }} } );

        // Dirichlet BC
        Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.5 );
        tPropLeaderDirichlet->set_parameters( { tImposedDisp } );

        // set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

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

        // set space dimension to CM, SP
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::DEVIATORIC );
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );

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
            int tNumDofP     = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDofHatP;
            fill_phat_Elast( tLeaderDofHatP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDofHatP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, tNumDofDisp + tNumDofP - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofDisp - 1 },
                    { tNumDofDisp, tNumDofDisp + tNumDofP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp + tNumDofP, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp + tNumDofP, tNumDofDisp + tNumDofP, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Thermo_Elastic_Dirichlet_Symmetric_Nitsche",
        "[moris],[fem],[IWG_Struc_Thermo_Elastic_Dirichlet_Symmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tDispDofTypes( 0 ), tTempDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderSelect = std::make_shared< fem::Property >();
    tPropLeaderSelect->set_val_function( tMValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderCTE = std::make_shared< fem::Property >();
    tPropLeaderCTE->set_parameters( {{{ 1.0 }}} );
    tPropLeaderCTE->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderTRef = std::make_shared< fem::Property >();
    tPropLeaderTRef->set_parameters( {{{ 1.25 }}} );
    tPropLeaderTRef->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderStrucLinIso->set_dof_type_list( tDofTypes, { "Displacement", "Temperature" } );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderCTE, "CTE" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderTRef, "ReferenceTemperature" );
    tCMLeaderStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPDirichletNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDofTypes );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderSelect, "Select" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set selection matrix parameters
        tPropLeaderSelect->set_parameters( { {{ static_cast< real >( iSpaceDim ) }} } );

        // Dirichlet BC
        Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.5 );
        tPropLeaderDirichlet->set_parameters( { tImposedDisp } );

        // set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

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

        // set space dimension to CM, SP
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::DEVIATORIC );
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );

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
            int tNumDofTemp  = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDofHatTemp;
            fill_phat_Elast( tLeaderDofHatTemp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDofHatTemp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, tNumDofDisp + tNumDofTemp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofDisp - 1 },
                    { tNumDofDisp, tNumDofDisp + tNumDofTemp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp + tNumDofTemp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp + tNumDofTemp, tNumDofDisp + tNumDofTemp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Thermo_Elastic_Dirichlet_Unsymmetric_Nitsche",
        "[moris],[fem],[IWG_Struc_Thermo_Elastic_Dirichlet_Unsymmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tDispDofTypes( 0 ), tTempDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderSelect = std::make_shared< fem::Property >();
    tPropLeaderSelect->set_val_function( tMValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderCTE = std::make_shared< fem::Property >();
    tPropLeaderCTE->set_parameters( {{{ 1.0 }}} );
    tPropLeaderCTE->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderTRef = std::make_shared< fem::Property >();
    tPropLeaderTRef->set_parameters( {{{ 1.25 }}} );
    tPropLeaderTRef->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderStrucLinIso->set_dof_type_list( tDofTypes, { "Displacement", "Temperature" } );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderCTE, "CTE" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderTRef, "ReferenceTemperature" );
    tCMLeaderStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPDirichletNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDofTypes );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderSelect, "Select" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set selection matrix parameters
        tPropLeaderSelect->set_parameters( { {{ static_cast< real >( iSpaceDim ) }} } );

        // Dirichlet BC
        Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.5 );
        tPropLeaderDirichlet->set_parameters( { tImposedDisp } );

        // set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

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

        // set space dimension to CM, SP
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::DEVIATORIC );
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );

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
            int tNumDofTemp  = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDofHatTemp;
            fill_phat_Elast( tLeaderDofHatTemp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDofHatTemp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, tNumDofDisp + tNumDofTemp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofDisp - 1 },
                    { tNumDofDisp, tNumDofDisp + tNumDofTemp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp + tNumDofTemp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp + tNumDofTemp, tNumDofDisp + tNumDofTemp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Linear_Elasticity_Damage_Dirichlet_Symmetric_Nitsche",
        "[moris],[fem],[IWG_Struc_Linear_Elasticity_Damage_Dirichlet_Symmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 2E-5;

    // define aperturbation relative size
    real tPerturbation = 1E-5;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tDisplDofTypes = { MSI::Dof_Type::UX };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes           = { tDisplDofTypes, tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // Damage setup - 0,2,0
    // Local equivalent strain - 0 - energy release rate
    // Damage law - 2 - smooth exponential
    // Smoothing law - 2 - smoothing
    moris::Cell< Matrix< DDRMat > > tDamageParameters = { //
            { { 0.0 } },                                      //
            { { 2.0, 1.0e-3, 10.0 } },                        //
            { { 0.0 } }
    };

    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( {{{ 10.0 }}} );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( {{{ 0.3 }}} );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderSelect = std::make_shared< fem::Property >();
    tPropLeaderSelect->set_val_function( tMValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeader->set_local_properties();
    tCMLeader->set_parameters( tDamageParameters );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 100.0 }} } );
    tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPDirichletNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( { tDisplDofTypes } );
    tIWG->set_dof_type_list( tDofTypes );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMLeader, "ElastLinIso" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderSelect, "Select" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
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

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
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

        // set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set selection matrix parameters
        tPropLeaderSelect->set_parameters( { {{ static_cast< real >( iSpaceDim ) }} } );

        // Dirichlet BC
        Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.5 );
        tPropLeaderDirichlet->set_parameters( { tImposedDisp } );

        // set space dimension to CM, SP
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeader->set_space_dim( iSpaceDim );

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

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 5.0e-3*tLeaderDOFHatNlEqStrain;
            Matrix< DDRMat > tLeaderDOFHatHistory;
            fill_phat_Elast( tLeaderDOFHatHistory, iSpaceDim, iInterpOrder );
            tLeaderDOFHatHistory = 5.0e-3*tLeaderDOFHatHistory;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDisplDofTypes );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatHistory );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDispl - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDispl, tNumDofDispl + tNumDofNlEqStrain - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofDispl + tNumDofNlEqStrain, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofDispl - 1 },
                    { tNumDofDispl, tNumDofDispl + tNumDofNlEqStrain - 1 },
                    { tNumDofDispl + tNumDofNlEqStrain, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/* END_TEST_CASE */

