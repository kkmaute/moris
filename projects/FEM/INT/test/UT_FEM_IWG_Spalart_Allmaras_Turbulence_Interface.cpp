/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.cpp
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
// FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"
// LINALG
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

// This UT tests the velocity interface IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order

TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Interface_Symmetric_Nitsche",
        "[IWG_Spalart_Allmaras_Turbulence_Interface_Symmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

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
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };

    // dof type list
    Vector< Vector< MSI::Dof_Type > > tVisDofTypes = { { MSI::Dof_Type::VISCOSITY } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes    = { tVisDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // create the leader properties
    std::shared_ptr< fem::Property > tPropLeaderViscosity = std::make_shared< fem::Property >();
    tPropLeaderViscosity->set_parameters( { { { 2.0 } } } );
    tPropLeaderViscosity->set_val_function( tConstValFunc );
    // tPropLeaderViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropLeaderViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropLeaderViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderWallDistance = std::make_shared< fem::Property >();
    tPropLeaderWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropLeaderWallDistance->set_val_function( tConstValFunc );

    // create the follower properties
    std::shared_ptr< fem::Property > tPropFollowerViscosity = std::make_shared< fem::Property >();
    tPropFollowerViscosity->set_parameters( { { { 1.0 } } } );
    tPropFollowerViscosity->set_val_function( tConstValFunc );
    // tPropFollowerViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropFollowerViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropFollowerViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerWallDistance = std::make_shared< fem::Property >();
    tPropFollowerWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropFollowerWallDistance->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMLeaderSATurbulence->set_dof_type_list( tDofTypes );
    tCMLeaderSATurbulence->set_property( tPropLeaderWallDistance, "WallDistance" );
    tCMLeaderSATurbulence->set_property( tPropLeaderViscosity, "KinViscosity" );
    tCMLeaderSATurbulence->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMFollowerSATurbulence->set_dof_type_list( tDofTypes );
    tCMFollowerSATurbulence->set_property( tPropFollowerWallDistance, "WallDistance" );
    tCMFollowerSATurbulence->set_property( tPropFollowerViscosity, "KinViscosity" );
    tCMFollowerSATurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::SPALART_ALLMARAS_NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { { { 1.0 } } } );
    tSPNitsche->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_constitutive_model( tCMFollowerSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_INTERFACE_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tVisDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 } };
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // space and time geometry interpolators
        //------------------------------------------------------------------------------
        // create a space geometry interpolation rule
        mtk::Interpolation_Rule tGIRule( tGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tCMLeaderSATurbulence->set_space_dim( iSpaceDim );
        tCMFollowerSATurbulence->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
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
            int tNumDofVis = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVis );

            // fill random coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatVis;
            fill_phat( tFollowerDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tFollowerFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVis, 2 * tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofVis - 1 },
                { tNumDofVis, 2 * tNumDofVis - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofVis, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofVis, 2 * tNumDofVis, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >        tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tLeaderFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager                         tFollowerFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderFIManager.mFI                     = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator = &tGI;
            tLeaderFIManager.mIGGeometryInterpolator = &tGI;
            tFollowerFIManager.mFI                      = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator  = &tGI;
            tFollowerFIManager.mIGGeometryInterpolator  = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager  = &tFollowerFIManager;

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
                if ( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Interface_Symmetric_Nitsche_Negative",
        "[IWG_Spalart_Allmaras_Turbulence_Interface_Symmetric_Nitsche_Negative]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

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
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };

    // dof type list
    Vector< Vector< MSI::Dof_Type > > tVisDofTypes = { { MSI::Dof_Type::VISCOSITY } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes    = { tVisDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // create the leader properties
    std::shared_ptr< fem::Property > tPropLeaderViscosity = std::make_shared< fem::Property >();
    tPropLeaderViscosity->set_parameters( { { { 2.0 } } } );
    tPropLeaderViscosity->set_val_function( tConstValFunc );
    // tPropLeaderViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropLeaderViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropLeaderViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderWallDistance = std::make_shared< fem::Property >();
    tPropLeaderWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropLeaderWallDistance->set_val_function( tConstValFunc );

    // create the follower properties
    std::shared_ptr< fem::Property > tPropFollowerViscosity = std::make_shared< fem::Property >();
    tPropFollowerViscosity->set_parameters( { { { 1.0 } } } );
    tPropFollowerViscosity->set_val_function( tConstValFunc );
    // tPropFollowerViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropFollowerViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropFollowerViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerWallDistance = std::make_shared< fem::Property >();
    tPropFollowerWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropFollowerWallDistance->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMLeaderSATurbulence->set_dof_type_list( tDofTypes );
    tCMLeaderSATurbulence->set_property( tPropLeaderWallDistance, "WallDistance" );
    tCMLeaderSATurbulence->set_property( tPropLeaderViscosity, "KinViscosity" );
    tCMLeaderSATurbulence->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMFollowerSATurbulence->set_dof_type_list( tDofTypes );
    tCMFollowerSATurbulence->set_property( tPropFollowerWallDistance, "WallDistance" );
    tCMFollowerSATurbulence->set_property( tPropFollowerViscosity, "KinViscosity" );
    tCMFollowerSATurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::SPALART_ALLMARAS_NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { { { 1.0 } } } );
    tSPNitsche->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_constitutive_model( tCMFollowerSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_INTERFACE_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tVisDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 } };
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // space and time geometry interpolators
        //------------------------------------------------------------------------------
        // create a space geometry interpolation rule
        mtk::Interpolation_Rule tGIRule( tGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tCMLeaderSATurbulence->set_space_dim( iSpaceDim );
        tCMFollowerSATurbulence->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
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
            int tNumDofVis = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );
            tLeaderDOFHatVis = -1.0 * tLeaderDOFHatVis;

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVis );

            // fill random coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatVis;
            fill_phat( tFollowerDOFHatVis, iSpaceDim, iInterpOrder );
            tFollowerDOFHatVis = -1.0 * tFollowerDOFHatVis;

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tFollowerFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVis, 2 * tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofVis - 1 },
                { tNumDofVis, 2 * tNumDofVis - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofVis, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofVis, 2 * tNumDofVis, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >        tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tLeaderFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager                         tFollowerFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderFIManager.mFI                     = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator = &tGI;
            tLeaderFIManager.mIGGeometryInterpolator = &tGI;
            tFollowerFIManager.mFI                      = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator  = &tGI;
            tFollowerFIManager.mIGGeometryInterpolator  = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager  = &tFollowerFIManager;

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
                if ( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Interface_Unsymmetric_Nitsche",
        "[IWG_Spalart_Allmaras_Turbulence_Interface_Unsymmetric_Nitsche]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

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
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };

    // dof type list
    Vector< Vector< MSI::Dof_Type > > tVisDofTypes = { { MSI::Dof_Type::VISCOSITY } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes    = { tVisDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // create the leader properties
    std::shared_ptr< fem::Property > tPropLeaderViscosity = std::make_shared< fem::Property >();
    tPropLeaderViscosity->set_parameters( { { { 2.0 } } } );
    tPropLeaderViscosity->set_val_function( tConstValFunc );
    // tPropLeaderViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropLeaderViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropLeaderViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderWallDistance = std::make_shared< fem::Property >();
    tPropLeaderWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropLeaderWallDistance->set_val_function( tConstValFunc );

    // create the follower properties
    std::shared_ptr< fem::Property > tPropFollowerViscosity = std::make_shared< fem::Property >();
    tPropFollowerViscosity->set_parameters( { { { 1.0 } } } );
    tPropFollowerViscosity->set_val_function( tConstValFunc );
    // tPropFollowerViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropFollowerViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropFollowerViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerWallDistance = std::make_shared< fem::Property >();
    tPropFollowerWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropFollowerWallDistance->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMLeaderSATurbulence->set_dof_type_list( tDofTypes );
    tCMLeaderSATurbulence->set_property( tPropLeaderWallDistance, "WallDistance" );
    tCMLeaderSATurbulence->set_property( tPropLeaderViscosity, "KinViscosity" );
    tCMLeaderSATurbulence->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMFollowerSATurbulence->set_dof_type_list( tDofTypes );
    tCMFollowerSATurbulence->set_property( tPropFollowerWallDistance, "WallDistance" );
    tCMFollowerSATurbulence->set_property( tPropFollowerViscosity, "KinViscosity" );
    tCMFollowerSATurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::SPALART_ALLMARAS_NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { { { 1.0 } } } );
    tSPNitsche->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_constitutive_model( tCMFollowerSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_INTERFACE_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tVisDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 } };
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // space and time geometry interpolators
        //------------------------------------------------------------------------------
        // create a space geometry interpolation rule
        mtk::Interpolation_Rule tGIRule( tGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tSPNitsche->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
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
            int tNumDofVis = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVis );

            // fill random coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatVis;
            fill_phat( tFollowerDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tFollowerFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVis, 2 * tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofVis - 1 },
                { tNumDofVis, 2 * tNumDofVis - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofVis, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofVis, 2 * tNumDofVis, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >        tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tLeaderFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager                         tFollowerFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderFIManager.mFI                     = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator = &tGI;
            tLeaderFIManager.mIGGeometryInterpolator = &tGI;
            tFollowerFIManager.mFI                      = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator  = &tGI;
            tFollowerFIManager.mIGGeometryInterpolator  = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager  = &tFollowerFIManager;

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
                if ( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Interface_Unsymmetric_Nitsche_Negative",
        "[IWG_Spalart_Allmaras_Turbulence_Interface_Unsymmetric_Nitsche_Negative]" )
{
    // define an epsilon environment
    real tEpsilon = 5E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

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
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };

    // dof type list
    Vector< Vector< MSI::Dof_Type > > tVisDofTypes = { { MSI::Dof_Type::VISCOSITY } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes    = { tVisDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // create the leader properties
    std::shared_ptr< fem::Property > tPropLeaderViscosity = std::make_shared< fem::Property >();
    tPropLeaderViscosity->set_parameters( { { { 2.0 } } } );
    tPropLeaderViscosity->set_val_function( tConstValFunc );
    // tPropLeaderViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropLeaderViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropLeaderViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderWallDistance = std::make_shared< fem::Property >();
    tPropLeaderWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropLeaderWallDistance->set_val_function( tConstValFunc );

    // create the follower properties
    std::shared_ptr< fem::Property > tPropFollowerViscosity = std::make_shared< fem::Property >();
    tPropFollowerViscosity->set_parameters( { { { 1.0 } } } );
    tPropFollowerViscosity->set_val_function( tConstValFunc );
    // tPropFollowerViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropFollowerViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropFollowerViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerWallDistance = std::make_shared< fem::Property >();
    tPropFollowerWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropFollowerWallDistance->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMLeaderSATurbulence->set_dof_type_list( tDofTypes );
    tCMLeaderSATurbulence->set_property( tPropLeaderWallDistance, "WallDistance" );
    tCMLeaderSATurbulence->set_property( tPropLeaderViscosity, "KinViscosity" );
    tCMLeaderSATurbulence->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMFollowerSATurbulence->set_dof_type_list( tDofTypes );
    tCMFollowerSATurbulence->set_property( tPropFollowerWallDistance, "WallDistance" );
    tCMFollowerSATurbulence->set_property( tPropFollowerViscosity, "KinViscosity" );
    tCMFollowerSATurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::SPALART_ALLMARAS_NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { { { 1.0 } } } );
    tSPNitsche->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_constitutive_model( tCMFollowerSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_INTERFACE_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tVisDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerSATurbulence, "SpalartAllmarasTurbulence", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 } };
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // space and time geometry interpolators
        //------------------------------------------------------------------------------
        // create a space geometry interpolation rule
        mtk::Interpolation_Rule tGIRule( tGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tSPNitsche->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
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
            int tNumDofVis = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );
            tLeaderDOFHatVis = -1.0 * tLeaderDOFHatVis;

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVis );

            // fill random coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatVis;
            fill_phat( tFollowerDOFHatVis, iSpaceDim, iInterpOrder );
            tFollowerDOFHatVis = -1.0 * tFollowerDOFHatVis;

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tFollowerFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVis, 2 * tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofVis - 1 },
                { tNumDofVis, 2 * tNumDofVis - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofVis, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofVis, 2 * tNumDofVis, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >        tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tLeaderFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager                         tFollowerFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderFIManager.mFI                     = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator = &tGI;
            tLeaderFIManager.mIGGeometryInterpolator = &tGI;
            tFollowerFIManager.mFI                      = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator  = &tGI;
            tFollowerFIManager.mIGGeometryInterpolator  = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager  = &tFollowerFIManager;

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
                if ( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
} /*END_TEST_CASE*/

