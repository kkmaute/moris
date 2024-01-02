/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Incompressible_NS_Interface_Nitsche.cpp
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
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"
//LINALG
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

// This UT tests the velocity interface IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order

TEST_CASE( "IWG_Incompressible_NS_Velocity_Interface_Symmetric_Nitsche",
        "[IWG_Incompressible_NS_Velocity_Interface_Symmetric_Nitsche]" )
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
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the leader properties
    std::shared_ptr< fem::Property > tPropLeaderViscosity = std::make_shared< fem::Property >();
    tPropLeaderViscosity->set_parameters( { {{ 1.0 }} } );
    tPropLeaderViscosity->set_val_function( tConstValFunc );
    //tPropLeaderViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropLeaderViscosity->set_val_function( tVXFIValFunc );
    //tPropLeaderViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderKinViscosity = std::make_shared< fem::Property >();
    tPropLeaderKinViscosity->set_val_function( tConstValFunc );
    tPropLeaderKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderDensity = std::make_shared< fem::Property >();
    tPropLeaderDensity->set_parameters( { {{ 1.0 }} } );
    tPropLeaderDensity->set_val_function( tConstValFunc );
    //tPropLeaderViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropLeaderViscosity->set_val_function( tVXFIValFunc );
    //tPropLeaderViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // create the follower properties
    std::shared_ptr< fem::Property > tPropFollowerViscosity = std::make_shared< fem::Property >();
    tPropFollowerViscosity->set_parameters( { {{ 1.0 }} } );
    tPropFollowerViscosity->set_val_function( tConstValFunc );
    //tPropFollowerViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropFollowerViscosity->set_val_function( tVXFIValFunc );
    //tPropFollowerViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerKinViscosity = std::make_shared< fem::Property >();
    tPropFollowerKinViscosity->set_val_function( tConstValFunc );
    tPropFollowerKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerDensity = std::make_shared< fem::Property >();
    tPropFollowerDensity->set_parameters( { {{ 1.0 }} } );
    tPropFollowerDensity->set_val_function( tConstValFunc );
    //tPropFollowerDensity->set_dof_type_list( { tVelDofTypes } );
    //tPropFollowerDensity->set_val_function( tVXFIValFunc );
    //tPropFollowerDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMLeaderTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMLeaderTurbulence->set_property( tPropLeaderViscosity, "Viscosity" );
    tCMLeaderTurbulence->set_property( tPropLeaderKinViscosity, "KinViscosity" );
    tCMLeaderTurbulence->set_property( tPropLeaderDensity, "Density" );
    tCMLeaderTurbulence->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMFollowerTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMFollowerTurbulence->set_property( tPropFollowerViscosity, "Viscosity" );
    tCMFollowerTurbulence->set_property( tPropFollowerKinViscosity, "KinViscosity" );
    tCMFollowerTurbulence->set_property( tPropFollowerDensity, "Density" );
    tCMFollowerTurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { {{ 1.0 }} } );
    tSPNitsche->set_property( tPropLeaderViscosity, "Material", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_property( tPropFollowerViscosity, "Material", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_INTERFACE_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tVelDofTypes );
    tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes }, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes }, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderTurbulence, "IncompressibleFluid", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerTurbulence, "IncompressibleFluid", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) )       = 0;
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )        = 1;
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                         { 1.0, 0.0 },
                         { 1.0, 1.0 },
                         { 0.0, 1.0 }};

               // set velocity dof types
               tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };

               // set viscosity property parameters
               tPropLeaderKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
               tPropFollowerKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
               break;
            }
            case 3 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0, 0.0 },
                         { 1.0, 0.0, 0.0 },
                         { 1.0, 1.0, 0.0 },
                         { 0.0, 1.0, 0.0 },
                         { 0.0, 0.0, 1.0 },
                         { 1.0, 0.0, 1.0 },
                         { 1.0, 1.0, 1.0 },
                         { 0.0, 1.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };

                // set viscosity property parameters
                tPropLeaderKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
                tPropFollowerKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
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
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tCMLeaderTurbulence->set_space_dim( iSpaceDim );
        tCMFollowerTurbulence->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
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
            int tNumDofVel = tNumCoeff * iSpaceDim;
            int tNumDofP   = tNumCoeff;
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator viscosity
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatVis );

            // fill coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatVel;
            fill_uhat( tFollowerDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tFollowerDOFHatP;
            fill_phat( tFollowerDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tFollowerDOFHatVis;
            fill_phat( tFollowerDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVel );

            // create the field interpolator pressure
            tFollowerFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tFollowerFIs( 1 )->set_coeff( tFollowerDOFHatP );

            // create the field interpolator pressure
            tFollowerFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tFollowerFIs( 2 )->set_coeff( tFollowerDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofP + tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 },
                    { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };
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
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
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

                // check that jacobian is symmetric for symmetric version using upper diagonal block
                 uint tNumRowBlock = tJacobian.n_rows()/iSpaceDim;
                 Matrix< DDRMat > tBlock = tJacobian({ 0, tNumRowBlock-1 },{ 0, tNumRowBlock-1 });

                 real tRelError = norm( tBlock - trans( tBlock ) ) / norm ( tBlock );
                 REQUIRE( tRelError < 1e-12 );

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
            tFollowerFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "IWG_Incompressible_NS_Velocity_Interface_Unsymmetric_Nitsche",
        "[IWG_Incompressible_NS_Velocity_Interface_Unsymmetric_Nitsche]" )
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
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the leader properties
    std::shared_ptr< fem::Property > tPropLeaderViscosity = std::make_shared< fem::Property >();
    tPropLeaderViscosity->set_parameters( { {{ 1.0 }} } );
    tPropLeaderViscosity->set_val_function( tConstValFunc );
    //tPropLeaderViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropLeaderViscosity->set_val_function( tVXFIValFunc );
    //tPropLeaderViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderKinViscosity = std::make_shared< fem::Property >();
    tPropLeaderKinViscosity->set_val_function( tConstValFunc );
    tPropLeaderKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderDensity = std::make_shared< fem::Property >();
    tPropLeaderDensity->set_parameters( { {{ 1.0 }} } );
    tPropLeaderDensity->set_val_function( tConstValFunc );
    //tPropLeaderViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropLeaderViscosity->set_val_function( tVXFIValFunc );
    //tPropLeaderViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // create the follower properties
    std::shared_ptr< fem::Property > tPropFollowerViscosity = std::make_shared< fem::Property >();
    tPropFollowerViscosity->set_parameters( { {{ 1.0 }} } );
    tPropFollowerViscosity->set_val_function( tConstValFunc );
    //tPropFollowerViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropFollowerViscosity->set_val_function( tVXFIValFunc );
    //tPropFollowerViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerKinViscosity = std::make_shared< fem::Property >();
    tPropFollowerKinViscosity->set_val_function( tConstValFunc );
    tPropFollowerKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerDensity = std::make_shared< fem::Property >();
    tPropFollowerDensity->set_parameters( { {{ 1.0 }} } );
    tPropFollowerDensity->set_val_function( tConstValFunc );
    //tPropFollowerDensity->set_dof_type_list( { tVelDofTypes } );
    //tPropFollowerDensity->set_val_function( tVXFIValFunc );
    //tPropFollowerDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMLeaderTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMLeaderTurbulence->set_property( tPropLeaderViscosity, "Viscosity" );
    tCMLeaderTurbulence->set_property( tPropLeaderKinViscosity, "KinViscosity" );
    tCMLeaderTurbulence->set_property( tPropLeaderDensity, "Density" );
    tCMLeaderTurbulence->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMFollowerTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMFollowerTurbulence->set_property( tPropFollowerViscosity, "Viscosity" );
    tCMFollowerTurbulence->set_property( tPropFollowerKinViscosity, "KinViscosity" );
    tCMFollowerTurbulence->set_property( tPropFollowerDensity, "Density" );
    tCMFollowerTurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { {{ 1.0 }} } );
    tSPNitsche->set_property( tPropLeaderViscosity, "Material", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_property( tPropFollowerViscosity, "Material", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_INTERFACE_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tVelDofTypes );
    tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes }, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes }, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderTurbulence, "IncompressibleFluid", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerTurbulence, "IncompressibleFluid", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) )       = 0;
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )        = 1;
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                         { 1.0, 0.0 },
                         { 1.0, 1.0 },
                         { 0.0, 1.0 }};

               // set velocity dof types
               tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };

               // set viscosity property parameters
               tPropLeaderKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
               tPropFollowerKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
               break;
            }
            case 3 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0, 0.0 },
                         { 1.0, 0.0, 0.0 },
                         { 1.0, 1.0, 0.0 },
                         { 0.0, 1.0, 0.0 },
                         { 0.0, 0.0, 1.0 },
                         { 1.0, 0.0, 1.0 },
                         { 1.0, 1.0, 1.0 },
                         { 0.0, 1.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };

                // set viscosity property parameters
                tPropLeaderKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
                tPropFollowerKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
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
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tCMLeaderTurbulence->set_space_dim( iSpaceDim );
        tCMFollowerTurbulence->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
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
            int tNumDofVel = tNumCoeff * iSpaceDim;
            int tNumDofP   = tNumCoeff;
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator viscosity
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatVis );

            // fill coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatVel;
            fill_uhat( tFollowerDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tFollowerDOFHatP;
            fill_phat( tFollowerDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tFollowerDOFHatVis;
            fill_phat( tFollowerDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVel );

            // create the field interpolator pressure
            tFollowerFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tFollowerFIs( 1 )->set_coeff( tFollowerDOFHatP );

            // create the field interpolator pressure
            tFollowerFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tFollowerFIs( 2 )->set_coeff( tFollowerDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofP + tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 },
                    { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };
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
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
}/*END_TEST_CASE*/
// This UT tests the pressure interface IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order

TEST_CASE( "IWG_Incompressible_NS_Pressure_Interface_Symmetric_Nitsche",
        "[IWG_Incompressible_NS_Pressure_Interface_Symmetric_Nitsche]" )
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
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< MSI::Dof_Type > tVisDofTypes  = { MSI::Dof_Type::VISCOSITY };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the leader properties
    std::shared_ptr< fem::Property > tPropLeaderViscosity = std::make_shared< fem::Property >();
    tPropLeaderViscosity->set_parameters( { {{ 1.0 }} } );
    tPropLeaderViscosity->set_val_function( tConstValFunc );
//    tPropLeaderViscosity->set_dof_type_list( { tVelDofTypes } );
//    tPropLeaderViscosity->set_val_function( tVXFIValFunc );
//    tPropLeaderViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderKinViscosity = std::make_shared< fem::Property >();
    tPropLeaderKinViscosity->set_val_function( tConstValFunc );
    tPropLeaderKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderDensity = std::make_shared< fem::Property >();
    tPropLeaderDensity->set_parameters( { {{ 1.0 }} } );
    tPropLeaderDensity->set_val_function( tConstValFunc );

//    tPropLeaderDensity->set_dof_type_list( { tVelDofTypes } );
//    tPropLeaderDensity->set_val_function( tVXFIValFunc );
//    tPropLeaderDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // create the follower properties
    std::shared_ptr< fem::Property > tPropFollowerViscosity = std::make_shared< fem::Property >();
    tPropFollowerViscosity->set_parameters( { {{ 1.0 }} } );
    tPropFollowerViscosity->set_val_function( tConstValFunc );
//    tPropFollowerViscosity->set_dof_type_list( { tVelDofTypes } );
//    tPropFollowerViscosity->set_val_function( tVXFIValFunc );
//    tPropFollowerViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerKinViscosity = std::make_shared< fem::Property >();
    tPropFollowerKinViscosity->set_val_function( tConstValFunc );
    tPropFollowerKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerDensity = std::make_shared< fem::Property >();
    tPropFollowerDensity->set_parameters( { {{ 1.0 }} } );
    tPropFollowerDensity->set_val_function( tConstValFunc );
//    tPropFollowerDensity->set_dof_type_list( { tVelDofTypes } );
//    tPropFollowerDensity->set_val_function( tVXFIValFunc );
//    tPropFollowerDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMLeaderFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMLeaderFluid->set_property( tPropLeaderViscosity, "Viscosity" );
    tCMLeaderFluid->set_property( tPropLeaderKinViscosity, "KinViscosity" );
    tCMLeaderFluid->set_property( tPropLeaderDensity, "Density" );
    tCMLeaderFluid->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMFollowerFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMFollowerFluid->set_property( tPropFollowerViscosity, "Viscosity" );
    tCMFollowerFluid->set_property( tPropFollowerKinViscosity, "KinViscosity" );
    tCMFollowerFluid->set_property( tPropFollowerDensity, "Density" );
    tCMFollowerFluid->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { {{ 1.0 }} } );
    tSPNitsche->set_property( tPropLeaderViscosity, "Material", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_property( tPropFollowerViscosity, "Material", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tPDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderFluid, "IncompressibleFluid", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerFluid, "IncompressibleFluid", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                         { 1.0, 0.0 },
                         { 1.0, 1.0 },
                         { 0.0, 1.0 }};

               // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };

                // set viscosity property parameters
                tPropLeaderKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
                tPropFollowerKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
               break;
            }
            case 3 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0, 0.0 },
                         { 1.0, 0.0, 0.0 },
                         { 1.0, 1.0, 0.0 },
                         { 0.0, 1.0, 0.0 },
                         { 0.0, 0.0, 1.0 },
                         { 1.0, 0.0, 1.0 },
                         { 1.0, 1.0, 1.0 },
                         { 0.0, 1.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };

                // set viscosity property parameters
                tPropLeaderKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
                tPropFollowerKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
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
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tCMLeaderFluid->set_space_dim( iSpaceDim );
        tCMFollowerFluid->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
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
            int tNumDofVel  = tNumCoeff * iSpaceDim;
            int tNumDofP    = tNumCoeff;
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator viscosity
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatVis );

            // fill coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatVel;
            fill_uhat( tFollowerDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tFollowerDOFHatP;
            fill_phat( tFollowerDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tFollowerDOFHatVis;
            fill_phat( tFollowerDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVel );

            // create the field interpolator pressure
            tFollowerFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tFollowerFIs( 1 )->set_coeff( tFollowerDOFHatP );

            // create the field interpolator pressure
            tFollowerFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tFollowerFIs( 2 )->set_coeff( tFollowerDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofP + tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 },
                    { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };
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
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "IWG_Incompressible_NS_Pressure_Interface_Unsymmetric_Nitsche",
        "[IWG_Incompressible_NS_Pressure_Interface_Unsymmetric_Nitsche]" )
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
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the leader properties
    std::shared_ptr< fem::Property > tPropLeaderViscosity = std::make_shared< fem::Property >();
    tPropLeaderViscosity->set_parameters( { {{ 1.0 }} } );
    tPropLeaderViscosity->set_val_function( tConstValFunc );
//    tPropLeaderViscosity->set_dof_type_list( { tVelDofTypes } );
//    tPropLeaderViscosity->set_val_function( tVXFIValFunc );
//    tPropLeaderViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderKinViscosity = std::make_shared< fem::Property >();
    tPropLeaderKinViscosity->set_val_function( tConstValFunc );
    tPropLeaderKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropLeaderDensity = std::make_shared< fem::Property >();
    tPropLeaderDensity->set_parameters( { {{ 1.0 }} } );
    tPropLeaderDensity->set_val_function( tConstValFunc );

//    tPropLeaderDensity->set_dof_type_list( { tVelDofTypes } );
//    tPropLeaderDensity->set_val_function( tVXFIValFunc );
//    tPropLeaderDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // create the follower properties
    std::shared_ptr< fem::Property > tPropFollowerViscosity = std::make_shared< fem::Property >();
    tPropFollowerViscosity->set_parameters( { {{ 1.0 }} } );
    tPropFollowerViscosity->set_val_function( tConstValFunc );
//    tPropFollowerViscosity->set_dof_type_list( { tVelDofTypes } );
//    tPropFollowerViscosity->set_val_function( tVXFIValFunc );
//    tPropFollowerViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerKinViscosity = std::make_shared< fem::Property >();
    tPropFollowerKinViscosity->set_val_function( tConstValFunc );
    tPropFollowerKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropFollowerDensity = std::make_shared< fem::Property >();
    tPropFollowerDensity->set_parameters( { {{ 1.0 }} } );
    tPropFollowerDensity->set_val_function( tConstValFunc );
//    tPropFollowerDensity->set_dof_type_list( { tVelDofTypes } );
//    tPropFollowerDensity->set_val_function( tVXFIValFunc );
//    tPropFollowerDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMLeaderFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMLeaderFluid->set_property( tPropLeaderViscosity, "Viscosity" );
    tCMLeaderFluid->set_property( tPropLeaderKinViscosity, "KinViscosity" );
    tCMLeaderFluid->set_property( tPropLeaderDensity, "Density" );
    tCMLeaderFluid->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMFollowerFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMFollowerFluid->set_property( tPropFollowerViscosity, "Viscosity" );
    tCMFollowerFluid->set_property( tPropFollowerKinViscosity, "KinViscosity" );
    tCMFollowerFluid->set_property( tPropFollowerDensity, "Density" );
    tCMFollowerFluid->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { {{ 1.0 }} } );
    tSPNitsche->set_property( tPropLeaderViscosity, "Material", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_property( tPropFollowerViscosity, "Material", mtk::Leader_Follower::FOLLOWER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tPDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderFluid, "IncompressibleFluid", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerFluid, "IncompressibleFluid", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                         { 1.0, 0.0 },
                         { 1.0, 1.0 },
                         { 0.0, 1.0 }};

               // set velocity dof types
               tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };

               // set viscosity property parameters
               tPropLeaderKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
               tPropFollowerKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
               break;
            }
            case 3 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0, 0.0 },
                         { 1.0, 0.0, 0.0 },
                         { 1.0, 1.0, 0.0 },
                         { 0.0, 1.0, 0.0 },
                         { 0.0, 0.0, 1.0 },
                         { 1.0, 0.0, 1.0 },
                         { 1.0, 1.0, 1.0 },
                         { 0.0, 1.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };

                // set viscosity property parameters
                tPropLeaderKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
                tPropFollowerKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
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
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tCMLeaderFluid->set_space_dim( iSpaceDim );
        tCMFollowerFluid->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
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
            int tNumDofVel  = tNumCoeff * iSpaceDim;
            int tNumDofP    = tNumCoeff;
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator viscosity
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatVis );

            // fill coefficients for leader FI
            Matrix< DDRMat > tFollowerDOFHatVel;;
            fill_uhat( tFollowerDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tFollowerDOFHatP;
            fill_phat( tFollowerDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tFollowerDOFHatVis;
            fill_phat( tFollowerDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVel );

            // create the field interpolator pressure
            tFollowerFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tFollowerFIs( 1 )->set_coeff( tFollowerDOFHatP );

            // create the field interpolator pressure
            tFollowerFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tFollowerFIs( 2 )->set_coeff( tFollowerDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofP + tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 },
                    { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };
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
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
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
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
}/*END_TEST_CASE*/

