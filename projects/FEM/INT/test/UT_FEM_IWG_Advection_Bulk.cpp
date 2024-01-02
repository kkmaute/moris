/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Advection_Bulk.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"
//MTK/src
#include "cl_MTK_Enums.hpp"
//LINALG/src
#include "op_equal_equal.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
//FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_MTK_Integrator.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Advection_Bulk", "[IWG_Advection_Bulk]" )
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
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tTempDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );
    //    tPropHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropHeatCapacity->set_val_function( tTEMPFIValFunc );
    //    tPropHeatCapacity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 10.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropDensity->set_val_function( tTEMPFIValFunc );
    //    tPropDensity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 5.0 }} } );
    tPropConductivity->set_val_function( tConstValFunc );
    //    tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropConductivity->set_val_function( tTEMPFIValFunc );
    //    tPropConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLoad = std::make_shared< fem::Property >();
    tPropLoad->set_parameters( { {{ 1.5 }} } );
    tPropLoad->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropAdvectionSource = std::make_shared< fem::Property >();
    tPropAdvectionSource->set_parameters( { {{ 2.0 * 1.5 }} } );
    tPropAdvectionSource->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMDiffusion =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMDiffusion->set_dof_type_list( { tTempDofTypes } );
    tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMDiffusion->set_property( tPropDensity,      "Density" );
    tCMDiffusion->set_property( tPropHeatCapacity, "HeatCapacity" );
    tCMDiffusion->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPSUPG =
            tSPFactory.create_SP( fem::Stabilization_Type::SUPG_ADVECTION );
    tSPSUPG->set_dof_type_list( { tVelDofTypes } );
    tSPSUPG->set_property( tPropConductivity,    "Conductivity", mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropDensity,         "Density",      mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropHeatCapacity,    "HeatCapacity", mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropAdvectionSource, "Source",       mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPSUPG->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tIWG->set_stabilization_parameter( tSPSUPG, "SUPG" );
    tIWG->set_property(tPropLoad, "Load", mtk::Leader_Follower::LEADER);

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                        { 1.0, 0.0 },
                        { 1.0, 1.0 },
                        { 0.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };
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
        tCMDiffusion->set_space_dim( iSpaceDim );
        tSPSUPG->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder   = tIntegrationOrders( iSpaceDim - 2 );

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
            int tNumDofTEMP = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_phat( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator temperature
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatTEMP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };

            // set size and fill the set Jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofTEMP,
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    tNumDofVel + tNumDofTEMP,
                    tNumDofVel + tNumDofTEMP,
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv,tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

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
}

TEST_CASE( "IWG_Advection_PhaseChange_Bulk", "[IWG_Advection_PhaseChange_Bulk]" )
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
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tTempDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );
    //    tPropHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropHeatCapacity->set_val_function( tTEMPFIValFunc );
    //    tPropHeatCapacity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 10.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropDensity->set_val_function( tTEMPFIValFunc );
    //    tPropDensity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 5.0 }} } );
    tPropConductivity->set_val_function( tConstValFunc );
    //    tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropConductivity->set_val_function( tTEMPFIValFunc );
    //    tPropConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLatentHeat = std::make_shared< fem::Property >();
    tPropLatentHeat->set_parameters( { {{ 22.0 }} } );
    tPropLatentHeat->set_val_function( tConstValFunc );

    // Phase change temperature chosen such that test operates within phase change regime
    std::shared_ptr< fem::Property > tPropPCTemp = std::make_shared< fem::Property >();
    tPropPCTemp->set_parameters( { {{ 32.468 }} } );
    tPropPCTemp->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropPhaseState = std::make_shared< fem::Property >();
    tPropPhaseState->set_parameters( { {{ 2.0 }} } );
    tPropPhaseState->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropPCconst = std::make_shared< fem::Property >();
    tPropPCconst->set_parameters( { {{ 0.6 }} } );
    tPropPCconst->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropLoad = std::make_shared< fem::Property >();
    tPropLoad->set_parameters( { {{ 1.5 }} } );
    tPropLoad->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropAdvectionSource = std::make_shared< fem::Property >();
    tPropAdvectionSource->set_parameters( { {{ 2.0 * 1.5 }} } );
    tPropAdvectionSource->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMDiffusion =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO_PC );
    tCMDiffusion->set_dof_type_list( { tTempDofTypes } );
    tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMDiffusion->set_property( tPropDensity,      "Density" );
    tCMDiffusion->set_property( tPropHeatCapacity, "HeatCapacity" );
    tCMDiffusion->set_property( tPropLatentHeat,   "LatentHeat" );
    tCMDiffusion->set_property( tPropPCTemp,       "PCTemp" );
    tCMDiffusion->set_property( tPropPhaseState,   "PhaseStateFunction" );
    tCMDiffusion->set_property( tPropPCconst,      "PhaseChangeConst" );
    tCMDiffusion->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPSUPG =
            tSPFactory.create_SP( fem::Stabilization_Type::SUPG_ADVECTION );
    tSPSUPG->set_dof_type_list( tDofTypes );
    tSPSUPG->set_property( tPropConductivity,    "Conductivity",       mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropDensity,         "Density",            mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropHeatCapacity,    "HeatCapacity",       mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropLatentHeat,      "LatentHeat",         mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropPCTemp,          "PCTemp",             mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropPhaseState,      "PhaseStateFunction", mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropPCconst,         "PhaseChangeConst",   mtk::Leader_Follower::LEADER );
    tSPSUPG->set_property( tPropAdvectionSource, "Source",             mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPSUPG->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tIWG->set_stabilization_parameter( tSPSUPG, "SUPG" );
    tIWG->set_property(tPropLoad, "Load", mtk::Leader_Follower::LEADER);

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                        { 1.0, 0.0 },
                        { 1.0, 1.0 },
                        { 0.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };
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
        tCMDiffusion->set_space_dim( iSpaceDim );
        tSPSUPG->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder   = tIntegrationOrders( iSpaceDim - 2 );

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
            int tNumDofTEMP = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_phat( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator temperature
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatTEMP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };

            // set size and fill the set Jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofTEMP,
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    tNumDofVel + tNumDofTEMP,
                    tNumDofVel + tNumDofTEMP,
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv,tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop over integartion points
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
}

TEST_CASE( "IWG_Advection_Bulk_YZBeta", "[IWG_Advection_Bulk_YZBeta]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-5;

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
    Vector< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tTempDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );
    //    tPropHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropHeatCapacity->set_val_function( tTEMPFIValFunc );
    //    tPropHeatCapacity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 10.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropDensity->set_val_function( tTEMPFIValFunc );
    //    tPropDensity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 5.0 }} } );
    tPropConductivity->set_val_function( tConstValFunc );
    //    tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropConductivity->set_val_function( tTEMPFIValFunc );
    //    tPropConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLoad = std::make_shared< fem::Property >();
    tPropLoad->set_parameters( { {{ 1.5 }} } );
    tPropLoad->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropAdvectionSource = std::make_shared< fem::Property >();
    tPropAdvectionSource->set_parameters( { {{ 2.0 * 1.5 }} } );
    tPropAdvectionSource->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropBeta = std::make_shared< fem::Property >();
    tPropBeta->set_parameters( { {{ 1.0 }} } );
    tPropBeta->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropRefState = std::make_shared< fem::Property >();
    tPropRefState->set_parameters( { {{ 1.0 }} } );
    tPropRefState->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMDiffusion =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMDiffusion->set_dof_type_list( { tTempDofTypes } );
    tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMDiffusion->set_property( tPropDensity,      "Density" );
    tCMDiffusion->set_property( tPropHeatCapacity, "HeatCapacity" );
    tCMDiffusion->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPYZBeta =
            tSPFactory.create_SP( fem::Stabilization_Type::YZBETA_ADVECTION );
    tSPYZBeta->set_dof_type_list( { tTempDofTypes } );
    tSPYZBeta->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tSPYZBeta->set_property( tPropBeta,      "Beta",           mtk::Leader_Follower::LEADER );
    tSPYZBeta->set_property( tPropRefState,  "ReferenceState", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPYZBeta->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tIWG->set_stabilization_parameter( tSPYZBeta, "YZBeta" );
    tIWG->set_property(tPropLoad, "Load", mtk::Leader_Follower::LEADER);

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                        { 1.0, 0.0 },
                        { 1.0, 1.0 },
                        { 0.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };
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
        tCMDiffusion->set_space_dim( iSpaceDim );
        tSPYZBeta->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder   = tIntegrationOrders( iSpaceDim - 2 );

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
            int tNumDofTEMP = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_phat( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator temperature
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatTEMP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };

            // set size and fill the set Jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofTEMP,
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    tNumDofVel + tNumDofTEMP,
                    tNumDofVel + tNumDofTEMP,
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv,tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

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
}

TEST_CASE( "IWG_Advection_Bulk_Crosswind", "[IWG_Advection_Bulk_Crosswind]" )
{
    // define an epsilon environment
    real tEpsilon = 5E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-5;

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
    Vector< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tTempDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );
    //    tPropHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropHeatCapacity->set_val_function( tTEMPFIValFunc );
    //    tPropHeatCapacity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 10.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropDensity->set_val_function( tTEMPFIValFunc );
    //    tPropDensity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 5.0 }} } );
    tPropConductivity->set_val_function( tConstValFunc );
    //    tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropConductivity->set_val_function( tTEMPFIValFunc );
    //    tPropConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLoad = std::make_shared< fem::Property >();
    tPropLoad->set_parameters( { {{ 1.5 }} } );
    tPropLoad->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropAdvectionSource = std::make_shared< fem::Property >();
    tPropAdvectionSource->set_parameters( { {{ 2.0 * 1.5 }} } );
    tPropAdvectionSource->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropBeta = std::make_shared< fem::Property >();
    tPropBeta->set_parameters( { {{ 1.0 }} } );
    tPropBeta->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropRefState = std::make_shared< fem::Property >();
    tPropRefState->set_parameters( { {{ 1.0 }} } );
    tPropRefState->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMDiffusion =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMDiffusion->set_dof_type_list( { tTempDofTypes } );
    tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMDiffusion->set_property( tPropDensity,      "Density" );
    tCMDiffusion->set_property( tPropHeatCapacity, "HeatCapacity" );
    tCMDiffusion->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPCrosswind =
            tSPFactory.create_SP( fem::Stabilization_Type::CROSSWIND );
    tSPCrosswind->set_parameters( { { { 0.7 } }, { { 1.0e-6 } } } );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPCrosswind->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tIWG->set_stabilization_parameter( tSPCrosswind, "DiffusionCrosswind" );
    tIWG->set_property(tPropLoad, "Load", mtk::Leader_Follower::LEADER);

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                        { 1.0, 0.0 },
                        { 1.0, 1.0 },
                        { 0.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };
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
        tCMDiffusion->set_space_dim( iSpaceDim );
        tSPCrosswind->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder   = tIntegrationOrders( iSpaceDim - 2 );

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
            int tNumDofTEMP = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_phat( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator temperature
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatTEMP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };

            // set size and fill the set Jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofTEMP,
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    tNumDofVel + tNumDofTEMP,
                    tNumDofVel + tNumDofTEMP,
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv,tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

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
}

TEST_CASE( "IWG_Advection_Bulk_Isotropic_Diffusion", "[IWG_Advection_Bulk_Isotropic_Diffusion]" )
{
    // define an epsilon environment
    real tEpsilon = 5E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-5;

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
    Vector< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };

    Vector< Vector< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tTempDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );
    //    tPropHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropHeatCapacity->set_val_function( tTEMPFIValFunc );
    //    tPropHeatCapacity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 10.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropDensity->set_val_function( tTEMPFIValFunc );
    //    tPropDensity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 5.0 }} } );
    tPropConductivity->set_val_function( tConstValFunc );
    //    tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropConductivity->set_val_function( tTEMPFIValFunc );
    //    tPropConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLoad = std::make_shared< fem::Property >();
    tPropLoad->set_parameters( { {{ 1.5 }} } );
    tPropLoad->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropAdvectionSource = std::make_shared< fem::Property >();
    tPropAdvectionSource->set_parameters( { {{ 2.0 * 1.5 }} } );
    tPropAdvectionSource->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropBeta = std::make_shared< fem::Property >();
    tPropBeta->set_parameters( { {{ 1.0 }} } );
    tPropBeta->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropRefState = std::make_shared< fem::Property >();
    tPropRefState->set_parameters( { {{ 1.0 }} } );
    tPropRefState->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMDiffusion =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMDiffusion->set_dof_type_list( { tTempDofTypes } );
    tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMDiffusion->set_property( tPropDensity,      "Density" );
    tCMDiffusion->set_property( tPropHeatCapacity, "HeatCapacity" );
    tCMDiffusion->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPCrosswind =
            tSPFactory.create_SP( fem::Stabilization_Type::CROSSWIND );
    tSPCrosswind->set_parameters( { { { 0.7 } }, { { 1.0e-6 } } } );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPCrosswind->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tIWG->set_stabilization_parameter( tSPCrosswind, "DiffusionIsotropic" );
    tIWG->set_property(tPropLoad, "Load", mtk::Leader_Follower::LEADER);

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

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

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                        { 1.0, 0.0 },
                        { 1.0, 1.0 },
                        { 0.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };
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
        tCMDiffusion->set_space_dim( iSpaceDim );
        tSPCrosswind->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder   = tIntegrationOrders( iSpaceDim - 2 );

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
            int tNumDofTEMP = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_phat( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator temperature
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatTEMP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };

            // set size and fill the set Jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofTEMP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofTEMP,
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    tNumDofVel + tNumDofTEMP,
                    tNumDofVel + tNumDofTEMP,
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv,tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

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
}

