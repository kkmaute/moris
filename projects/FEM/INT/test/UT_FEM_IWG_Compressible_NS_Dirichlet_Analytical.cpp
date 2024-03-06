/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Compressible_NS_Dirichlet_Analytical.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_IWG_Compressible_NS_Dirichlet_Nitsche.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
//LINALG/src
#include "op_equal_equal.hpp"
//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_MM_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Fields_for_NS_Compressible_UT_for_BD.cpp"
#include "FEM_Test_Proxy/cl_FEM_Flux_Matrix_Reference_Values_Nitsche.cpp"

#include "fn_trans.hpp"
#include "fn_FEM_Check.hpp"

// debug - output to hdf5
#include "paths.hpp"
#include "HDF5_Tools.hpp"
#include "fn_FEM_Check.hpp"
#include "fn_sqrtmat.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Compressible_NS_Dirichlet_Analytical",
        "[IWG_Compressible_NS_Dirichlet_Analytical]" )
{
    // define an epsilon environment (max allowed relative error)
    real tEpsilon = 1.0E-10;

    // define absolute tolerance accepted as numerical error
    real tAbsTol = 1.0E-12;

    // weight for the GLS term, use 0 for turning GLS off, use 1 for weak form + regular GLS formulation
    real tUpwindFactor = 0.0;

    // penalty factors
    real  tPpenatly =   0.0;
    real tUXpenatly = 100.0;
    real tUYpenatly = 100.0;
    real  tTpenatly =   0.0;

    // prescribed values
    real  tPpresc =   7.8366834e+04;
    real tUXpresc =   0.0;
    real tUYpresc =   0.0;
    real  tTpresc = 273.0;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // dof type list
    Vector< MSI::Dof_Type > tPressureDof = { MSI::Dof_Type::P };
    Vector< MSI::Dof_Type > tVelocityDof = { MSI::Dof_Type::VX };
    Vector< MSI::Dof_Type > tTempDof     = { MSI::Dof_Type::TEMP };

    Vector< Vector< MSI::Dof_Type > > tDofTypes         = { tPressureDof, tVelocityDof, tTempDof };
    Vector< Vector< MSI::Dof_Type > > tResidualDofTypes = tDofTypes;

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties

    // dynamic viscosity
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    real tMu = 1.716e-5;
    tPropViscosity->set_parameters( { {{ tMu }} } );
    tPropViscosity->set_val_function( tConstValFunc );

    // isochoric heat capacity
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 0.718e3 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );

    // specific gas constant
    std::shared_ptr< fem::Property > tPropGasConstant = std::make_shared< fem::Property >();
    tPropGasConstant->set_parameters( { {{ 287.058 }} } );
    tPropGasConstant->set_val_function( tConstValFunc );

    // thermal conductivity
    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 24.35e-3 }} } );
    tPropConductivity->set_val_function( tConstValFunc );

    // Body heat load
    std::shared_ptr< fem::Property > tPropHeatLoad = std::make_shared< fem::Property >();
    tPropHeatLoad->set_parameters( { {{ 1.0e+5 }} } );
    tPropHeatLoad->set_val_function( Func_Heat_Load_Distribution );

    //------------------------------------------------------------------------------
    // prescribed values

    // prescribe pressure
    std::shared_ptr< fem::Property > tPropPrescPres = std::make_shared< fem::Property >();
    tPropPrescPres->set_parameters( { {{ tPpresc }} } );
    tPropPrescPres->set_val_function( tConstValFunc );

    // prescribed velocity
    std::shared_ptr< fem::Property > tPropPrescVel = std::make_shared< fem::Property >();
    tPropPrescVel->set_parameters( { {{ tUXpresc },{ tUYpresc }} } );
    tPropPrescVel->set_val_function( tConstValFunc );

    // velocity-select-matrix
    std::shared_ptr< fem::Property > tPropVelSelect = std::make_shared< fem::Property >();
    tPropVelSelect->set_parameters( { {{ 0.0, 0.0 },{ 0.0, 1.0 }} } );
    tPropVelSelect->set_val_function( tConstValFunc );

    // prescribed temperature
    std::shared_ptr< fem::Property > tPropPrescTemp = std::make_shared< fem::Property >();
    tPropPrescTemp->set_parameters( { {{ tTpresc }} } );
    tPropPrescTemp->set_val_function( tConstValFunc );

    // upwinding weight
    std::shared_ptr< fem::Property > tPropUpwind = std::make_shared< fem::Property >();
    tPropUpwind->set_parameters( { {{ tUpwindFactor }} } );
    tPropUpwind->set_val_function( tConstValFunc );

    // empty property
    std::shared_ptr< fem::Property > tPropEmpty = nullptr;

    // define material model and assign properties
    fem::MM_Factory tMMFactory;

    std::shared_ptr< fem::Material_Model > tMMFluid =
            tMMFactory.create_MM( fem::Material_Type::PERFECT_GAS );
    tMMFluid->set_dof_type_list( {tPressureDof, tTempDof } );
    tMMFluid->set_property( tPropHeatCapacity, "IsochoricHeatCapacity" );
    tMMFluid->set_property( tPropGasConstant,  "SpecificGasConstant" );

    // define constitutive model and assign properties
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_COMPRESSIBLE_NEWTONIAN );
    tCMLeaderFluid->set_dof_type_list( {tPressureDof, tVelocityDof, tTempDof } );
    tCMLeaderFluid->set_property( tPropViscosity,    "DynamicViscosity" );
    tCMLeaderFluid->set_property( tPropConductivity, "ThermalConductivity" );
    tCMLeaderFluid->set_material_model( tMMFluid, "ThermodynamicMaterialModel" );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::COMPRESSIBLE_DIRICHLET_NITSCHE );
    tSPNitsche->set_property( tPropViscosity, "DynamicViscosity", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_property( tPropConductivity, "ThermalConductivity", mtk::Leader_Follower::LEADER );
    Vector< Matrix< DDRMat > > tNitscheParams = { {{tPpenatly},{tUXpenatly/tMu},{tUYpenatly/tMu},{tTpenatly}} };
    tSPNitsche->set_parameters( tNitscheParams );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE );

    tIWG->set_residual_dof_type( tResidualDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropViscosity,    "DynamicViscosity" );
    tIWG->set_property( tPropConductivity, "ThermalConductivity" );
    tIWG->set_material_model( tMMFluid, "FluidMM" );
    tIWG->set_constitutive_model( tCMLeaderFluid, "FluidCM" );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitschePenaltyParameter" );

    // tIWG->set_property( tPropPrescPres, "PrescribedDof1" );
    // //tIWG->set_property( tPropPrescVel,  "PrescribedVelocity" );
    // //tIWG->set_property( tPropPrescVel,  "tPropVelSelect" );
    // //tIWG->set_property( tPropPrescTemp, "PrescribedDof3" );
    // tIWG->set_property( tPropUpwind,    "PressureUpwind" );

    //------------------------------------------------------------------------------
    // set a fem set pointer

    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tMMFluid->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tCMLeaderFluid->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tIWG->set_fem_set( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tMMFluid->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tMMFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )   = 0;
    tMMFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 1;

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )     = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )     = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

    // set number of spatial dimensions
    uint iSpaceDim = 2;

    // set geometry type
    tGeometryType = mtk::Geometry_Type::QUAD;

    // set velocity dof types
    tVelocityDof = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

    // set space dimension to CM
    tMMFluid->set_space_dim( iSpaceDim );
    tCMLeaderFluid->set_space_dim( iSpaceDim );
    tSPNitsche->set_space_dim( iSpaceDim );

    // // create normal for Nitsche
    // Matrix< DDRMat > tNormal = {{3.0/5.0},{4.0/5.0}};

    // create normal for Pressure
    Matrix< DDRMat > tNormal = {{-1.0},{0.0}};

    // assign normal to IWG
    tIWG->set_normal( tNormal );

    // set interpolation order
    // uint iInterpOrder = 2;

    // create an interpolation order
    mtk::Interpolation_Order tGIInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

    //------------------------------------------------------------------------------
    // prepare element information

    // initialize matrices to be filled
    Matrix< DDRMat >tXHat;
    Matrix< DDRMat >tTHat;
    Matrix< DDRMat >tLeaderDOFHatP;
    Matrix< DDRMat >tLeaderDOFHatVel;
    Matrix< DDRMat >tLeaderDOFHatTemp;

    // fill in data for element size
    fill_data_Nitsche_element( tXHat, tTHat );

    // fill DoF values
    fill_Nitsche_PHat( tLeaderDOFHatP );
    //tLeaderDOFHatP = trans( tLeaderDOFHatP );
    fill_Nitsche_UHat( tLeaderDOFHatVel );
    fill_Nitsche_THat( tLeaderDOFHatTemp );
    //tLeaderDOFHatTemp = trans( tLeaderDOFHatTemp );

    //------------------------------------------------------------------------------
    // space and time geometry interpolators
    // create a space geometry interpolation rule
    mtk::Interpolation_Rule tGIRule( tGeometryType,
            mtk::Interpolation_Type::LAGRANGE,
            tGIInterpolationOrder,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

    // set the coefficients xHat, tHat
    tGI.set_coeff( tXHat, tTHat );

    //------------------------------------------------------------------------------
    // integration points
    // get an integration order
    mtk::Integration_Order tIntegrationOrder = mtk::Integration_Order::QUAD_3x3;

    // create an integration rule
    mtk::Integration_Rule tIntegrationRule(
            tGeometryType,
            mtk::Integration_Type::GAUSS,
            tIntegrationOrder,
            mtk::Geometry_Type::LINE,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_3 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
    tIntegrator.get_weights( tIntegWeights );

    //------------------------------------------------------------------------------
    // field interpolators
    // create an interpolation order
    mtk::Interpolation_Order tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

    // number of dof for interpolation order
    uint tNumCoeff = 18;

    // get number of dof per type
    int tNumDofP  = tNumCoeff;
    int tNumDofVel  = tNumCoeff * iSpaceDim;
    int tNumDofTemp = tNumCoeff;
    int tTotalNumDof = tNumDofP + tNumDofVel + tNumDofTemp;

    //create a space time interpolation rule
    mtk::Interpolation_Rule tFIRule (
            tGeometryType,
            mtk::Interpolation_Type::LAGRANGE,
            tInterpolationOrder,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a cell of field interpolators for IWG
    Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

    // create the field interpolator density
    tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPressureDof );
    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatP );

    // create the field interpolator velocity
    tLeaderFIs( 1 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDof );
    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVel );

    // create the field interpolator pressure
    tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDof );
    tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatTemp );

    // set size and fill the set residual assembly map
    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofP - 1 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofP, tNumDofP + tNumDofVel - 1 } };
    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tTotalNumDof - 1 } };

    // set size and fill the set jacobian assembly map
    Matrix< DDSMat > tJacAssembly = {
            { 0, tNumDofP - 1 },
            { tNumDofP, tNumDofP + tNumDofVel - 1 },
            { tNumDofP + tNumDofVel, tTotalNumDof - 1 } };
    tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
    tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
    tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;

    // set size and init the set residual and jacobian
    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size(
            tTotalNumDof,
            1,
            0.0 );
    tIWG->mSet->mJacobian.set_size(
            tTotalNumDof,
            tTotalNumDof,
            0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested leader dof type
    tIWG->mLeaderSideInfo.mRequestedGlobalDofTypes = tDofTypes;

    // create a field interpolator manager
    Vector< Vector< enum gen::PDV_Type > > tDummyDv;
    Vector< Vector< enum mtk::Field_Type > > tDummyField;
    Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tLeaderFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIWG->mSet->mLeaderFIManager = &tFIManager;

    // set IWG field interpolator manager
    tIWG->set_field_interpolator_manager( &tFIManager );

    // set the interpolator manager to the set
    tCMLeaderFluid->mSet->mLeaderFIManager = &tFIManager;

    // set IWG field interpolator manager
    tCMLeaderFluid->set_field_interpolator_manager( &tFIManager );

    // create evaluation point xi, tau
    Matrix< DDRMat > tParamPoint = {
            { -1.0 },
            { +0.57735026918962 },
            { +0.57735026918962 } };

    // set integration point
    tCMLeaderFluid->mSet->mLeaderFIManager->set_space_time( tParamPoint );
    tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

    // check evaluation of the residual for IWG
    //------------------------------------------------------------------------------
    // reset residual & jacobian
    tIWG->mSet->mResidual( 0 ).fill( 0.0 );
    tIWG->mSet->mJacobian.fill( 0.0 );

    // get the child
    fem::IWG_Compressible_NS_Dirichlet_Nitsche * tChildIWG = dynamic_cast< fem::IWG_Compressible_NS_Dirichlet_Nitsche * > ( tIWG.get() );

    //------------------------------------------------------------------------------
    // Test Upwind Term

    // Set properties
    tIWG->set_property( tPropPrescPres, "PrescribedDof1" );
    tIWG->set_property( tPropUpwind,    "PressureUpwind" );

    // reset IWG
    tIWG->reset_eval_flags();
    tIWG->set_field_interpolator_manager( &tFIManager );

    // Check Upwind Term
    std::cout << "\nUpwind Term \n" << std::flush;
    Matrix< DDRMat > tResUpwind = -1.0 * tChildIWG->W_trans() * tChildIWG->UpwindOperator() * tChildIWG->jump();
    Matrix< DDRMat > tResUpwind_ref = get_reference_ResPrsUpwind();
    REQUIRE( fem::check( tResUpwind, tResUpwind_ref, tEpsilon, true, true, tAbsTol ) );

    //------------------------------------------------------------------------------
    // set eval point on Nitsche boundary

    // set integration point
    tParamPoint = {
            { +0.57735026918962 },
            { -1.0 },
            { +0.57735026918962 } };
    tCMLeaderFluid->mSet->mLeaderFIManager->set_space_time( tParamPoint );
    tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

    // assign normal to IWG
    tNormal = {{0.0},{-1.0}};
    tIWG->set_normal( tNormal );

    //------------------------------------------------------------------------------
    // Test Penalty Term

    // Set properties
    tIWG->set_property( tPropEmpty,     "PrescribedDof1" );
    tIWG->set_property( tPropPrescVel,  "PrescribedVelocity" );
    tIWG->set_property( tPropVelSelect, "SelectVelocity" );
    tIWG->set_property( tPropEmpty,     "PressureUpwind" );

    // reset IWG
    tIWG->reset_eval_flags();
    tIWG->set_field_interpolator_manager( &tFIManager );

    // Check Penalty Term
    std::cout << "\nPenalty Term \n" << std::flush;
    Matrix< DDRMat > tDiagSP = {
            { 0.0,        0.0,        0.0, 0.0 },
            { 0.0, tUXpenatly,        0.0, 0.0 },
            { 0.0,        0.0, tUYpenatly, 0.0 },
            { 0.0,        0.0,        0.0, 0.0 } };
    Matrix< DDRMat > tResPen = tChildIWG->W_trans() * tDiagSP * tChildIWG->jump();
    Matrix< DDRMat > tResPen_ref = get_reference_ResVelNitschePen();
    REQUIRE( fem::check( tResPen, tResPen_ref, tEpsilon, true, true, tAbsTol ) );

    //------------------------------------------------------------------------------
    // Test Consistency Term

    // Set properties
    tIWG->set_property( tPropEmpty,     "PrescribedDof1" );
    tIWG->set_property( tPropPrescVel,  "PrescribedVelocity" );
    tIWG->set_property( tPropVelSelect, "SelectVelocity" );
    tIWG->set_property( tPropEmpty,     "PressureUpwind" );

    // reset IWG
    tIWG->reset_eval_flags();
    tIWG->set_field_interpolator_manager( &tFIManager );

    // Consistency Term
    std::cout << "\nConsistency Term \n" << std::flush;
    Matrix< DDRMat > tResCon = -1.0 * tChildIWG->W_trans() * tChildIWG->select_matrix() * tChildIWG->Traction();
    Matrix< DDRMat > tResCon_ref = get_reference_ResVelNitscheCon();
    REQUIRE( fem::check( tResCon, tResCon_ref, tEpsilon, true, true, tAbsTol ) );

    //------------------------------------------------------------------------------
    // Test Adjoint Term

    // Set properties
    tIWG->set_property( tPropEmpty,     "PrescribedDof1" );
    tIWG->set_property( tPropPrescVel,  "PrescribedVelocity" );
    tIWG->set_property( tPropVelSelect, "SelectVelocity" );
    tIWG->set_property( tPropEmpty,     "PressureUpwind" );

    // reset IWG
    tIWG->reset_eval_flags();
    tIWG->set_field_interpolator_manager( &tFIManager );

    // Adjoint Term
    std::cout << "\nAdjoint Term \n" << std::flush;
    Matrix< DDRMat > tResAdj = tChildIWG->TestTraction() * tChildIWG->jump();
    Matrix< DDRMat > tResAdj_ref = get_reference_ResVelNitscheAdj();
    REQUIRE( fem::check( tResAdj, tResAdj_ref, tEpsilon, true, true, tAbsTol ) );

    //------------------------------------------------------------------------------
    // Test Whole Nitsche Residual

    // Set properties
    tIWG->set_property( tPropEmpty,     "PrescribedDof1" );
    tIWG->set_property( tPropPrescVel,  "PrescribedVelocity" );
    tIWG->set_property( tPropVelSelect, "SelectVelocity" );
    tIWG->set_property( tPropEmpty,     "PressureUpwind" );

    // reset IWG
    tIWG->reset_eval_flags();
    tIWG->set_field_interpolator_manager( &tFIManager );

    // Complete Residual
    std::cout << "\nComplete Residual \n" << std::flush;
    tIWG->compute_residual( 1.0 );
    Matrix< DDRMat > tRes = tIWG->mSet->get_residual()( 0 );
    Matrix< DDRMat > tRes_ref = get_reference_ResNitsche();
    REQUIRE( fem::check( tRes, tRes_ref, tEpsilon, true, true, tAbsTol ) );

    //------------------------------------------------------------------------------

    // clean up
    tLeaderFIs.clear();

}/*END_TEST_CASE*/

