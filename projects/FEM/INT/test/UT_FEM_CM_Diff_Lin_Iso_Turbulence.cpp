/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM_Diff_Lin_Iso_Turbulence.cpp
 *
 */

#include "catch.hpp"

#define protected public
#define private   public
// FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private

//LINALG/src
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "fn_FEM_Check.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"
#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "CM_Diff_Lin_Iso_Turbulence", "[CM_Diff_Lin_Iso_Turbulence]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-5;

    // define a perturbation relative size
    real tPerturbation = 1.0E-6;

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
    Vector< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tTempDofTypes, tVisDofTypes };

    // create the properties
    std::shared_ptr< fem::Property > tPropKinViscosity = std::make_shared< fem::Property >();
    tPropKinViscosity->set_parameters( { { { 1.0 } } } );
    tPropKinViscosity->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropTurbPrandtl = std::make_shared< fem::Property >();
    tPropTurbPrandtl->set_parameters( { { { 0.7 } } } );
    tPropTurbPrandtl->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { { { 2.0 } } } );
    tPropConductivity->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { { { 3.0 } } } );
    tPropHeatCapacity->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO_TURBULENCE_SPALART_ALLMARAS );
    tCMLeader->set_dof_type_list( { tTempDofTypes, tVisDofTypes } );
    tCMLeader->set_property( tPropKinViscosity, "KinematicViscosity" );
    tCMLeader->set_property( tPropTurbPrandtl, "TurbulentPrandtl" );
    tCMLeader->set_property( tPropDensity, "Density" );
    tCMLeader->set_property( tPropConductivity, "Conductivity" );
    tCMLeader->set_property( tPropHeatCapacity, "HeatCapacity" );

    tCMLeader->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 0;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // set size and populate the set leader dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 0;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

    // standard diffusion for comparison
    std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiff =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMLeaderDiff->set_dof_type_list( { tTempDofTypes } );
    tCMLeaderDiff->set_property( tPropDensity,      "Density" );
    tCMLeaderDiff->set_property( tPropConductivity, "Conductivity" );
    tCMLeaderDiff->set_property( tPropHeatCapacity, "HeatCapacity" );

    tCMLeaderDiff->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet2 = new fem::Set();
    tCMLeaderDiff->set_set_pointer( static_cast< fem::Set* >( tSet2 ) );

    // set size for the set EqnObjDofTypeList
    tCMLeaderDiff->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeaderDiff->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderDiff->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 0;

    // set size and populate the set leader dof type map
    tCMLeaderDiff->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderDiff->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )        = 0;

    // build global dof type list
    tCMLeaderDiff->get_global_dof_type_list();

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( 1, 1, 10.0 );

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
        mtk::Interpolation_Rule tGIRule(
                tGeometryType,
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

        // set space dimensions for property, CM and SP
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeaderDiff->set_space_dim( iSpaceDim );

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
                    mtk::Integration_Order::BAR_2 );

            // create an integrator
            mtk::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule (
                    tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatTemp;
            fill_phat( tLeaderDOFHatTemp, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatViscosity;
            fill_phat( tLeaderDOFHatViscosity, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatTemp );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatViscosity );

            // create a field interpolator manager
            Vector< Vector< gen::PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeader->mSet->mLeaderFIManager     = &tFIManager;
            tCMLeaderDiff->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMLeader->set_field_interpolator_manager( &tFIManager );
            tCMLeaderDiff->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMLeader->reset_eval_flags();
                tCMLeaderDiff->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tCMLeaderDiff->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // populate the requested leader dof type for CM
                Vector< Vector< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeader->get_global_dof_type_list();
                Vector< Vector< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes2 =
                                       tCMLeaderDiff->get_global_dof_type_list();

                // populate the test leader dof type for CM
                Vector< Vector< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();
                Vector< Vector< MSI::Dof_Type > > tLeaderDofTypes2 =
                                       tCMLeaderDiff->get_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    const Vector< MSI::Dof_Type >& tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // get turbulence specific CM
                    CM_Diffusion_Linear_Isotropic_Turbulence* tCMLeaderPtr =
                            dynamic_cast< CM_Diffusion_Linear_Isotropic_Turbulence* >( tCMLeader.get() );

                    // effective conductivity
                    //------------------------------------------------------------------------------
                    // evaluate analytical derivative
                    Matrix< DDRMat > tdeffconddu = tCMLeaderPtr->deffconddu( tDofDerivative );

                    // evaluate derivative by FD
                    Matrix< DDRMat > tdeffcondduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::EFF_COND,
                            tdeffcondduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check deffconddu analytical and FD match " << std::endl;
                    // print( tdeffconddu, "tdeffconddu" );
                    // print( tdeffcondduFD, "tdeffcondduFD" );
                    bool tCheckEffCond = fem::check( tdeffconddu, tdeffcondduFD, tEpsilon );
                    REQUIRE( tCheckEffCond );

                    // effective conductivity space der
                    //------------------------------------------------------------------------------
                    // evaluate analytical derivative
                    Matrix< DDRMat > tdeffconddxdu = tCMLeaderPtr->deffconddxdu( tDofDerivative, 1 );

                    // evaluate derivative by FD
                    Matrix< DDRMat > tdeffconddxduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::EFF_COND_SPACE_DER,
                            tdeffconddxduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check tdeffconddxdu analytical and FD match " << std::endl;
                    // print( tdeffconddxdu, "tdeffconddxdu" );
                    // print( tdeffconddxduFD, "tdeffconddxduFD" );
                    bool tCheckEffCondSpaceDer = fem::check( tdeffconddxdu, tdeffconddxduFD, tEpsilon );
                    REQUIRE( tCheckEffCondSpaceDer );

                    // test effective conductivity
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    Vector< Vector< MSI::Dof_Type > > tTestDofTypes = { tTempDofTypes };
                    for ( uint iTestDof = 0; iTestDof < tTestDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        const Vector< MSI::Dof_Type >& tDofTest = tTestDofTypes( iTestDof );
                        // std::cout << "tDofTest " << static_cast< uint >( tDofTest( 0 ) ) << std::endl;

                        // evaluate analytical derivative
                        Matrix< DDRMat > tdtesteffconddu = tCMLeaderPtr->dtesteffconddu(
                                tDofDerivative,
                                tDofTest );

                        // evaluate derivative by FD
                        Matrix< DDRMat > tdtesteffcondduFD;
                        tCMLeader->eval_derivative_FD(
                                CM_Request_Type::TEST_EFF_COND,
                                tdtesteffcondduFD,
                                tDofDerivative,
                                tPerturbation,
                                tDofTest,
                                tNormal,
                                tJump );

                        // check that analytical and FD match
                        // std::cout << "check tdtesteffconddu analytical and FD match " << std::endl;
                        // print( tdtesteffconddu, "tdtesteffconddu" );
                        // print( tdtesteffcondduFD, "tdtesteffcondduFD" );
                        bool tCheckTestEffCond = fem::check( tdtesteffconddu, tdtesteffcondduFD, tEpsilon );
                        REQUIRE( tCheckTestEffCond );
                    }

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMLeader->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::FLUX,
                            tdfluxduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    bool tCheckFluxFluid = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckFluxFluid );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMLeader->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TRACTION,
                            tdtractionduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    bool tCheckTraction = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckTraction );

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for ( uint iTestDof = 0; iTestDof < tTestDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        const Vector< MSI::Dof_Type >& tDofTest = tTestDofTypes( iTestDof );

                        // evaluate dtesttractiondu
                        Matrix< DDRMat > tdtesttractiondu = tCMLeader->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMLeader->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD,
                                tPerturbation,
                                tNormal );

                        // check that analytical and FD match
                        bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                        REQUIRE( tCheckTestTraction );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMLeader->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::DIV_FLUX,
                            tddivfluxduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    bool tCheckDivFlux = fem::check( tddivfluxdu, tddivfluxduFD, tEpsilon );
                    REQUIRE( tCheckDivFlux );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/

