/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM_Fluid.cpp
 *
 */

#include "catch.hpp"

#define protected public
#define private   public
//FEM/INT/src
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

using namespace moris;
using namespace fem;

TEST_CASE( "CM_Fluid", "[CM_Fluid]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-6;

    // define a perturbation relative size
    real tPerturbation = 2.0E-4;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = { tVelDofTypes( 0 ), tPDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    //tPropViscosity->set_parameters( { {{ 1.0 }}, {{ 0.0 },{0.0}} } );
    tPropViscosity->set_val_function( tConstValFunc );
    tPropViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );
    //tPropViscosity->set_dof_type_list( { tPDofTypes } );
    //tPropViscosity->set_val_function( tPFIValFunc );
    //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //tPropDensity->set_dof_type_list( { tPDofTypes } );
    //tPropDensity->set_val_function( tPFIValFunc );
    //tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMLeaderFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) } );
    tCMLeaderFluid->set_property( tPropViscosity, "Viscosity" );
    tCMLeaderFluid->set_property( tPropDensity, "Density" );
    tCMLeaderFluid->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMLeaderFluid->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeaderFluid->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeaderFluid->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMLeaderFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;

    // set size and populate the set leader dof type map
    tCMLeaderFluid->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderFluid->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMLeaderFluid->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;

    // build global dof type list
    tCMLeaderFluid->get_global_dof_type_list();

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( iSpaceDim, 1, 10.0 );

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

                // set viscosity space derivative parameters
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{ 0.0 },{0.0}} } );

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

                // set viscosity space derivative parameters
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );

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
        tCMLeaderFluid->set_space_dim( iSpaceDim );

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
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeaderFluid->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMLeaderFluid->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMLeaderFluid->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeaderFluid->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // populate the requested leader dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeaderFluid->get_global_dof_type_list();

                // populate the test leader dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeaderFluid->get_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMLeaderFluid->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMLeaderFluid->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckFluxFluid = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckFluxFluid );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate dstraindu
                    Matrix< DDRMat > tdstraindu = tCMLeaderFluid->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeaderFluid->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckStrainFluid = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckStrainFluid );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMLeaderFluid->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMLeaderFluid->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                    // check that analytical and FD match
                    bool tCheckTractionFluid = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckTractionFluid );

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for( uint iTestDof = 0; iTestDof < tLeaderDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( iTestDof );

                        // evaluate dtesttractiondu
                        Matrix< DDRMat > tdtesttractiondu = tCMLeaderFluid->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tJump,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMLeaderFluid->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD,
                                tPerturbation,
                                tNormal,
                                tJump );

                        // check that analytical and FD match
                        bool tCheckTestTractionFluid = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                        REQUIRE( tCheckTestTractionFluid );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMLeaderFluid->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMLeaderFluid->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivFluxFluid = fem::check( tddivfluxdu, tddivfluxduFD, tEpsilon );
                    REQUIRE( tCheckDivFluxFluid );

                    // div strain
                    //------------------------------------------------------------------------------
                    // evaluate ddivstraindu
                    Matrix< DDRMat > tddivstraindu = tCMLeaderFluid->ddivstraindu( tDofDerivative );

                    // evaluate ddivstraindu by FD
                    Matrix< DDRMat > tddivstrainduFD;
                    tCMLeaderFluid->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivStrainFluid = fem::check( tddivstraindu, tddivstrainduFD, tEpsilon );
                    REQUIRE( tCheckDivStrainFluid );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "CM_Laminar_With_Turbulence", "[CM_Laminar_With_Turbulence]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-6;

    // define a perturbation relative size
    real tPerturbation = 2.0E-4;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

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
    moris::Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_val_function( tConstValFunc );
    tPropViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );
    //tPropViscosity->set_dof_type_list( { tPDofTypes } );
    //tPropViscosity->set_val_function( tPFIValFunc );
    //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropKinViscosity = std::make_shared< fem::Property >();
    tPropKinViscosity->set_val_function( tConstValFunc );
    tPropKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //tPropDensity->set_dof_type_list( { tPDofTypes } );
    //tPropDensity->set_val_function( tPFIValFunc );
    //tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderLaminar =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMLeaderLaminar->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) } );
    tCMLeaderLaminar->set_property( tPropViscosity, "Viscosity" );
    tCMLeaderLaminar->set_property( tPropDensity, "Density" );
    tCMLeaderLaminar->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMLeaderTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMLeaderTurbulence->set_property( tPropViscosity, "Viscosity" );
    tCMLeaderTurbulence->set_property( tPropKinViscosity, "KinViscosity" );
    tCMLeaderTurbulence->set_property( tPropDensity, "Density" );
    tCMLeaderTurbulence->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMLeaderLaminar->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tCMLeaderTurbulence->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeaderTurbulence->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set leader dof type map
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P  ) )        = 1;
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // build global dof type list
    tCMLeaderLaminar->get_global_dof_type_list();
    tCMLeaderTurbulence->get_global_dof_type_list();

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( iSpaceDim, 1, 10.0 );

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
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );

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
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );

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

        // set space dimensions for property, CM and SP
        tCMLeaderLaminar->set_space_dim( iSpaceDim );
        tCMLeaderTurbulence->set_space_dim( iSpaceDim );

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

            // fill random coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );
            tLeaderDOFHatVis = 0.0 * tLeaderDOFHatVis;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator viscosity
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatVis );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeaderLaminar->mSet->mLeaderFIManager = &tFIManager;
            tCMLeaderTurbulence->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMLeaderLaminar->set_field_interpolator_manager( &tFIManager );
            tCMLeaderTurbulence->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMLeaderTurbulence->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeaderLaminar->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tCMLeaderTurbulence->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // populate the requested leader dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeaderLaminar->get_global_dof_type_list();

                // populate the test leader dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeaderLaminar->get_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMLeaderTurbulence->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMLeaderLaminar->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckFluxTurbulence = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckFluxTurbulence );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate dstraindu
                    Matrix< DDRMat > tdstraindu = tCMLeaderTurbulence->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeaderLaminar->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckStrainTurbulence = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckStrainTurbulence );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMLeaderTurbulence->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMLeaderLaminar->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                    // check that analytical and FD match
                    bool tCheckTractionTurbulence = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckTractionTurbulence );

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for( uint iTestDof = 0; iTestDof < tLeaderDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( iTestDof );

                        // evaluate dtesttractiondu
                        Matrix< DDRMat > tdtesttractiondu = tCMLeaderTurbulence->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tJump,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMLeaderLaminar->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD,
                                tPerturbation,
                                tNormal,
                                tJump );

                        // check that analytical and FD match
                        bool tCheckTestTractionTurbulence = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                        REQUIRE( tCheckTestTractionTurbulence );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMLeaderTurbulence->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMLeaderLaminar->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivFluxTurbulence = fem::check( tddivfluxdu, tddivfluxduFD, tEpsilon );
                    REQUIRE( tCheckDivFluxTurbulence );

                    // div strain
                    //------------------------------------------------------------------------------
                    // evaluate ddivstraindu
                    Matrix< DDRMat > tddivstraindu = tCMLeaderTurbulence->ddivstraindu( tDofDerivative );

                    // evaluate ddivstraindu by FD
                    Matrix< DDRMat > tddivstrainduFD;
                    tCMLeaderLaminar->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivStrainTurbulence = fem::check( tddivstraindu, tddivstrainduFD, tEpsilon );
                    REQUIRE( tCheckDivStrainTurbulence );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "CM_Laminar_Turbulence_Only", "[CM_Laminar_Turbulence_Only]" )
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
    moris::Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    //tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunc );
    tPropViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );
    //tPropViscosity->set_dof_type_list( { tPDofTypes } );
    //tPropViscosity->set_val_function( tPFIValFunc );
    //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropKinViscosity = std::make_shared< fem::Property >();
    tPropKinViscosity->set_val_function( tConstValFunc );
    tPropKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //tPropDensity->set_dof_type_list( { tPDofTypes } );
    //tPropDensity->set_val_function( tPFIValFunc );
    //tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderLaminar =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMLeaderLaminar->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) } );
    tCMLeaderLaminar->set_property( tPropViscosity, "Viscosity" );
    tCMLeaderLaminar->set_property( tPropDensity, "Density" );
    tCMLeaderLaminar->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMLeaderTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMLeaderTurbulence->set_property( tPropViscosity, "Viscosity" );
    tCMLeaderTurbulence->set_property( tPropKinViscosity, "KinViscosity" );
    tCMLeaderTurbulence->set_property( tPropDensity, "Density" );
    tCMLeaderTurbulence->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMLeaderLaminar->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tCMLeaderTurbulence->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeaderTurbulence->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set leader dof type map
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P  ) )        = 1;
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // build global dof type list
    tCMLeaderLaminar->get_global_dof_type_list();
    tCMLeaderTurbulence->get_global_dof_type_list();

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( iSpaceDim, 1, 10.0 );

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
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );

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
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );

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

        // set space dimensions for property, CM and SP
        tCMLeaderLaminar->set_space_dim( iSpaceDim );
        tCMLeaderTurbulence->set_space_dim( iSpaceDim );

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

            // fill random coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator viscosity
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatVis );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeaderLaminar->mSet->mLeaderFIManager = &tFIManager;
            tCMLeaderTurbulence->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMLeaderLaminar->set_field_interpolator_manager( &tFIManager );
            tCMLeaderTurbulence->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMLeaderLaminar->reset_eval_flags();
                tCMLeaderTurbulence->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeaderLaminar->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tCMLeaderTurbulence->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // populate the requested leader dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeaderTurbulence->get_global_dof_type_list();

                // populate the test leader dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeaderTurbulence->get_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMLeaderTurbulence->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMLeaderTurbulence->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                    if( tDofDerivative( 0 ) == MSI::Dof_Type::VX || tDofDerivative( 0 ) == MSI::Dof_Type::P )
                    {
                        // get the laminar part of dfluxdu
                        Matrix< DDRMat > tdfluxduLaminar = tCMLeaderLaminar->dFluxdDOF( tDofDerivative );

                        // get turbulent contribution
                        tdfluxdu   = tdfluxdu - tdfluxduLaminar;
                        tdfluxduFD = tdfluxduFD - tdfluxduLaminar;
                    }

                    // check that analytical and FD match
                    bool tCheckFluxTurbulence = fem::check(
                            tdfluxdu,
                            tdfluxduFD,
                            tEpsilon );
                    REQUIRE( tCheckFluxTurbulence );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate dstraindu
                    Matrix< DDRMat > tdstraindu = tCMLeaderTurbulence->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeaderTurbulence->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                    if( tDofDerivative( 0 ) == MSI::Dof_Type::VX || tDofDerivative( 0 ) == MSI::Dof_Type::P )
                    {
                        // evaluate laminar part of dstraindu
                        Matrix< DDRMat > tdstrainduLaminar = tCMLeaderLaminar->dStraindDOF( tDofDerivative );

                        // get turbulent part
                        tdstraindu = tdstraindu - tdstrainduLaminar;
                        tdstrainduFD = tdstrainduFD - tdstrainduLaminar;
                    }

                    // check that analytical and FD match
                    bool tCheckStrainTurbulence = fem::check(
                            tdstraindu,
                            tdstrainduFD,
                            tEpsilon );
                    REQUIRE( tCheckStrainTurbulence );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMLeaderTurbulence->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMLeaderTurbulence->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                    if( tDofDerivative( 0 ) == MSI::Dof_Type::VX || tDofDerivative( 0 ) == MSI::Dof_Type::P )
                    {
                        // evaluate laminar part of dtractiondu
                        Matrix< DDRMat > tdtractionduLaminar = tCMLeaderLaminar->dTractiondDOF( tDofDerivative, tNormal );

                        // get turbulent contribution
                        tdtractiondu = tdtractiondu - tdtractionduLaminar;
                        tdtractionduFD = tdtractionduFD - tdtractionduLaminar;
                    }

                    // check that analytical and FD match
                    bool tCheckTractionTurbulence = fem::check(
                            tdtractiondu,
                            tdtractionduFD,
                            tEpsilon );
                    REQUIRE( tCheckTractionTurbulence );

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for( uint iTestDof = 0; iTestDof < tLeaderDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( iTestDof );

                        // evaluate dtesttractiondu
                        Matrix< DDRMat > tdtesttractiondu = tCMLeaderTurbulence->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tJump,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMLeaderTurbulence->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD,
                                tPerturbation,
                                tNormal,
                                tJump );

                        if( tDofTest( 0 ) == MSI::Dof_Type::VX || tDofTest( 0 ) == MSI::Dof_Type::P )
                        {
                            if( tDofDerivative( 0 ) == MSI::Dof_Type::VX || tDofDerivative( 0 ) == MSI::Dof_Type::P )
                            {
                                Matrix< DDRMat > tdtesttractionduLaminar = tCMLeaderLaminar->dTestTractiondDOF(
                                        tDofDerivative,
                                        tNormal,
                                        tJump,
                                        tDofTest );

                                tdtesttractiondu   = tdtesttractiondu - tdtesttractionduLaminar;
                                tdtesttractionduFD = tdtesttractionduFD - tdtesttractionduLaminar;
                            }
                        }

                        // check that analytical and FD match
                        bool tCheckTestTractionTurbulence = fem::check(
                                tdtesttractiondu,
                                tdtesttractionduFD,
                                tEpsilon );
                        REQUIRE( tCheckTestTractionTurbulence );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMLeaderTurbulence->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMLeaderTurbulence->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                    if( tDofDerivative( 0 ) == MSI::Dof_Type::VX || tDofDerivative( 0 ) == MSI::Dof_Type::P )
                    {
                        // evaluate the laminar part of ddivfluxdu
                        Matrix< DDRMat > tddivfluxduLaminar = tCMLeaderLaminar->ddivfluxdu( tDofDerivative );

                        // get turbulent contribution
                        tddivfluxdu = tddivfluxdu - tddivfluxduLaminar;
                        tddivfluxduFD = tddivfluxduFD - tddivfluxduLaminar;
                    }

                    // check that analytical and FD match
                    bool tCheckDivFluxTurbulence = fem::check(
                            tddivfluxdu,
                            tddivfluxduFD,
                            tEpsilon );
                    REQUIRE( tCheckDivFluxTurbulence );

                    // div strain
                    //------------------------------------------------------------------------------
                    // evaluate ddivstraindu
                    Matrix< DDRMat > tddivstraindu = tCMLeaderTurbulence->ddivstraindu( tDofDerivative );

                    // evaluate ddivstraindu by FD
                    Matrix< DDRMat > tddivstrainduFD;
                    tCMLeaderTurbulence->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                    if( tDofDerivative( 0 ) == MSI::Dof_Type::VX || tDofDerivative( 0 ) == MSI::Dof_Type::P )
                    {
                        // evaluate the laminar part of ddivstraindu
                        Matrix< DDRMat > tddivstrainduLaminar = tCMLeaderLaminar->ddivstraindu( tDofDerivative );

                        // get the turbulent contribution
                        tddivstraindu = tddivstraindu - tddivstrainduLaminar;
                        tddivstrainduFD = tddivstrainduFD - tddivstrainduLaminar;
                    }

                    // check that analytical and FD match
                    bool tCheckDivStrainTurbulence = fem::check(
                            tddivstraindu,
                            tddivstrainduFD,
                            tEpsilon );
                    REQUIRE( tCheckDivStrainTurbulence );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "CM_Fluid_Turbulence", "[CM_Fluid_Turbulence]" )
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
    moris::Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunc );
    tPropViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );
    //tPropViscosity->set_dof_type_list( { tPDofTypes } );
    //tPropViscosity->set_val_function( tPFIValFunc );
    //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropKinViscosity = std::make_shared< fem::Property >();
    tPropKinViscosity->set_val_function( tConstValFunc );
    tPropKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //tPropDensity->set_dof_type_list( { tPDofTypes } );
    //tPropDensity->set_val_function( tPFIValFunc );
    //tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMLeaderTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMLeaderTurbulence->set_property( tPropViscosity, "Viscosity" );
    tCMLeaderTurbulence->set_property( tPropKinViscosity, "KinViscosity" );
    tCMLeaderTurbulence->set_property( tPropDensity, "Density" );
    tCMLeaderTurbulence->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMLeaderTurbulence->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeaderTurbulence->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tCMLeaderTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set leader dof type map
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P  ) )        = 1;
    tCMLeaderTurbulence->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // build global dof type list
    tCMLeaderTurbulence->get_global_dof_type_list();

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( iSpaceDim, 1, 10.0 );

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
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );

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
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );

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

        // set space dimensions for property, CM and SP
        tCMLeaderTurbulence->set_space_dim( iSpaceDim );

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

            // fill random coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator viscosity
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatVis );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeaderTurbulence->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMLeaderTurbulence->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMLeaderTurbulence->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeaderTurbulence->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // populate the requested leader dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeaderTurbulence->get_global_dof_type_list();

                // populate the test leader dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeaderTurbulence->get_dof_type_list();
                //moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes = { tVelDofTypes };

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMLeaderTurbulence->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMLeaderTurbulence->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckFluxTurbulence = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckFluxTurbulence );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate dstraindu
                    Matrix< DDRMat > tdstraindu = tCMLeaderTurbulence->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeaderTurbulence->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckStrainTurbulence = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckStrainTurbulence );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMLeaderTurbulence->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMLeaderTurbulence->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                    // check that analytical and FD match
                    bool tCheckTractionTurbulence = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckTractionTurbulence );

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for( uint iTestDof = 0; iTestDof < tLeaderDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( iTestDof );

                        // evaluate dtesttractiondu
                        Matrix< DDRMat > tdtesttractiondu = tCMLeaderTurbulence->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tJump,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMLeaderTurbulence->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD,
                                tPerturbation,
                                tNormal,
                                tJump );

                        // check that analytical and FD match
                        bool tCheckTestTractionTurbulence = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                        REQUIRE( tCheckTestTractionTurbulence );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMLeaderTurbulence->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMLeaderTurbulence->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivFluxTurbulence = fem::check( tddivfluxdu, tddivfluxduFD, tEpsilon );
                    REQUIRE( tCheckDivFluxTurbulence );

                    // div strain
                    //------------------------------------------------------------------------------
                    // evaluate ddivstraindu
                    Matrix< DDRMat > tddivstraindu = tCMLeaderTurbulence->ddivstraindu( tDofDerivative );

                    // evaluate ddivstraindu by FD
                    Matrix< DDRMat > tddivstrainduFD;
                    tCMLeaderTurbulence->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivStrainTurbulence = fem::check( tddivstraindu, tddivstrainduFD, tEpsilon );
                    REQUIRE( tCheckDivStrainTurbulence );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/

