/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM_Diff_Lin_Iso_Turbulence_Smagorinsky.cpp
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
#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence_Spalart_Allmaras.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "CM_Diff_Lin_Iso_Turbulence_Spallart_Allmaras", "[CM_Diff_Lin_Iso_Turbulence_Spallart_Allmaras]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-6;

    // define a perturbation relative size
    real tPerturbation = 1.0E-5;

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
    Vector< Vector< MSI::Dof_Type > > tTempDofTypes = { { MSI::Dof_Type::TEMP } };
    Vector< Vector< MSI::Dof_Type > > tVisDofTypes  = { { MSI::Dof_Type::VISCOSITY } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = { tTempDofTypes( 0 ), tVisDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropKinViscosity = std::make_shared< fem::Property >();
    tPropKinViscosity->set_parameters( { { { 0.7 / 2.0 } } } );
    tPropKinViscosity->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropTurbPrandtl = std::make_shared< fem::Property >();
    tPropTurbPrandtl->set_parameters( { { { 0.7 } } } );
    tPropTurbPrandtl->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { { { 2.0 } } } );
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
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropDensity, "Density" );
    tCMLeader->set_property( tPropConductivity, "Conductivity" );
    tCMLeader->set_property( tPropHeatCapacity, "HeatCapacity" );
    tCMLeader->set_property( tPropTurbPrandtl, "TurbulentPrandtl" );
    tCMLeader->set_property( tPropKinViscosity, "KinematicViscosity" );

    tCMLeader->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // set size and populate the set leader dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

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
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatTemp );

            // create the field interpolator velocity
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVis );

            // create a field interpolator manager
            Vector< Vector< gen::PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeader->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMLeader->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMLeader->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // populate the requested leader dof type for CM
                Vector< Vector< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeader->get_global_dof_type_list();

                // populate the test leader dof type for CM
                Vector< Vector< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    const Vector< MSI::Dof_Type >& tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );
                    // std::cout << "tDofDerivative " << static_cast< uint >( tDofDerivative( 0 ) ) << std::endl;

                    // get turbulence specific CM
                    CM_Diffusion_Linear_Isotropic_Turbulence_Spalart_Allmaras* tCMLeaderPtr =
                            dynamic_cast< CM_Diffusion_Linear_Isotropic_Turbulence_Spalart_Allmaras* >( tCMLeader.get() );

                    // chi
                    //------------------------------------------------------------------------------
                    // evaluate analytical derivative
                    Matrix< DDRMat > tdchidu = tCMLeaderPtr->dchidu( tDofDerivative );

                    // evaluate derivative by FD
                    Matrix< DDRMat > tdchiduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::CHI,
                            tdchiduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check tdchidu analytical and FD match " << std::endl;
                    // print( tdchidu, "tdchidu" );
                    // print( tdchiduFD, "tdchiduFD" );
                    bool tCheckChi = fem::check( tdchidu, tdchiduFD, tEpsilon );
                    REQUIRE( tCheckChi );

                    // chi space derivative
                    //------------------------------------------------------------------------------
                    // evaluate analytical derivative
                    Matrix< DDRMat > tdchidxdu = tCMLeaderPtr->dchidxdu( tDofDerivative, 1 );

                    // evaluate derivative by FD
                    Matrix< DDRMat > tdchidxduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::CHI_SPACE_DER,
                            tdchidxduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check tdchidxdu analytical and FD match " << std::endl;
                    // print( tdchidxdu, "tdchidxdu" );
                    // print( tdchidxduFD, "tdchidxduFD" );
                    bool tCheckChiSpaceDer = fem::check( tdchidxdu, tdchidxduFD, tEpsilon );
                    REQUIRE( tCheckChiSpaceDer );

                    // fv1
                    //------------------------------------------------------------------------------
                    // evaluate analytical derivative
                    Matrix< DDRMat > tdfv1du = tCMLeaderPtr->dfv1du( tDofDerivative );

                    // evaluate derivative by FD
                    Matrix< DDRMat > tdfv1duFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::FV1,
                            tdfv1duFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check tdfv1du analytical and FD match " << std::endl;
                    // print( tdfv1du, "tdfv1du" );
                    // print( tdfv1duFD, "tdfv1duFD" );
                    bool tCheckFv1 = fem::check( tdfv1du, tdfv1duFD, tEpsilon );
                    REQUIRE( tCheckFv1 );

                    // fv1 space derivative
                    //------------------------------------------------------------------------------
                    // evaluate analytical derivative
                    Matrix< DDRMat > tdfv1dxdu = tCMLeaderPtr->dfv1dxdu( tDofDerivative, 1 );

                    // evaluate derivative by FD
                    Matrix< DDRMat > tdfv1dxduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::FV1_SPACE_DER,
                            tdfv1dxduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check tdfv1dxdu analytical and FD match " << std::endl;
                    // print( tdfv1dxdu, "tdfv1dxdu" );
                    // print( tdfv1dxduFD, "tdfv1dxduFD" );
                    bool tCheckFv1SpaceDer = fem::check( tdfv1dxdu, tdfv1dxduFD, tEpsilon );
                    REQUIRE( tCheckFv1SpaceDer );

                    // turbulent dynamic viscosity
                    //------------------------------------------------------------------------------
                    // evaluate dturbdynviscdu
                    Matrix< DDRMat > tdturbdynviscdu = tCMLeaderPtr->dturbdynviscdu( tDofDerivative );

                    // evaluate dturbdynviscdu by FD
                    Matrix< DDRMat > tdturbdynviscduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TURB_DYN_VISC,
                            tdturbdynviscduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check tdturbdynviscdu analytical and FD match " << std::endl;
                    // print( tdturbdynviscdu, "tdturbdynviscdu" );
                    // print( tdturbdynviscduFD, "tdturbdynviscduFD" );
                    bool tCheckTurbDynVisc = fem::check( tdturbdynviscdu, tdturbdynviscduFD, tEpsilon );
                    REQUIRE( tCheckTurbDynVisc );

                    // turbulent dynamic viscosity space der
                    //------------------------------------------------------------------------------
                    // evaluate dturbdynviscdxdu
                    Matrix< DDRMat > tdturbdynviscdxdu = tCMLeaderPtr->dturbdynviscdxdu( tDofDerivative, 1 );

                    // evaluate dturbdynviscdxdu by FD
                    Matrix< DDRMat > tdturbdynviscdxduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TURB_DYN_VISC_SPACE_DER,
                            tdturbdynviscdxduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check tdturbdynviscdxdu analytical and FD match " << std::endl;
                    // print( tdturbdynviscdxdu, "tdturbdynviscdxdu" );
                    // print( tdturbdynviscdxduFD, "tdturbdynviscdxduFD" );
                    bool tCheckTurbDynViscSpaceDer = fem::check( tdturbdynviscdxdu, tdturbdynviscdxduFD, tEpsilon );
                    REQUIRE( tCheckTurbDynViscSpaceDer );

                    // test turbulent dynamic viscosity
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    Vector< Vector< MSI::Dof_Type > > tTestDofTypes = { tTempDofTypes( 0 ) };
                    for ( uint iTestDof = 0; iTestDof < tTestDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        const Vector< MSI::Dof_Type >& tDofTest = tTestDofTypes( iTestDof );
                        // std::cout << "tDofTest " << static_cast< uint >( tDofTest( 0 ) ) << std::endl;

                        // evaluate dtestturbdynviscdu
                        Matrix< DDRMat > tdtestturbdynviscdu = tCMLeaderPtr->dtestturbdynviscdu(
                                tDofDerivative,
                                tDofTest );

                        // evaluate dtestturbdynviscduFD by FD
                        Matrix< DDRMat > tdtestturbdynviscduFD;
                        tCMLeader->eval_derivative_FD(
                                CM_Request_Type::TEST_TURB_DYN_VISC,
                                tdtestturbdynviscduFD,
                                tDofDerivative,
                                tPerturbation,
                                tDofTest,
                                tNormal,
                                tJump );

                        // check that analytical and FD match
                        // std::cout << "check tdtestturbdynviscdu analytical and FD match " << std::endl;
                        // print( tdtestturbdynviscdu, "tdtestturbdynviscdu" );
                        // print( tdtestturbdynviscduFD, "tdtestturbdynviscduFD" );
                        bool tCheckTestTurbDynVisc = fem::check( tdtestturbdynviscdu, tdtestturbdynviscduFD, tEpsilon );
                        REQUIRE( tCheckTestTurbDynVisc );
                    }
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/
