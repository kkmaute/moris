/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM_Struc_NonLinear_Neo_Hookean.cpp
 *
 */

#include "catch.hpp"

#define protected public
#define private public
// FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private

// LINALG/src
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
// FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "fn_FEM_Check.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"

using namespace moris;
using namespace fem;

TEST_CASE( "CM_Struc_NonLinear_Neo_Hookean",
        "[CM_Struc_NonLinear_Neo_Hookean]" )
{
    // define an epsilon environment
    real tEpsilon = 5.0E-4;

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
    Vector< Vector< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = tDispDofTypes;

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );
    tPropEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );
    tPropNu->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN );
    tCMLeaderStrucLinIso->set_dof_type_list( { tDispDofTypes } );
    tCMLeaderStrucLinIso->set_property( tPropEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tCMLeaderStrucLinIso->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeaderStrucLinIso->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeaderStrucLinIso->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderStrucLinIso->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tCMLeaderStrucLinIso->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderStrucLinIso->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // build global dof type list
    tCMLeaderStrucLinIso->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( iSpaceDim, 1, 10.0 );

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

                // set velocity dof types
                tDispDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };

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

                // set velocity dof types
                tDispDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } };

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
        Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimensions for property, CM and SP
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );

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

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule(
                    tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat_Elast( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );

            tLeaderDOFHatVel = tLeaderDOFHatVel * 0.01;

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( 0.1 * tLeaderDOFHatVel );

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >        tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeaderStrucLinIso->mSet->mLeaderFIManager = &tFIManager;

            // set CM field interpolator manager
            tCMLeaderStrucLinIso->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset CM evaluation flags
                tCMLeaderStrucLinIso->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeaderStrucLinIso->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // populate the requested leader dof type for CM
                Vector< Vector< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeaderStrucLinIso->get_global_dof_type_list();

                // populate the test leader dof type for CM
                Vector< Vector< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeaderStrucLinIso->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Vector< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate Lagrange Green strain
                    //                    const Matrix< DDRMat > & tLGStrain =
                    //                            tCMLeaderStrucLinIso->strain( CM_Function_Type::LAGRANGIAN );

                    //                    const Matrix< DDRMat > & tLGTestStrain =
                    //                            tCMLeaderStrucLinIso->testStrain( CM_Function_Type::LAGRANGIAN );
                    //                    print(tLGTestStrain,"tLGTestStrain");

                    // evaluate dLGStraindu
                    Matrix< DDRMat > tdLGStraindu =
                            tCMLeaderStrucLinIso->dStraindDOF(
                                    tDofDerivative,
                                    CM_Function_Type::LAGRANGIAN );

                    // evaluate dLGStraindu by FD
                    Matrix< DDRMat > tdLGStrainduFD;
                    tCMLeaderStrucLinIso->eval_dStraindDOF_FD(
                            tDofDerivative,
                            tdLGStrainduFD,
                            tPerturbation,
                            fem::FDScheme_Type::POINT_5,
                            CM_Function_Type::LAGRANGIAN );

                    // check that analytical and FD match
                    bool tCheckLGStrain = fem::check( tdLGStraindu, tdLGStrainduFD, tEpsilon );
                    REQUIRE( tCheckLGStrain );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate PK1
                    //                    const Matrix< DDRMat > & tPK1Stress =
                    //                            tCMLeaderStrucLinIso->flux( CM_Function_Type::PK1 );
                    //                    print( tPK1Stress, "tPK1Stress" );

                    // evaluate PK2
                    //                    const Matrix< DDRMat > & tPK2Stress =
                    //                            tCMLeaderStrucLinIso->flux( CM_Function_Type::PK2 );
                    //                    print( tPK2Stress, "tPK2Stress" );

                    // evaluate Cauchy
                    //                    const Matrix< DDRMat > & tCauchyStress =
                    //                            tCMLeaderStrucLinIso->flux( CM_Function_Type::CAUCHY );
                    //                    print( tCauchyStress, "tCauchyStress" );

                    // evaluate dPK2du
                    Matrix< DDRMat > tdPK2Stressdu =
                            tCMLeaderStrucLinIso->dFluxdDOF(
                                    tDofDerivative,
                                    CM_Function_Type::PK2 );
                    //                    print(tdPK2Stressdu,"tdPK2Stressdu");

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdPK2StressduFD;
                    tCMLeaderStrucLinIso->eval_dFluxdDOF_FD(
                            tDofDerivative,
                            tdPK2StressduFD,
                            tPerturbation,
                            fem::FDScheme_Type::POINT_5,
                            CM_Function_Type::PK2 );
                    //                    print(tdPK2StressduFD,"tdPK2StressduFD");

                    // check that analytical and FD match
                    bool tCheckPK2 = fem::check( tdPK2Stressdu, tdPK2StressduFD, tEpsilon );
                    REQUIRE( tCheckPK2 );

                    // evaluate dPK1du
                    Matrix< DDRMat > tdPK1Stressdu =
                            tCMLeaderStrucLinIso->dFluxdDOF(
                                    tDofDerivative,
                                    CM_Function_Type::PK1 );
                    //                    print(tdPK1Stressdu,"tdPK1Stressdu");

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdPK1StressduFD;
                    tCMLeaderStrucLinIso->eval_dFluxdDOF_FD(
                            tDofDerivative,
                            tdPK1StressduFD,
                            tPerturbation,
                            fem::FDScheme_Type::POINT_5,
                            CM_Function_Type::PK1 );
                    //                    print(tdPK1StressduFD,"tdPK1StressduFD");

                    // check that analytical and FD match
                    bool tCheckPK1 = fem::check( tdPK1Stressdu, tdPK1StressduFD, tEpsilon );
                    REQUIRE( tCheckPK1 );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate PK1 traction
                    //                    const Matrix< DDRMat > & tPK1Traction =
                    //                            tCMLeaderStrucLinIso->traction( tNormal, CM_Function_Type::PK1 );
                    //                    print( tPK1Traction, "tPK1Traction" );

                    // evaluate PK2 traction
                    //                    const Matrix< DDRMat > & tPK2Traction =
                    //                            tCMLeaderStrucLinIso->traction( tNormal, CM_Function_Type::PK2 );
                    //                    print( tPK2Traction, "tPK2Traction" );

                    //                    // evaluate Cauchy
                    //                    const Matrix< DDRMat > & tCauchyTraction =
                    //                            tCMLeaderStrucLinIso->traction( tNormal, CM_Function_Type::CAUCHY );
                    ////                    print( tCauchyTraction, "tCauchyTraction" );

                    //                    // evaluate dPK2Tractiondu
                    Matrix< DDRMat > tdPK2Tractiondu =
                            tCMLeaderStrucLinIso->dTractiondDOF( tDofDerivative, tNormal, CM_Function_Type::PK2 );

                    // evaluate dPK2Tractiondu by FD
                    Matrix< DDRMat > tdPK2TractionduFD;
                    tCMLeaderStrucLinIso->eval_dtractiondu_FD(
                            tDofDerivative,
                            tdPK2TractionduFD,
                            tPerturbation,
                            tNormal,
                            fem::FDScheme_Type::POINT_5,
                            CM_Function_Type::PK2 );

                    // check that analytical and FD match
                    bool tCheckPK2Traction = fem::check( tdPK2Tractiondu, tdPK2TractionduFD, tEpsilon );
                    REQUIRE( tCheckPK2Traction );

                    // evaluate dPK1Tractiondu
                    Matrix< DDRMat > tdPK1Tractiondu =
                            tCMLeaderStrucLinIso->dTractiondDOF( tDofDerivative, tNormal, CM_Function_Type::PK1 );

                    // evaluate dPK1Tractiondu by FD
                    Matrix< DDRMat > tdPK1TractionduFD;
                    tCMLeaderStrucLinIso->eval_dtractiondu_FD(
                            tDofDerivative,
                            tdPK1TractionduFD,
                            tPerturbation,
                            tNormal,
                            fem::FDScheme_Type::POINT_5,
                            CM_Function_Type::PK1 );

                    // check that analytical and FD match
                    bool tCheckPK1Traction = fem::check( tdPK1Tractiondu, tdPK1TractionduFD, tEpsilon );
                    REQUIRE( tCheckPK1Traction );

                    //                    // test traction
                    //                    //------------------------------------------------------------------------------
                    //
                    //                    // loop over test dof type
                    for ( uint iTestDof = 0; iTestDof < tLeaderDofTypes.size(); iTestDof++ )
                    {

                        //                        const Matrix< DDRMat > & tLGdTestStraindDOF =
                        //                                tCMLeaderStrucLinIso->dTestStraindDOF( iTestDof, CM_Function_Type::LAGRANGIAN );

                        //                        // get the test dof type
                        Vector< MSI::Dof_Type > tDofTest = tLeaderDofTypes( iTestDof );
                        //
                        //                        // evaluate PK1 test traction
                        //                        const Matrix< DDRMat >& tPK1TestTraction =
                        //                                tCMLeaderStrucLinIso->testTraction( tNormal, tDofTest, CM_Function_Type::PK1 );
                        //                        print( tPK1TestTraction, "tPK1TestTraction" );
                        //
                        //                        // evaluate PK2 test traction
                        //                        const Matrix< DDRMat >& tPK2TestTraction =
                        //                                tCMLeaderStrucLinIso->testTraction( tNormal, tDofTest, CM_Function_Type::PK2 );
                        //                        print( tPK2TestTraction, "tPK2TestTraction" );
                        //
                        ////                        // evaluate Cauchy test traction
                        ////                        const Matrix< DDRMat > & tCauchyTestTraction =
                        ////                                tCMLeaderStrucLinIso->testTraction( tDofTest, tNormal, CM_Function_Type::CAUCHY );
                        ////                        print( tCauchyTestTraction, "tCauchyTestTraction" );
                        //
                        // evaluate dPK2testtractiondu
                        Matrix< DDRMat > tdPK2testtractiondu =
                                tCMLeaderStrucLinIso->dTestTractiondDOF(
                                        tDofDerivative,
                                        tNormal,
                                        tJump,
                                        tDofTest,
                                        CM_Function_Type::PK2 );
                        //
                        //                        // evaluate dPK2testtractiondu by FD
                        Matrix< DDRMat > dPK2testtractionduFD;
                        tCMLeaderStrucLinIso->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                dPK2testtractionduFD,
                                tPerturbation,
                                tNormal,
                                tJump,
                                fem::FDScheme_Type::POINT_5,
                                CM_Function_Type::PK2 );

                        //                        print(tdPK2testtractiondu-dPK2testtractionduFD,"dif");
                        //
                        //                        // check that analytical and FD match
                        bool tCheckPK2TestTraction = fem::check( tdPK2testtractiondu, dPK2testtractionduFD, tEpsilon );
                        REQUIRE( tCheckPK2TestTraction );
                        //
                        //                        // evaluate dPK1testtractiondu
                        //                        Matrix< DDRMat > tdPK1testtractiondu =
                        //                                tCMLeaderStrucLinIso->dTestTractiondDOF(
                        //                                tDofDerivative,
                        //                                tNormal,
                        //                                tJump,
                        //                                tDofTest,
                        //                                CM_Function_Type::PK1 );
                        //                        print(tdPK1testtractiondu,"tdPK1testtractiondu");
                        //
                        //                        // evaluate dPK1testtractiondu by FD
                        Matrix< DDRMat > dPK1testtractionduFD;
                        tCMLeaderStrucLinIso->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                dPK1testtractionduFD,
                                tPerturbation,
                                tNormal,
                                tJump,
                                fem::FDScheme_Type::POINT_5,
                                CM_Function_Type::PK1 );
                        //                        print(dPK1testtractionduFD,"dPK1testtractionduFD");

                        //                        print(tdPK1testtractiondu-dPK1testtractionduFD,"dif");
                        //
                        //                        // check that analytical and FD match
                        //                        bool tCheckPK1TestTraction = fem::check( tdPK1testtractiondu, dPK1testtractionduFD, tEpsilon );
                        //                        REQUIRE( tCheckPK1TestTraction );
                    }
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/
