/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM_Struc_Linear_Isotropic_Damage.cpp
 *
 */

#include "catch.hpp"

#define protected public
#define private   public
// FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp"
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
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"

using namespace moris;
using namespace fem;

TEST_CASE( "CM_Struc_Linear_Isotropic_Damage_020_Growth", "[CM_Struc_Lin_Iso_Dam_020_Growth]" )
{
    // Damage setup - 0,2,0
    // Local equivalent strain - 0 - energy release rate
    // Damage law - 2 - smooth exponential
    // Smoothing law - 0 - no smoothing
    moris::Cell< Matrix< DDRMat > > tDamageParameters = { //
        { { 0.0 } },                                      //
        { { 2.0, 1.0e-3, 10.0 } },                        //
        { { 0.0 } }
    };

    // define an epsilon environment
    real tEpsilon = 1.0E-5;

    // define a perturbation relative size
    real tPerturbation = 1.0E-5;

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
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };
    Matrix< DDRMat > tNumHalfCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes       = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes           = { tDispDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropNu, "PoissonRatio" );
    tCMLeader->set_local_properties();
    tCMLeader->set_parameters( tDamageParameters );

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // set size and populate the set master dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // std::cout << "Space dimension " << iSpaceDim << std::endl;

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

        // set space dimensions for property, CM, and SP
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // std::cout << "Interpolation order " << iInterpOrder << std::endl;

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

            // create a cell of field interpolators current
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // get number of coefficients
            uint tNumHalfCoeffCurrent = tNumHalfCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
            uint tNumCoeffCurrent     = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // create the field interpolator for displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 1.0e-3 * tLeaderDOFHatNlEqStrain;
            tLeaderDOFHatNlEqStrain( { tNumHalfCoeffCurrent, tNumCoeffCurrent - 1 }, { 0, 0 } ) =
                    2.0 * tLeaderDOFHatNlEqStrain( { 0, tNumHalfCoeffCurrent - 1 }, { 0, 0 } );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            // fill coefficients for master FI
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeader->mSet->mLeaderFIManager = &tFIManager;

            // set CM field interpolator manager
            tCMLeader->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // std::cout << "Integration point " << iGP << " / " << tNumGPs << std::endl;

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // reset CM evaluation flags
                tCMLeader->reset_eval_flags();

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeader->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // std::cout << "Requested dof derivative point " << iGP << " / " << tNumGPs << std::endl;

                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // cast constitutive model base class pointer to elasticity damage constitutive model
                    CM_Struc_Linear_Isotropic_Damage* tCMLeaderPtr =
                            dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMLeader.get() );

                    // equivalent strain
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > teqstrain = tCMLeaderPtr->equivalent_strain();
                    // print( teqstrain, "teqstrain" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdeqstraindu = tCMLeaderPtr->dEqStraindu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdeqstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::EQSTRAIN,
                            tdeqstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump,
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdeqstraindu analytical and FD match " << std::endl;
                    // print( tdeqstraindu, "tdeqstraindu" );
                    // print( tdeqstrainduFD, "tdeqstrainduFD" );
                    bool tCheckdEqStraindu = fem::check( tdeqstraindu, tdeqstrainduFD, tEpsilon );
                    REQUIRE( tCheckdEqStraindu );

                    // nonlocal equivalent strain history
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > thistory = tCMLeaderPtr->history();
                    // print( thistory, "thistory" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdhistorydu = tCMLeaderPtr->dHistorydu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdhistoryduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::HISTORY,
                            tdhistoryduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdhistorydu analytical and FD match " << std::endl;
                    // print( tdhistorydu, "tdhistorydu" );
                    // print( tdhistoryduFD, "tdhistoryduFD" );
                    bool tCheckdHistorydu = fem::check( tdhistorydu, tdhistoryduFD, tEpsilon );
                    REQUIRE( tCheckdHistorydu );

                    // damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tdamage = tCMLeaderPtr->damage();
                    // print( tdamage, "tdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tddamagedu = tCMLeaderPtr->dDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tddamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::DAMAGE,
                            tddamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tddamagedu, "tddamagedu" );
                    // print( tddamageduFD, "tddamageduFD" );
                    bool tCheckdDamagedu = fem::check( tddamagedu, tddamageduFD, tEpsilon );
                    REQUIRE( tCheckdDamagedu );

                    // smooth damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tsmoothdamage = tCMLeaderPtr->smooth_damage();
                    // print( tsmoothdamage, "tsmoothdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tdsmoothdamagedu = tCMLeaderPtr->dSmoothDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tdsmoothdamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::SMOOTH_DAMAGE,
                            tdsmoothdamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tdsmoothdamagedu, "tdsmoothdamagedu" );
                    // print( tdsmoothdamageduFD, "tdsmoothdamageduFD" );
                    bool tCheckdSmoothDamagedu = fem::check( tdsmoothdamagedu, tdsmoothdamageduFD, tEpsilon );
                    REQUIRE( tCheckdSmoothDamagedu );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate flux
                    Matrix< DDRMat > tflux = tCMLeader->flux();
                    // print( tflux, "tflux" );

                    // evaluate derivative of flux
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
                    // std::cout << "check tdfluxdu analytical and FD match " << std::endl;
                    // print( tdfluxdu, "tdfluxdu" );
                    // print( tdfluxduFD, "tdfluxduFD" );
                    bool tCheckdFluxdu = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckdFluxdu );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate traction
                    Matrix< DDRMat > ttraction = tCMLeader->traction( tNormal );
                    // print( ttraction, "ttraction" );

                    // evaluate derivative of traction
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
                    // std::cout << "check tdtractiondu analytical and FD match " << std::endl;
                    // print( tdtractiondu, "tdtractiondu" );
                    // print( tdtractionduFD, "tdtractionduFD" );
                    bool tCheckdTractiondu = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckdTractiondu );

                    // test traction -- only displacement as test dof type !!!
                    //------------------------------------------------------------------------------
                    // get the test dof type
                    Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( 0 );

                    // evaluate test traction
                    Matrix< DDRMat > ttesttraction = tCMLeader->testTraction(
                            tNormal,
                            tDofTest );
                    // print( ttesttraction, "ttesttraction" );

                    // evaluate derivative test traction
                    Matrix< DDRMat > tdtesttractiondu = tCMLeader->dTestTractiondDOF(
                            tDofDerivative,
                            tNormal,
                            tJump,
                            tDispDofTypes( 0 ) );

                    // evaluate dtesttractiondu by FD
                    Matrix< DDRMat > tdtesttractionduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TEST_TRACTION,
                            tdtesttractionduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDispDofTypes( 0 ),
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdtesttractiondu analytical and FD match " << std::endl;
                    // print( tdtesttractiondu, "tdtesttractiondu" );
                    // print( tdtesttractionduFD, "tdtesttractionduFD" );
                    bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                    REQUIRE( tCheckTestTraction );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate strain
                    Matrix< DDRMat > tstrain = tCMLeader->strain();
                    // print( tstrain, "tstrain" );

                    // evaluate derivative of strain
                    Matrix< DDRMat > tdstraindu = tCMLeader->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdstraindu analytical and FD match " << std::endl;
                    // print( tdstraindu, "tdstraindu" );
                    // print( tdstrainduFD, "tdstrainduFD" );
                    bool tCheckdStraindu = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckdStraindu );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "CM_Struc_Linear_Isotropic_Damage_020_No_Growth", "[CM_Struc_Lin_Iso_Dam_020_No_Growth]" )
{
    // Damage setup - 0,2,0
    // Local equivalent strain - 0 - energy release rate
    // Damage law - 2 - smooth exponential
    // Smoothing law - 0 - no smoothing
    moris::Cell< Matrix< DDRMat > > tDamageParameters = { //
        { { 0.0 } },                                      //
        { { 2.0, 1.0e-3, 10.0 } },                        //
        { { 0.0 } }
    };

    // define an epsilon environment
    real tEpsilon = 1.0E-5;

    // define a perturbation relative size
    real tPerturbation = 1.0E-5;

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
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs     = { { 8, 18, 32 }, { 16, 54, 128 } };
    Matrix< DDRMat > tNumHalfCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes       = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes           = { tDispDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropNu, "PoissonRatio" );
    tCMLeader->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // set size and populate the set master dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // std::cout << "Space dimension " << iSpaceDim << std::endl;

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

        // set space dimensions for property, CM, and SP
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeader->set_parameters( tDamageParameters );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // std::cout << "Interpolation order " << iInterpOrder << std::endl;

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

            // create a cell of field interpolators current
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // get number of coefficients
            uint tNumHalfCoeffCurrent = tNumHalfCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
            uint tNumCoeffCurrent     = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // create the field interpolator for displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 1.0e-3 * tLeaderDOFHatNlEqStrain;
            tLeaderDOFHatNlEqStrain( { tNumHalfCoeffCurrent, tNumCoeffCurrent - 1 }, { 0, 0 } ) =
                    0.5 * tLeaderDOFHatNlEqStrain( { 0, tNumHalfCoeffCurrent - 1 }, { 0, 0 } );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            // fill coefficients for master FI
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeader->mSet->mLeaderFIManager = &tFIManager;

            // set CM field interpolator manager
            tCMLeader->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // std::cout << "Integration point " << iGP << " / " << tNumGPs << std::endl;

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // reset CM evaluation flags
                tCMLeader->reset_eval_flags();

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeader->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // cast constitutive model base class pointer to elasticity damage constitutive model
                    CM_Struc_Linear_Isotropic_Damage* tCMLeaderPtr =
                            dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMLeader.get() );

                    // equivalent strain
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > teqstrain = tCMLeaderPtr->equivalent_strain();
                    // print( teqstrain, "teqstrain" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdeqstraindu = tCMLeaderPtr->dEqStraindu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdeqstrainduFD;
                    tCMLeader->eval_derivative_FD( CM_Request_Type::EQSTRAIN, tdeqstrainduFD, tDofDerivative, tPerturbation, tDofDerivative, tNormal, tJump, FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdeqstraindu analytical and FD match " << std::endl;
                    // print( tdeqstraindu, "tdeqstraindu" );
                    // print( tdeqstrainduFD, "tdeqstrainduFD" );
                    bool tCheckdEqStraindu = fem::check( tdeqstraindu, tdeqstrainduFD, tEpsilon );
                    REQUIRE( tCheckdEqStraindu );

                    // nonlocal equivalent strain history
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > thistory = tCMLeaderPtr->history();
                    // print( thistory, "thistory" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdhistorydu = tCMLeaderPtr->dHistorydu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdhistoryduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::HISTORY,
                            tdhistoryduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdhistorydu analytical and FD match " << std::endl;
                    // print( tdhistorydu, "tdhistorydu" );
                    // print( tdhistoryduFD, "tdhistoryduFD" );
                    bool tCheckdHistorydu = fem::check( tdhistorydu, tdhistoryduFD, tEpsilon );
                    REQUIRE( tCheckdHistorydu );

                    // damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tdamage = tCMLeaderPtr->damage();
                    // print( tdamage, "tdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tddamagedu = tCMLeaderPtr->dDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tddamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::DAMAGE,
                            tddamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tddamagedu, "tddamagedu" );
                    // print( tddamageduFD, "tddamageduFD" );
                    bool tCheckdDamagedu = fem::check( tddamagedu, tddamageduFD, tEpsilon );
                    REQUIRE( tCheckdDamagedu );

                    // smooth damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tsmoothdamage = tCMLeaderPtr->smooth_damage();
                    // print( tsmoothdamage, "tsmoothdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tdsmoothdamagedu = tCMLeaderPtr->dSmoothDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tdsmoothdamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::SMOOTH_DAMAGE,
                            tdsmoothdamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tdsmoothdamagedu, "tdsmoothdamagedu" );
                    // print( tdsmoothdamageduFD, "tdsmoothdamageduFD" );
                    bool tCheckdSmoothDamagedu = fem::check( tdsmoothdamagedu, tdsmoothdamageduFD, tEpsilon );
                    REQUIRE( tCheckdSmoothDamagedu );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate flux
                    Matrix< DDRMat > tflux = tCMLeader->flux();
                    // print( tflux, "tflux" );

                    // evaluate derivative of flux
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
                    // std::cout << "check tdfluxdu analytical and FD match " << std::endl;
                    // print( tdfluxdu, "tdfluxdu" );
                    // print( tdfluxduFD, "tdfluxduFD" );
                    bool tCheckdFluxdu = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckdFluxdu );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate traction
                    Matrix< DDRMat > ttraction = tCMLeader->traction( tNormal );
                    // print( ttraction, "ttraction" );

                    // evaluate derivative of traction
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
                    // std::cout << "check tdtractiondu analytical and FD match " << std::endl;
                    // print( tdtractiondu, "tdtractiondu" );
                    // print( tdtractionduFD, "tdtractionduFD" );
                    bool tCheckdTractiondu = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckdTractiondu );

                    // test traction -- only displacement as test dof type !!!
                    //------------------------------------------------------------------------------
                    // get the test dof type
                    Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( 0 );

                    // evaluate test traction
                    Matrix< DDRMat > ttesttraction = tCMLeader->testTraction(
                            tNormal,
                            tDofTest );
                    // print( ttesttraction, "ttesttraction" );

                    // evaluate derivative test traction
                    Matrix< DDRMat > tdtesttractiondu = tCMLeader->dTestTractiondDOF(
                            tDofDerivative,
                            tNormal,
                            tJump,
                            tDispDofTypes( 0 ) );

                    // evaluate dtesttractiondu by FD
                    Matrix< DDRMat > tdtesttractionduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TEST_TRACTION,
                            tdtesttractionduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDispDofTypes( 0 ),
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdtesttractiondu analytical and FD match " << std::endl;
                    // print( tdtesttractiondu, "tdtesttractiondu" );
                    // print( tdtesttractionduFD, "tdtesttractionduFD" );
                    bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                    REQUIRE( tCheckTestTraction );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate strain
                    Matrix< DDRMat > tstrain = tCMLeader->strain();
                    // print( tstrain, "tstrain" );

                    // evaluate derivative of strain
                    Matrix< DDRMat > tdstraindu = tCMLeader->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdstraindu analytical and FD match " << std::endl;
                    // print( tdstraindu, "tdstraindu" );
                    // print( tdstrainduFD, "tdstrainduFD" );
                    bool tCheckdStraindu = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckdStraindu );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "CM_Struc_Linear_Isotropic_Damage_101_Growth", "[CM_Struc_Lin_Iso_Dam_101_Growth]" )
{
    // Damage setup - 1,0,1
    // Local equivalent strain - 1 - tensile/compressive strength
    // Damage law - 0- linear
    // Smoothing law - 1 - ks smoothing
    moris::Cell< Matrix< DDRMat > > tDamageParameters = { //
        { { 1.0, 4.0 } },                                 //
        { { 0.0, 1.0e-3, 10.0 } },                        //
        { { 1.0, 7.0 } }
    };

    // define an epsilon environment
    real tEpsilon = 4.0E-5;

    // define a perturbation relative size
    real tPerturbation = 1.0E-5;

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
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs     = { { 8, 18, 32 }, { 16, 54, 128 } };
    Matrix< DDRMat > tNumHalfCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes       = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes           = { tDispDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropNu, "PoissonRatio" );
    tCMLeader->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // set size and populate the set master dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // std::cout << "Space dimension " << iSpaceDim << std::endl;

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

        // set space dimensions for property, CM, and SP
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeader->set_parameters( tDamageParameters );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // std::cout << "Interpolation order " << iInterpOrder << std::endl;

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

            // create a cell of field interpolators current
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // get number of coefficients
            uint tNumHalfCoeffCurrent = tNumHalfCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
            uint tNumCoeffCurrent     = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // create the field interpolator for displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 1.0e-3 * tLeaderDOFHatNlEqStrain;
            tLeaderDOFHatNlEqStrain( { tNumHalfCoeffCurrent, tNumCoeffCurrent - 1 }, { 0, 0 } ) =
                    2.0 * tLeaderDOFHatNlEqStrain( { 0, tNumHalfCoeffCurrent - 1 }, { 0, 0 } );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            // fill coefficients for master FI
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeader->mSet->mLeaderFIManager = &tFIManager;

            // set CM field interpolator manager
            tCMLeader->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // std::cout << "Integration point " << iGP << " / " << tNumGPs << std::endl;

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // reset CM evaluation flags
                tCMLeader->reset_eval_flags();

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeader->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // std::cout << "Requested dof derivative point " << iRequestedDof << " / " << tRequestedLeaderGlobalDofTypes.size() << std::endl;

                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // cast constitutive model base class pointer to elasticity damage constitutive model
                    CM_Struc_Linear_Isotropic_Damage* tCMLeaderPtr =
                            dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMLeader.get() );

                    // equivalent strain
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > teqstrain = tCMLeaderPtr->equivalent_strain();
                    // print( teqstrain, "teqstrain" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdeqstraindu = tCMLeaderPtr->dEqStraindu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdeqstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::EQSTRAIN,
                            tdeqstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump,
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdeqstraindu analytical and FD match " << std::endl;
                    // print( tdeqstraindu, "tdeqstraindu" );
                    // print( tdeqstrainduFD, "tdeqstrainduFD" );
                    bool tCheckdEqStraindu = fem::check( tdeqstraindu, tdeqstrainduFD, tEpsilon );
                    REQUIRE( tCheckdEqStraindu );

                    // nonlocal equivalent strain history
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > thistory = tCMLeaderPtr->history();
                    // print( thistory, "thistory" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdhistorydu = tCMLeaderPtr->dHistorydu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdhistoryduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::HISTORY,
                            tdhistoryduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdhistorydu analytical and FD match " << std::endl;
                    // print( tdhistorydu, "tdhistorydu" );
                    // print( tdhistoryduFD, "tdhistoryduFD" );
                    bool tCheckdHistorydu = fem::check( tdhistorydu, tdhistoryduFD, tEpsilon );
                    REQUIRE( tCheckdHistorydu );

                    // damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tdamage = tCMLeaderPtr->damage();
                    // print( tdamage, "tdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tddamagedu = tCMLeaderPtr->dDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tddamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::DAMAGE,
                            tddamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tddamagedu, "tddamagedu" );
                    // print( tddamageduFD, "tddamageduFD" );
                    bool tCheckdDamagedu = fem::check( tddamagedu, tddamageduFD, tEpsilon );
                    REQUIRE( tCheckdDamagedu );

                    // smooth damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tsmoothdamage = tCMLeaderPtr->smooth_damage();
                    // print( tsmoothdamage, "tsmoothdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tdsmoothdamagedu = tCMLeaderPtr->dSmoothDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tdsmoothdamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::SMOOTH_DAMAGE,
                            tdsmoothdamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tdsmoothdamagedu, "tdsmoothdamagedu" );
                    // print( tdsmoothdamageduFD, "tdsmoothdamageduFD" );
                    bool tCheckdSmoothDamagedu = fem::check( tdsmoothdamagedu, tdsmoothdamageduFD, tEpsilon );
                    REQUIRE( tCheckdSmoothDamagedu );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate flux
                    Matrix< DDRMat > tflux = tCMLeader->flux();
                    // print( tflux, "tflux" );

                    // evaluate derivative of flux
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
                    // std::cout << "check tdfluxdu analytical and FD match " << std::endl;
                    // print( tdfluxdu, "tdfluxdu" );
                    // print( tdfluxduFD, "tdfluxduFD" );
                    bool tCheckdFluxdu = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckdFluxdu );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate traction
                    Matrix< DDRMat > ttraction = tCMLeader->traction( tNormal );
                    // print( ttraction, "ttraction" );

                    // evaluate derivative of traction
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
                    // std::cout << "check tdtractiondu analytical and FD match " << std::endl;
                    // print( tdtractiondu, "tdtractiondu" );
                    // print( tdtractionduFD, "tdtractionduFD" );
                    bool tCheckdTractiondu = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckdTractiondu );

                    // test traction -- only displacement as test dof type !!!
                    //------------------------------------------------------------------------------
                    // get the test dof type
                    Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( 0 );

                    // evaluate test traction
                    Matrix< DDRMat > ttesttraction = tCMLeader->testTraction(
                            tNormal,
                            tDofTest );
                    // print( ttesttraction, "ttesttraction" );

                    // evaluate derivative test traction
                    Matrix< DDRMat > tdtesttractiondu = tCMLeader->dTestTractiondDOF(
                            tDofDerivative,
                            tNormal,
                            tJump,
                            tDispDofTypes( 0 ) );

                    // evaluate dtesttractiondu by FD
                    Matrix< DDRMat > tdtesttractionduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TEST_TRACTION,
                            tdtesttractionduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDispDofTypes( 0 ),
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdtesttractiondu analytical and FD match " << std::endl;
                    // print( tdtesttractiondu, "tdtesttractiondu" );
                    // print( tdtesttractionduFD, "tdtesttractionduFD" );
                    bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                    REQUIRE( tCheckTestTraction );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate strain
                    Matrix< DDRMat > tstrain = tCMLeader->strain();
                    // print( tstrain, "tstrain" );

                    // evaluate derivative of strain
                    Matrix< DDRMat > tdstraindu = tCMLeader->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdstraindu analytical and FD match " << std::endl;
                    // print( tdstraindu, "tdstraindu" );
                    // print( tdstrainduFD, "tdstrainduFD" );
                    bool tCheckdStraindu = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckdStraindu );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "CM_Struc_Linear_Isotropic_Damage_101_Threshold", "[CM_Struc_Lin_Iso_Dam_101_Threshold]" )
{
    // Damage setup - 1,0,1
    // Local equivalent strain - 1 - tensile/compressive strength
    // Damage law - 0- linear
    // Smoothing law - 1 - ks smoothing
    moris::Cell< Matrix< DDRMat > > tDamageParameters = { //
        { { 1.0, 4.0 } },                                 //
        { { 0.0, 1.0, 10.0 } },                           //
        { { 1.0, 7.0 } }
    };

    // define an epsilon environment
    real tEpsilon = 2.0E-5;

    // define a perturbation relative size
    real tPerturbation = 1.0E-5;

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
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs     = { { 8, 18, 32 }, { 16, 54, 128 } };
    Matrix< DDRMat > tNumHalfCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes       = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes           = { tDispDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropNu, "PoissonRatio" );
    tCMLeader->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // set size and populate the set master dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // std::cout << "Space dimension " << iSpaceDim << std::endl;

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

        // set space dimensions for property, CM, and SP
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeader->set_parameters( tDamageParameters );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // std::cout << "Interpolation order " << iInterpOrder << std::endl;

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

            // create a cell of field interpolators current
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // get number of coefficients
            uint tNumHalfCoeffCurrent = tNumHalfCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
            uint tNumCoeffCurrent     = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // create the field interpolator for displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 1.0e-3 * tLeaderDOFHatNlEqStrain;
            tLeaderDOFHatNlEqStrain( { tNumHalfCoeffCurrent, tNumCoeffCurrent - 1 }, { 0, 0 } ) =
                    2.0 * tLeaderDOFHatNlEqStrain( { 0, tNumHalfCoeffCurrent - 1 }, { 0, 0 } );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            // fill coefficients for master FI
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeader->mSet->mLeaderFIManager = &tFIManager;

            // set CM field interpolator manager
            tCMLeader->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // std::cout << "Integration point " << iGP << " / " << tNumGPs << std::endl;

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // reset CM evaluation flags
                tCMLeader->reset_eval_flags();

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeader->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // std::cout << "Requested dof derivative point " << iRequestedDof << " / " << tRequestedLeaderGlobalDofTypes.size() << std::endl;

                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // cast constitutive model base class pointer to elasticity damage constitutive model
                    CM_Struc_Linear_Isotropic_Damage* tCMLeaderPtr =
                            dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMLeader.get() );

                    // equivalent strain
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > teqstrain = tCMLeaderPtr->equivalent_strain();
                    // print( teqstrain, "teqstrain" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdeqstraindu = tCMLeaderPtr->dEqStraindu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdeqstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::EQSTRAIN,
                            tdeqstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump,
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdeqstraindu analytical and FD match " << std::endl;
                    // print( tdeqstraindu, "tdeqstraindu" );
                    // print( tdeqstrainduFD, "tdeqstrainduFD" );
                    bool tCheckdEqStraindu = fem::check( tdeqstraindu, tdeqstrainduFD, tEpsilon );
                    REQUIRE( tCheckdEqStraindu );

                    // nonlocal equivalent strain history
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > thistory = tCMLeaderPtr->history();
                    // print( thistory, "thistory" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdhistorydu = tCMLeaderPtr->dHistorydu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdhistoryduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::HISTORY,
                            tdhistoryduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdhistorydu analytical and FD match " << std::endl;
                    // print( tdhistorydu, "tdhistorydu" );
                    // print( tdhistoryduFD, "tdhistoryduFD" );
                    bool tCheckdHistorydu = fem::check( tdhistorydu, tdhistoryduFD, tEpsilon );
                    REQUIRE( tCheckdHistorydu );

                    // damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tdamage = tCMLeaderPtr->damage();
                    // print( tdamage, "tdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tddamagedu = tCMLeaderPtr->dDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tddamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::DAMAGE,
                            tddamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tddamagedu, "tddamagedu" );
                    // print( tddamageduFD, "tddamageduFD" );
                    bool tCheckdDamagedu = fem::check( tddamagedu, tddamageduFD, tEpsilon );
                    REQUIRE( tCheckdDamagedu );

                    // smooth damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tsmoothdamage = tCMLeaderPtr->smooth_damage();
                    // print( tsmoothdamage, "tsmoothdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tdsmoothdamagedu = tCMLeaderPtr->dSmoothDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tdsmoothdamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::SMOOTH_DAMAGE,
                            tdsmoothdamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tdsmoothdamagedu, "tdsmoothdamagedu" );
                    // print( tdsmoothdamageduFD, "tdsmoothdamageduFD" );
                    bool tCheckdSmoothDamagedu = fem::check( tdsmoothdamagedu, tdsmoothdamageduFD, tEpsilon );
                    REQUIRE( tCheckdSmoothDamagedu );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate flux
                    Matrix< DDRMat > tflux = tCMLeader->flux();
                    // print( tflux, "tflux" );

                    // evaluate derivative of flux
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
                    // std::cout << "check tdfluxdu analytical and FD match " << std::endl;
                    // print( tdfluxdu, "tdfluxdu" );
                    // print( tdfluxduFD, "tdfluxduFD" );
                    bool tCheckdFluxdu = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckdFluxdu );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate traction
                    Matrix< DDRMat > ttraction = tCMLeader->traction( tNormal );
                    // print( ttraction, "ttraction" );

                    // evaluate derivative of traction
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
                    // std::cout << "check tdtractiondu analytical and FD match " << std::endl;
                    // print( tdtractiondu, "tdtractiondu" );
                    // print( tdtractionduFD, "tdtractionduFD" );
                    bool tCheckdTractiondu = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckdTractiondu );

                    // test traction -- only displacement as test dof type !!!
                    //------------------------------------------------------------------------------
                    // get the test dof type
                    Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( 0 );

                    // evaluate test traction
                    Matrix< DDRMat > ttesttraction = tCMLeader->testTraction(
                            tNormal,
                            tDofTest );
                    // print( ttesttraction, "ttesttraction" );

                    // evaluate derivative test traction
                    Matrix< DDRMat > tdtesttractiondu = tCMLeader->dTestTractiondDOF(
                            tDofDerivative,
                            tNormal,
                            tJump,
                            tDispDofTypes( 0 ) );

                    // evaluate dtesttractiondu by FD
                    Matrix< DDRMat > tdtesttractionduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TEST_TRACTION,
                            tdtesttractionduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDispDofTypes( 0 ),
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdtesttractiondu analytical and FD match " << std::endl;
                    // print( tdtesttractiondu, "tdtesttractiondu" );
                    // print( tdtesttractionduFD, "tdtesttractionduFD" );
                    bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                    REQUIRE( tCheckTestTraction );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate strain
                    Matrix< DDRMat > tstrain = tCMLeader->strain();
                    // print( tstrain, "tstrain" );

                    // evaluate derivative of strain
                    Matrix< DDRMat > tdstraindu = tCMLeader->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdstraindu analytical and FD match " << std::endl;
                    // print( tdstraindu, "tdstraindu" );
                    // print( tdstrainduFD, "tdstrainduFD" );
                    bool tCheckdStraindu = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckdStraindu );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "CM_Struc_Linear_Isotropic_Damage_112_Growth", "[CM_Struc_Lin_Iso_Dam_112_Growth]" )
{
    // Damage setup - 1,0,1
    // Local equivalent strain - 1 - tensile/compressive strength
    // Damage law - 0- linear
    // Smoothing law - 1 - ks smoothing
    moris::Cell< Matrix< DDRMat > > tDamageParameters = { //
        { { 1.0, 4.0 } },                                 //
        { { 1.0, 1.0e-3, 0.95, 100.0 } },                 //
        { { 2.0, 7.0 } }
    };

    // define an epsilon environment
    real tEpsilon = 1.0E-4;

    // define a perturbation relative size
    real tPerturbation = 1.0E-5;

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
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs     = { { 8, 18, 32 }, { 16, 54, 128 } };
    Matrix< DDRMat > tNumHalfCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes       = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes           = { tDispDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropNu, "PoissonRatio" );
    tCMLeader->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // set size and populate the set master dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // std::cout << "Space dimension " << iSpaceDim << std::endl;

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

        // set space dimensions for property, CM, and SP
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeader->set_parameters( tDamageParameters );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // std::cout << "Interpolation order " << iInterpOrder << std::endl;

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

            // create a cell of field interpolators current
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // get number of coefficients
            uint tNumHalfCoeffCurrent = tNumHalfCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
            uint tNumCoeffCurrent     = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // create the field interpolator for displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 1.0e-3 * tLeaderDOFHatNlEqStrain;
            tLeaderDOFHatNlEqStrain( { tNumHalfCoeffCurrent, tNumCoeffCurrent - 1 }, { 0, 0 } ) =
                    2.0 * tLeaderDOFHatNlEqStrain( { 0, tNumHalfCoeffCurrent - 1 }, { 0, 0 } );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            // fill coefficients for master FI
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeader->mSet->mLeaderFIManager = &tFIManager;

            // set CM field interpolator manager
            tCMLeader->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // std::cout << "Integration point " << iGP << " / " << tNumGPs << std::endl;

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // reset CM evaluation flags
                tCMLeader->reset_eval_flags();

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeader->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // std::cout << "Requested dof derivative point " << iRequestedDof << " / " << tRequestedLeaderGlobalDofTypes.size() << std::endl;

                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // cast constitutive model base class pointer to elasticity damage constitutive model
                    CM_Struc_Linear_Isotropic_Damage* tCMLeaderPtr =
                            dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMLeader.get() );

                    // equivalent strain
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > teqstrain = tCMLeaderPtr->equivalent_strain();
                    // print( teqstrain, "teqstrain" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdeqstraindu = tCMLeaderPtr->dEqStraindu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdeqstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::EQSTRAIN,
                            tdeqstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump,
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdeqstraindu analytical and FD match " << std::endl;
                    // print( tdeqstraindu, "tdeqstraindu" );
                    // print( tdeqstrainduFD, "tdeqstrainduFD" );
                    bool tCheckdEqStraindu = fem::check( tdeqstraindu, tdeqstrainduFD, tEpsilon );
                    REQUIRE( tCheckdEqStraindu );

                    // nonlocal equivalent strain history
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > thistory = tCMLeaderPtr->history();
                    // print( thistory, "thistory" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdhistorydu = tCMLeaderPtr->dHistorydu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdhistoryduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::HISTORY,
                            tdhistoryduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdhistorydu analytical and FD match " << std::endl;
                    // print( tdhistorydu, "tdhistorydu" );
                    // print( tdhistoryduFD, "tdhistoryduFD" );
                    bool tCheckdHistorydu = fem::check( tdhistorydu, tdhistoryduFD, tEpsilon );
                    REQUIRE( tCheckdHistorydu );

                    // damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tdamage = tCMLeaderPtr->damage();
                    // print( tdamage, "tdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tddamagedu = tCMLeaderPtr->dDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tddamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::DAMAGE,
                            tddamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tddamagedu, "tddamagedu" );
                    // print( tddamageduFD, "tddamageduFD" );
                    bool tCheckdDamagedu = fem::check( tddamagedu, tddamageduFD, tEpsilon );
                    REQUIRE( tCheckdDamagedu );

                    // smooth damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tsmoothdamage = tCMLeaderPtr->smooth_damage();
                    // print( tsmoothdamage, "tsmoothdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tdsmoothdamagedu = tCMLeaderPtr->dSmoothDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tdsmoothdamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::SMOOTH_DAMAGE,
                            tdsmoothdamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tdsmoothdamagedu, "tdsmoothdamagedu" );
                    // print( tdsmoothdamageduFD, "tdsmoothdamageduFD" );
                    bool tCheckdSmoothDamagedu = fem::check( tdsmoothdamagedu, tdsmoothdamageduFD, tEpsilon );
                    REQUIRE( tCheckdSmoothDamagedu );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate flux
                    Matrix< DDRMat > tflux = tCMLeader->flux();
                    // print( tflux, "tflux" );

                    // evaluate derivative of flux
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
                    // std::cout << "check tdfluxdu analytical and FD match " << std::endl;
                    // print( tdfluxdu, "tdfluxdu" );
                    // print( tdfluxduFD, "tdfluxduFD" );
                    bool tCheckdFluxdu = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckdFluxdu );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate traction
                    Matrix< DDRMat > ttraction = tCMLeader->traction( tNormal );
                    // print( ttraction, "ttraction" );

                    // evaluate derivative of traction
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
                    // std::cout << "check tdtractiondu analytical and FD match " << std::endl;
                    // print( tdtractiondu, "tdtractiondu" );
                    // print( tdtractionduFD, "tdtractionduFD" );
                    bool tCheckdTractiondu = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckdTractiondu );

                    // test traction -- only displacement as test dof type !!!
                    //------------------------------------------------------------------------------
                    // get the test dof type
                    Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( 0 );

                    // evaluate test traction
                    Matrix< DDRMat > ttesttraction = tCMLeader->testTraction(
                            tNormal,
                            tDofTest );
                    // print( ttesttraction, "ttesttraction" );

                    // evaluate derivative test traction
                    Matrix< DDRMat > tdtesttractiondu = tCMLeader->dTestTractiondDOF(
                            tDofDerivative,
                            tNormal,
                            tJump,
                            tDispDofTypes( 0 ) );

                    // evaluate dtesttractiondu by FD
                    Matrix< DDRMat > tdtesttractionduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TEST_TRACTION,
                            tdtesttractionduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDispDofTypes( 0 ),
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdtesttractiondu analytical and FD match " << std::endl;
                    // print( tdtesttractiondu, "tdtesttractiondu" );
                    // print( tdtesttractionduFD, "tdtesttractionduFD" );
                    bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                    REQUIRE( tCheckTestTraction );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate strain
                    Matrix< DDRMat > tstrain = tCMLeader->strain();
                    // print( tstrain, "tstrain" );

                    // evaluate derivative of strain
                    Matrix< DDRMat > tdstraindu = tCMLeader->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdstraindu analytical and FD match " << std::endl;
                    // print( tdstraindu, "tdstraindu" );
                    // print( tdstrainduFD, "tdstrainduFD" );
                    bool tCheckdStraindu = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckdStraindu );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "CM_Struc_Linear_Isotropic_Damage_112_Threshold", "[CM_Struc_Lin_Iso_Dam_112_Threshold]" )
{
    // Damage setup - 1,0,1
    // Local equivalent strain - 1 - tensile/compressive strength
    // Damage law - 1 - exponential
    // Smoothing law - 2 - corrected ks smoothing
    moris::Cell< Matrix< DDRMat > > tDamageParameters = { //
        { { 1.0, 4.0 } },                                 //
        { { 1.0, 1.0, 0.95, 100.0 } },                    //
        { { 2.0, 7.0 } }
    };

    // define an epsilon environment
    real tEpsilon = 1.0E-4;

    // define a perturbation relative size
    real tPerturbation = 1.0E-5;

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
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs     = { { 8, 18, 32 }, { 16, 54, 128 } };
    Matrix< DDRMat > tNumHalfCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes       = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes           = { tDispDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropNu, "PoissonRatio" );
    tCMLeader->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // set size and populate the set master dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // std::cout << "Space dimension " << iSpaceDim << std::endl;

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

        // set space dimensions for property, CM, and SP
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );
        tCMLeader->set_parameters( tDamageParameters );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // std::cout << "Interpolation order " << iInterpOrder << std::endl;

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

            // create a cell of field interpolators current
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // get number of coefficients
            uint tNumHalfCoeffCurrent = tNumHalfCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
            uint tNumCoeffCurrent     = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // create the field interpolator for displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 1.0e-3 * tLeaderDOFHatNlEqStrain;
            tLeaderDOFHatNlEqStrain( { tNumHalfCoeffCurrent, tNumCoeffCurrent - 1 }, { 0, 0 } ) =
                    2.0 * tLeaderDOFHatNlEqStrain( { 0, tNumHalfCoeffCurrent - 1 }, { 0, 0 } );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            // fill coefficients for master FI
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum ge::PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMLeader->mSet->mLeaderFIManager = &tFIManager;

            // set CM field interpolator manager
            tCMLeader->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // std::cout << "Integration point " << iGP << " / " << tNumGPs << std::endl;

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // reset CM evaluation flags
                tCMLeader->reset_eval_flags();

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
                        tCMLeader->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // std::cout << "Requested dof derivative point " << iRequestedDof << " / " << tRequestedLeaderGlobalDofTypes.size() << std::endl;

                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // cast constitutive model base class pointer to elasticity damage constitutive model
                    CM_Struc_Linear_Isotropic_Damage* tCMLeaderPtr =
                            dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMLeader.get() );

                    // equivalent strain
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > teqstrain = tCMLeaderPtr->equivalent_strain();
                    // print( teqstrain, "teqstrain" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdeqstraindu = tCMLeaderPtr->dEqStraindu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdeqstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::EQSTRAIN,
                            tdeqstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump,
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdeqstraindu analytical and FD match " << std::endl;
                    // print( tdeqstraindu, "tdeqstraindu" );
                    // print( tdeqstrainduFD, "tdeqstrainduFD" );
                    bool tCheckdEqStraindu = fem::check( tdeqstraindu, tdeqstrainduFD, tEpsilon );
                    REQUIRE( tCheckdEqStraindu );

                    // nonlocal equivalent strain history
                    //------------------------------------------------------------------------------
                    // evaluate equivalent strain
                    Matrix< DDRMat > thistory = tCMLeaderPtr->history();
                    // print( thistory, "thistory" );

                    // evaluate derivative of equivalent strain
                    Matrix< DDRMat > tdhistorydu = tCMLeaderPtr->dHistorydu( tDofDerivative );

                    // evaluate derivative of equivalent strain by FD
                    Matrix< DDRMat > tdhistoryduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::HISTORY,
                            tdhistoryduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tdhistorydu analytical and FD match " << std::endl;
                    // print( tdhistorydu, "tdhistorydu" );
                    // print( tdhistoryduFD, "tdhistoryduFD" );
                    bool tCheckdHistorydu = fem::check( tdhistorydu, tdhistoryduFD, tEpsilon );
                    REQUIRE( tCheckdHistorydu );

                    // damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tdamage = tCMLeaderPtr->damage();
                    // print( tdamage, "tdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tddamagedu = tCMLeaderPtr->dDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tddamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::DAMAGE,
                            tddamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tddamagedu, "tddamagedu" );
                    // print( tddamageduFD, "tddamageduFD" );
                    bool tCheckdDamagedu = fem::check( tddamagedu, tddamageduFD, tEpsilon );
                    REQUIRE( tCheckdDamagedu );

                    // smooth damage
                    //------------------------------------------------------------------------------
                    // evaluate damage
                    Matrix< DDRMat > tsmoothdamage = tCMLeaderPtr->smooth_damage();
                    // print( tsmoothdamage, "tsmoothdamage" );

                    // evaluate derivative of damage wrt to dof
                    Matrix< DDRMat > tdsmoothdamagedu = tCMLeaderPtr->dSmoothDamagedu( tDofDerivative );

                    // evaluate derivative of damage wrt to dof by FD
                    Matrix< DDRMat > tdsmoothdamageduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::SMOOTH_DAMAGE,
                            tdsmoothdamageduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_3_CENTRAL );

                    // check that analytical and FD match
                    // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
                    // print( tdsmoothdamagedu, "tdsmoothdamagedu" );
                    // print( tdsmoothdamageduFD, "tdsmoothdamageduFD" );
                    bool tCheckdSmoothDamagedu = fem::check( tdsmoothdamagedu, tdsmoothdamageduFD, tEpsilon );
                    REQUIRE( tCheckdSmoothDamagedu );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate flux
                    Matrix< DDRMat > tflux = tCMLeader->flux();
                    // print( tflux, "tflux" );

                    // evaluate derivative of flux
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
                    // std::cout << "check tdfluxdu analytical and FD match " << std::endl;
                    // print( tdfluxdu, "tdfluxdu" );
                    // print( tdfluxduFD, "tdfluxduFD" );
                    bool tCheckdFluxdu = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckdFluxdu );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate traction
                    Matrix< DDRMat > ttraction = tCMLeader->traction( tNormal );
                    // print( ttraction, "ttraction" );

                    // evaluate derivative of traction
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
                    // std::cout << "check tdtractiondu analytical and FD match " << std::endl;
                    // print( tdtractiondu, "tdtractiondu" );
                    // print( tdtractionduFD, "tdtractionduFD" );
                    bool tCheckdTractiondu = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckdTractiondu );

                    // test traction -- only displacement as test dof type !!!
                    //------------------------------------------------------------------------------
                    // get the test dof type
                    Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( 0 );

                    // evaluate test traction
                    Matrix< DDRMat > ttesttraction = tCMLeader->testTraction(
                            tNormal,
                            tDofTest );
                    // print( ttesttraction, "ttesttraction" );

                    // evaluate derivative test traction
                    Matrix< DDRMat > tdtesttractiondu = tCMLeader->dTestTractiondDOF(
                            tDofDerivative,
                            tNormal,
                            tJump,
                            tDispDofTypes( 0 ) );

                    // evaluate dtesttractiondu by FD
                    Matrix< DDRMat > tdtesttractionduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::TEST_TRACTION,
                            tdtesttractionduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDispDofTypes( 0 ),
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdtesttractiondu analytical and FD match " << std::endl;
                    // print( tdtesttractiondu, "tdtesttractiondu" );
                    // print( tdtesttractionduFD, "tdtesttractionduFD" );
                    bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                    REQUIRE( tCheckTestTraction );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate strain
                    Matrix< DDRMat > tstrain = tCMLeader->strain();
                    // print( tstrain, "tstrain" );

                    // evaluate derivative of strain
                    Matrix< DDRMat > tdstraindu = tCMLeader->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,
                            tNormal,
                            tJump );

                    // check that analytical and FD match
                    // std::cout << "check tdstraindu analytical and FD match " << std::endl;
                    // print( tdstraindu, "tdstraindu" );
                    // print( tdstrainduFD, "tdstrainduFD" );
                    bool tCheckdStraindu = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckdStraindu );
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/

// FIXME zero test to be reintroduced
// TEST_CASE( "CM_Struc_Linear_Isotropic_Damage_112_Zero", "[CM_Struc_Lin_Iso_Dam_112_Zero]" )
//{
//     // Damage setup - 1,0,1
//     // Local equivalent strain - 1 - tensile/compressive strength
//     // Damage law - 0 - linear
//     // Smoothing law - 1 - ks smoothing
//     moris::Cell< Matrix< DDRMat > > tDamageParameters = { //
//         { { 1.0, 4.0 } },                                 //
//         { { 1.0, 1.0e-3, 0.95, 100.0 } },                 //
//         { { 2.0, 7.0 } }
//     };
//
//     // define an epsilon environment
//     real tEpsilon = 1.0E-4;
//
//     // define a perturbation relative size
//     real tPerturbation = 1.0E-5;
//
//     // init geometry inputs
//     //------------------------------------------------------------------------------
//     // create geometry type
//     mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;
//
//     // create space coeff xHat
//     Matrix< DDRMat > tXHat;
//
//     // create list of interpolation orders
//     moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
//         mtk::Interpolation_Order::LINEAR,
//         mtk::Interpolation_Order::QUADRATIC,
//         mtk::Interpolation_Order::CUBIC
//     };
//
//     // create list of integration orders
//     moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
//         mtk::Integration_Order::QUAD_2x2,
//         mtk::Integration_Order::HEX_2x2x2
//     };
//
//     // create list with number of coeffs
//     Matrix< DDRMat > tNumCoeffs     = { { 8, 18, 32 }, { 16, 54, 128 } };
//     Matrix< DDRMat > tNumHalfCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };
//
//     // dof type list
//     moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes       = { { MSI::Dof_Type::UX } };
//     moris::Cell< moris::Cell< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
//     moris::Cell< moris::Cell< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
//     moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes           = { tDispDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };
//
//     // create the properties
//     std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
//     tPropEMod->set_parameters( { { { 1.0 } } } );
//
//     std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
//     tPropNu->set_parameters( { { { 0.3 } } } );
//
//     // define constitutive models
//     fem::CM_Factory tCMFactory;
//
//     std::shared_ptr< fem::Constitutive_Model > tCMLeader =
//             tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
//     tCMLeader->set_dof_type_list( tDofTypes );
//     tCMLeader->set_property( tPropEMod, "YoungsModulus" );
//     tCMLeader->set_property( tPropNu, "PoissonRatio" );
//     tCMLeader->set_local_properties();
//
//     // set a fem set pointer
//     MSI::Equation_Set* tSet = new fem::Set();
//     tCMLeader->set_set_pointer( static_cast< fem::Set* >( tSet ) );
//
//     // set size for the set EqnObjDofTypeList
//     tCMLeader->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );
//
//     // set size and populate the set dof type map
//     tCMLeader->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
//     tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
//     tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
//     tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;
//
//     // set size and populate the set master dof type map
//     tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
//     tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
//     tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
//     tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;
//
//     // build global dof type list
//     tCMLeader->get_global_dof_type_list();
//
//     // loop on the space dimension
//     for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
//     {
//         // std::cout << "Space dimension " << iSpaceDim << std::endl;
//
//         // create and set normal
//         Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
//         tNormal = tNormal / norm( tNormal );
//
//         // create the jump
//         Matrix< DDRMat > tJump( iSpaceDim, 1, 10.0 );
//
//         // switch on space dimension
//         switch ( iSpaceDim )
//         {
//             case 2:
//             {
//                 // set geometry type
//                 tGeometryType = mtk::Geometry_Type::QUAD;
//
//                 // fill space coeff xHat
//                 tXHat = { { 0.0, 0.0 },
//                     { 1.0, 0.0 },
//                     { 1.0, 1.0 },
//                     { 0.0, 1.0 } };
//
//                 // set velocity dof types
//                 tDispDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };
//
//                 break;
//             }
//             case 3:
//             {
//                 // set geometry type
//                 tGeometryType = mtk::Geometry_Type::HEX;
//
//                 // fill space coeff xHat
//                 tXHat = { { 0.0, 0.0, 0.0 },
//                     { 1.0, 0.0, 0.0 },
//                     { 1.0, 1.0, 0.0 },
//                     { 0.0, 1.0, 0.0 },
//                     { 0.0, 0.0, 1.0 },
//                     { 1.0, 0.0, 1.0 },
//                     { 1.0, 1.0, 1.0 },
//                     { 0.0, 1.0, 1.0 } };
//
//                 // set velocity dof types
//                 tDispDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } };
//
//                 break;
//             }
//             default:
//             {
//                 MORIS_ERROR( false, " QUAD or HEX only." );
//                 break;
//             }
//         }
//
//         // space and time geometry interpolators
//         //------------------------------------------------------------------------------
//         // create a space geometry interpolation rule
//         mtk::Interpolation_Rule tGIRule(
//                 tGeometryType,
//                 mtk::Interpolation_Type::LAGRANGE,
//                 mtk::Interpolation_Order::LINEAR,
//                 mtk::Interpolation_Type::LAGRANGE,
//                 mtk::Interpolation_Order::LINEAR );
//
//         // create a space time geometry interpolator
//         Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );
//
//         // create time coeff tHat
//         Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };
//
//         // set the coefficients xHat, tHat
//         tGI.set_coeff( tXHat, tTHat );
//
//         // set space dimensions for property, CM, and SP
//         tCMLeader->set_space_dim( iSpaceDim );
//         tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );
//         tCMLeader->set_parameters( tDamageParameters );
//
//         // loop on the interpolation order
//         for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
//         {
//             // std::cout << "Interpolation order " << iInterpOrder << std::endl;
//
//             // integration points
//             //------------------------------------------------------------------------------
//             // get an integration order
//             mtk::Integration_Order tIntegrationOrder = tIntegrationOrders( iSpaceDim - 2 );
//
//             // create an integration rule
//             mtk::Integration_Rule tIntegrationRule(
//                     tGeometryType,
//                     mtk::Integration_Type::GAUSS,
//                     tIntegrationOrder,
//                     mtk::Geometry_Type::LINE,
//                     mtk::Integration_Type::GAUSS,
//                     mtk::Integration_Order::BAR_2 );
//
//             // create an integrator
//             mtk::Integrator tIntegrator( tIntegrationRule );
//
//             // get integration points
//             Matrix< DDRMat > tIntegPoints;
//             tIntegrator.get_points( tIntegPoints );
//
//             // field interpolators
//             //------------------------------------------------------------------------------
//             // create an interpolation order
//             mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );
//
//             // create a space time interpolation rule
//             mtk::Interpolation_Rule tFIRule(
//                     tGeometryType,
//                     mtk::Interpolation_Type::LAGRANGE,
//                     tInterpolationOrder,
//                     mtk::Interpolation_Type::LAGRANGE,
//                     mtk::Interpolation_Order::LINEAR );
//
//             // create a cell of field interpolators current
//             Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );
//
//             // get number of coefficients
//             uint tNumHalfCoeffCurrent = tNumHalfCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
//             uint tNumCoeffCurrent     = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
//
//             // create the field interpolator for displacement
//             tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
//             // fill coefficients for master FI
//             Matrix< DDRMat > tLeaderDOFHatDispl;
//             fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
//             tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );
//
//             // create the field interpolator for nonlocal equivalent strain
//             tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
//             // fill coefficients for master FI
//             Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
//             fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
//             tLeaderDOFHatNlEqStrain = 0.0 * tLeaderDOFHatNlEqStrain;
//             tLeaderDOFHatNlEqStrain( { tNumHalfCoeffCurrent, tNumCoeffCurrent - 1 }, { 0, 0 } ) =
//                     2.0 * tLeaderDOFHatNlEqStrain( { 0, tNumHalfCoeffCurrent - 1 }, { 0, 0 } );
//             tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );
//
//             // create the field interpolator for nonlocal equivalent strain
//             tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
//             // fill coefficients for master FI
//             tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatNlEqStrain );
//
//             // create a field interpolator manager
//             moris::Cell< moris::Cell< enum ge::PDV_Type > >        tDummyDv;
//             moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
//             Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
//
//             // populate the field interpolator manager
//             tFIManager.mFI                     = tLeaderFIs;
//             tFIManager.mIPGeometryInterpolator = &tGI;
//             tFIManager.mIGGeometryInterpolator = &tGI;
//
//             // set the interpolator manager to the set
//             tCMLeader->mSet->mLeaderFIManager = &tFIManager;
//
//             // set CM field interpolator manager
//             tCMLeader->set_field_interpolator_manager( &tFIManager );
//
//             uint tNumGPs = tIntegPoints.n_cols();
//             for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
//             {
//                 // std::cout << "Integration point " << iGP << " / " << tNumGPs << std::endl;
//
//                 // create evaluation point xi, tau
//                 Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );
//
//                 // set integration point
//                 tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );
//
//                 // reset CM evaluation flags
//                 tCMLeader->reset_eval_flags();
//
//                 // populate the requested master dof type for CM
//                 moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedLeaderGlobalDofTypes =
//                         tCMLeader->get_global_dof_type_list();
//
//                 // populate the test master dof type for CM
//                 moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
//                         tCMLeader->get_dof_type_list();
//
//                 // loop over requested dof type
//                 for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
//                 {
//                     // std::cout << "Requested dof derivative point " << iRequestedDof << " / " << tRequestedLeaderGlobalDofTypes.size() << std::endl;
//
//                     // derivative dof type
//                     Cell< MSI::Dof_Type > tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );
//
//                     // cast constitutive model base class pointer to elasticity damage constitutive model
//                     CM_Struc_Linear_Isotropic_Damage* tCMLeaderPtr =
//                             dynamic_cast< CM_Struc_Linear_Isotropic_Damage* >( tCMLeader.get() );
//
//                     // equivalent strain
//                     //------------------------------------------------------------------------------
//                     // evaluate equivalent strain
//                     Matrix< DDRMat > teqstrain = tCMLeaderPtr->equivalent_strain();
//                     print( teqstrain, "teqstrain" );
//
//                     // evaluate derivative of equivalent strain
//                     Matrix< DDRMat > tdeqstraindu = tCMLeaderPtr->dEqStraindu( tDofDerivative );
//
//                     // evaluate derivative of equivalent strain by FD
//                     Matrix< DDRMat > tdeqstrainduFD;
//                     tCMLeader->eval_derivative_FD(
//                             CM_Request_Type::EQSTRAIN,
//                             tdeqstrainduFD,
//                             tDofDerivative,
//                             tPerturbation,
//                             tDofDerivative,
//                             tNormal,
//                             tJump,
//                             FDScheme_Type::POINT_3_CENTRAL );
//
//                     // check that analytical and FD match
//                     // std::cout << "check tdeqstraindu analytical and FD match " << std::endl;
//                     // print( tdeqstraindu, "tdeqstraindu" );
//                     // print( tdeqstrainduFD, "tdeqstrainduFD" );
//                     bool tCheckdEqStraindu = fem::check( tdeqstraindu, tdeqstrainduFD, tEpsilon );
//                     REQUIRE( tCheckdEqStraindu );
//
//                     // nonlocal equivalent strain history
//                     //------------------------------------------------------------------------------
//                     // evaluate equivalent strain
//                     Matrix< DDRMat > thistory = tCMLeaderPtr->history();
//                     print( thistory, "thistory" );
//
//                     // evaluate derivative of equivalent strain
//                     Matrix< DDRMat > tdhistorydu = tCMLeaderPtr->dHistorydu( tDofDerivative );
//
//                     // evaluate derivative of equivalent strain by FD
//                     Matrix< DDRMat > tdhistoryduFD;
//                     tCMLeader->eval_derivative_FD(
//                             CM_Request_Type::HISTORY,
//                             tdhistoryduFD,
//                             tDofDerivative,
//                             tPerturbation,
//                             tDofDerivative,    // dummy
//                             tNormal,           // dummy
//                             tJump,             // dummy
//                             FDScheme_Type::POINT_3_CENTRAL );
//
//                     // check that analytical and FD match
//                     // std::cout << "check tdhistorydu analytical and FD match " << std::endl;
//                     // print( tdhistorydu, "tdhistorydu" );
//                     // print( tdhistoryduFD, "tdhistoryduFD" );
//                     bool tCheckdHistorydu = fem::check( tdhistorydu, tdhistoryduFD, tEpsilon );
//                     REQUIRE( tCheckdHistorydu );
//
//                     // damage
//                     //------------------------------------------------------------------------------
//                     // evaluate damage
//                     Matrix< DDRMat > tdamage = tCMLeaderPtr->damage();
//                     // print( tdamage, "tdamage" );
//
//                     // evaluate derivative of damage wrt to dof
//                     Matrix< DDRMat > tddamagedu = tCMLeaderPtr->dDamagedu( tDofDerivative );
//
//                     // evaluate derivative of damage wrt to dof by FD
//                     Matrix< DDRMat > tddamageduFD;
//                     tCMLeader->eval_derivative_FD(
//                             CM_Request_Type::DAMAGE,
//                             tddamageduFD,
//                             tDofDerivative,
//                             tPerturbation,
//                             tDofDerivative,    // dummy
//                             tNormal,           // dummy
//                             tJump,             // dummy
//                             FDScheme_Type::POINT_3_CENTRAL );
//
//                     // check that analytical and FD match
//                     // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
//                     // print( tddamagedu, "tddamagedu" );
//                     // print( tddamageduFD, "tddamageduFD" );
//                     bool tCheckdDamagedu = fem::check( tddamagedu, tddamageduFD, tEpsilon );
//                     REQUIRE( tCheckdDamagedu );
//
//                     // smooth damage
//                     //------------------------------------------------------------------------------
//                     // evaluate damage
//                     Matrix< DDRMat > tsmoothdamage = tCMLeaderPtr->smooth_damage();
//                     // print( tsmoothdamage, "tsmoothdamage" );
//
//                     // evaluate derivative of damage wrt to dof
//                     Matrix< DDRMat > tdsmoothdamagedu = tCMLeaderPtr->dSmoothDamagedu( tDofDerivative );
//
//                     // evaluate derivative of damage wrt to dof by FD
//                     Matrix< DDRMat > tdsmoothdamageduFD;
//                     tCMLeader->eval_derivative_FD(
//                             CM_Request_Type::SMOOTH_DAMAGE,
//                             tdsmoothdamageduFD,
//                             tDofDerivative,
//                             tPerturbation,
//                             tDofDerivative,    // dummy
//                             tNormal,           // dummy
//                             tJump,             // dummy
//                             FDScheme_Type::POINT_3_CENTRAL );
//
//                     // check that analytical and FD match
//                     // std::cout << "check tCheckdDamagedu analytical and FD match " << std::endl;
//                     // print( tdsmoothdamagedu, "tdsmoothdamagedu" );
//                     // print( tdsmoothdamageduFD, "tdsmoothdamageduFD" );
//                     bool tCheckdSmoothDamagedu = fem::check( tdsmoothdamagedu, tdsmoothdamageduFD, tEpsilon );
//                     REQUIRE( tCheckdSmoothDamagedu );
//
//                     // flux
//                     //------------------------------------------------------------------------------
//                     // evaluate flux
//                     Matrix< DDRMat > tflux = tCMLeader->flux();
//                     // print( tflux, "tflux" );
//
//                     // evaluate derivative of flux
//                     Matrix< DDRMat > tdfluxdu = tCMLeader->dFluxdDOF( tDofDerivative );
//
//                     // evaluate dfluxdu by FD
//                     Matrix< DDRMat > tdfluxduFD;
//                     tCMLeader->eval_derivative_FD(
//                             CM_Request_Type::FLUX,
//                             tdfluxduFD,
//                             tDofDerivative,
//                             tPerturbation,
//                             tDofDerivative,    // dummy
//                             tNormal,           // dummy
//                             tJump );           // dummy
//
//                     // check that analytical and FD match
//                     // std::cout << "check tdfluxdu analytical and FD match " << std::endl;
//                     // print( tdfluxdu, "tdfluxdu" );
//                     // print( tdfluxduFD, "tdfluxduFD" );
//                     bool tCheckdFluxdu = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
//                     REQUIRE( tCheckdFluxdu );
//
//                     // traction
//                     //------------------------------------------------------------------------------
//                     // evaluate traction
//                     Matrix< DDRMat > ttraction = tCMLeader->traction( tNormal );
//                     // print( ttraction, "ttraction" );
//
//                     // evaluate derivative of traction
//                     Matrix< DDRMat > tdtractiondu = tCMLeader->dTractiondDOF( tDofDerivative, tNormal );
//
//                     // evaluate dtractiondu by FD
//                     Matrix< DDRMat > tdtractionduFD;
//                     tCMLeader->eval_derivative_FD(
//                             CM_Request_Type::TRACTION,
//                             tdtractionduFD,
//                             tDofDerivative,
//                             tPerturbation,
//                             tDofDerivative,    // dummy
//                             tNormal,
//                             tJump );
//
//                     // check that analytical and FD match
//                     // std::cout << "check tdtractiondu analytical and FD match " << std::endl;
//                     // print( tdtractiondu, "tdtractiondu" );
//                     // print( tdtractionduFD, "tdtractionduFD" );
//                     bool tCheckdTractiondu = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
//                     REQUIRE( tCheckdTractiondu );
//
//                     // test traction -- only displacement as test dof type !!!
//                     //------------------------------------------------------------------------------
//                     // get the test dof type
//                     Cell< MSI::Dof_Type > tDofTest = tLeaderDofTypes( 0 );
//
//                     // evaluate test traction
//                     Matrix< DDRMat > ttesttraction = tCMLeader->testTraction(
//                             tNormal,
//                             tDofTest );
//                     // print( ttesttraction, "ttesttraction" );
//
//                     // evaluate derivative test traction
//                     Matrix< DDRMat > tdtesttractiondu = tCMLeader->dTestTractiondDOF(
//                             tDofDerivative,
//                             tNormal,
//                             tJump,
//                             tDispDofTypes( 0 ) );
//
//                     // evaluate dtesttractiondu by FD
//                     Matrix< DDRMat > tdtesttractionduFD;
//                     tCMLeader->eval_derivative_FD(
//                             CM_Request_Type::TEST_TRACTION,
//                             tdtesttractionduFD,
//                             tDofDerivative,
//                             tPerturbation,
//                             tDispDofTypes( 0 ),
//                             tNormal,
//                             tJump );
//
//                     // check that analytical and FD match
//                     // std::cout << "check tdtesttractiondu analytical and FD match " << std::endl;
//                     // print( tdtesttractiondu, "tdtesttractiondu" );
//                     // print( tdtesttractionduFD, "tdtesttractionduFD" );
//                     bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
//                     REQUIRE( tCheckTestTraction );
//
//                     // strain
//                     //------------------------------------------------------------------------------
//                     // evaluate strain
//                     Matrix< DDRMat > tstrain = tCMLeader->strain();
//                     // print( tstrain, "tstrain" );
//
//                     // evaluate derivative of strain
//                     Matrix< DDRMat > tdstraindu = tCMLeader->dStraindDOF( tDofDerivative );
//
//                     // evaluate dstraindu by FD
//                     Matrix< DDRMat > tdstrainduFD;
//                     tCMLeader->eval_derivative_FD(
//                             CM_Request_Type::STRAIN,
//                             tdstrainduFD,
//                             tDofDerivative,
//                             tPerturbation,
//                             tDofDerivative,
//                             tNormal,
//                             tJump );
//
//                     // check that analytical and FD match
//                     // std::cout << "check tdstraindu analytical and FD match " << std::endl;
//                     // print( tdstraindu, "tdstraindu" );
//                     // print( tdstrainduFD, "tdstrainduFD" );
//                     bool tCheckdStraindu = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
//                     REQUIRE( tCheckdStraindu );
//                 }
//             }
//             // clean up
//             tLeaderFIs.clear();
//         }
//     }
// } /*END_TEST_CASE*/
