/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM_Struct_Linear_MT.cpp
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
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Diffusion_UT.cpp"

using namespace moris;
using namespace fem;

void
tDOFFunction( moris::Matrix< moris::DDRMat >&          aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    aPropMatrix = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->val()( 0 );
}

void
tDOFFunctionDer( moris::Matrix< moris::DDRMat >&       aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{

    aPropMatrix = aFIManager->get_field_interpolators_for_type( MSI::Dof_Type::PHID )->N();
}

TEST_CASE( "CM_Struc_Linear_MT", "[CM_Struc_Lin_MT]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-4;

    // define a perturbation relative size
    real tPerturbation = 2.0E-7;

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
    Vector< Vector< MSI::Dof_Type > > tDispDofTypes      = { { MSI::Dof_Type::UX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes          = { { MSI::Dof_Type::UX }, { MSI::Dof_Type::PHID } };
    Vector< Vector< MSI::Dof_Type > > tSecondaryDofTypes = { { MSI::Dof_Type::PHID } };

    // create a list of property names used in the MT model
    Vector< std::string > tPropertyNames = { "YoungsModulusMatrix", "PoissonRatioMatrix", "YoungsModulusFiber", "PoissonRatioFiber", "VolumeFraction", "AspectRatio", "OrientationInPlane", "OrientationOutPlane" };

    // initialize a cell of values for the defined propoetis
    Vector< real > tValues = { 1.0, 0.2, 2.0, 0.2, 0.1, 10, 1.0, 1.0 };

    // initialize a propet cell that will be used in the CM later
    Vector< std::shared_ptr< fem::Property > > tPropertyCell;

    // loop over the property names to create shared pointer properties
    for ( uint iCounter = 0; iCounter < tPropertyNames.size(); iCounter++ )
    {
        // create a shared pointer prop and set the evolution function
        std::shared_ptr< fem::Property > iProp = std::make_shared< fem::Property >();
        iProp->set_parameters( { { { tValues( iCounter ) } } } );
        iProp->set_val_function( tConstValFunc_Elast );

        if ( tPropertyNames( iCounter ) == "OrientationOutPlane" )
        {
            iProp->set_parameters( { { { tValues( iCounter ) } } } );
            iProp->set_val_function( tDOFFunction );
            iProp->set_dof_type_list( { { MSI::Dof_Type::PHID } } );
            iProp->set_dof_derivative_functions( { tDOFFunctionDer } );
        }
        // store the the created propety
        tPropertyCell.push_back( iProp );
    }

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_MT );
    tCMLeaderStrucLinIso->set_dof_type_list( { tDispDofTypes } );

    // set the properties of the CM using names and shared pointers created earlier
    for ( uint iCounter = 0; iCounter < tPropertyNames.size(); iCounter++ )
    {
        tCMLeaderStrucLinIso->set_property( tPropertyCell( iCounter ), tPropertyNames( iCounter ) );
    }

    tCMLeaderStrucLinIso->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tCMLeaderStrucLinIso->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMLeaderStrucLinIso->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeaderStrucLinIso->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderStrucLinIso->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 0;
    tCMLeaderStrucLinIso->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::PHID ) ) = 1;

    // set size and populate the set leader dof type map
    tCMLeaderStrucLinIso->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderStrucLinIso->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 0;
    tCMLeaderStrucLinIso->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::PHID ) ) = 1;

    map< enum MSI::Dof_Type, std::string > tDofTypeToNameMap = MSI::get_dof_type_name_map();

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
                tDispDofTypes      = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };
                tSecondaryDofTypes = { { MSI::Dof_Type::PHID } };
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
                tDispDofTypes      = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } };
                tSecondaryDofTypes = { { MSI::Dof_Type::PHID } };

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
        tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );

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

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator phi_d field
            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0;
            fill_that( tUHat0, iSpaceDim, iInterpOrder );
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tSecondaryDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tUHat0 );

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

            // uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < 1; iGP++ )
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
                    // evaluate dstraindu
                    Matrix< DDRMat > tdstraindu = tCMLeaderStrucLinIso->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMLeaderStrucLinIso->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdstrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump );           // dummy

                    // check that analytical and FD match
                    // std::cout << "check tdstraindu analytical and FD match " << std::endl;
                    // print( tdstraindu, "tdstraindu" );
                    // print( tdstrainduFD, "tdstrainduFD" );
                    bool tCheckdStraindu = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckdStraindu );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMLeaderStrucLinIso->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMLeaderStrucLinIso->eval_derivative_FD(
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
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMLeaderStrucLinIso->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMLeaderStrucLinIso->eval_derivative_FD(
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

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for ( uint iTestDof = 0; iTestDof < tLeaderDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Vector< MSI::Dof_Type > tDofTest = tLeaderDofTypes( iTestDof );

                        // evaluate derivative test traction
                        Matrix< DDRMat > tdtesttractiondu = tCMLeaderStrucLinIso->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tJump,
                                tDispDofTypes( 0 ) );

                        // evaluate dtesttractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMLeaderStrucLinIso->eval_derivative_FD(
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
                    }
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/

