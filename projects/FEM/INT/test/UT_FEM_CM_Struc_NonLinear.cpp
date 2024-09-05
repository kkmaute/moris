/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM_Struc_NonLinear.cpp
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

TEST_CASE( "CM_Struc_NonLinear",
        "[CM_Struc_NonLinear]" )
{
    // define an epsilon environment
    real tEpsilon = 5.0e-5;

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

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
    tCMLeader->set_dof_type_list( { tDispDofTypes } );
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
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

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
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRAIN );

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

            // create a Vector of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( 0.1 * tLeaderDOFHatVel );

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >             tDummyDv;
            Vector< Vector< mtk::Field_Type > >                tDummyField;
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
                // reset CM evaluation flags
                tCMLeader->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                //                // change in volume
                //                //------------------------------------------------------------------------------
                //                real tJ = tCMLeader->volume_change_jacobian();
                //                real tDG = tCMLeader->strain( DEFORMATION_GRADIENT );
                //

                // populate the requested leader dof type for CM
                Vector< Vector< MSI::Dof_Type > >
                        tRequestedLeaderGlobalDofTypes =
                                tCMLeader->get_global_dof_type_list();

                // populate the test leader dof type for CM
                Vector< Vector< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    const Vector< MSI::Dof_Type >& tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // strain
                    //------------------------------------------------------------------------------
                    //                    // evaluate the deformation gradient
                    //                    const Matrix< DDRMat >& tDGStrain = tCMLeader->strain( CM_Function_Type::DEFORMATION_GRADIENT );
                    //                    print( tDGStrain, "tDGStrain" );
                    //
                    //                    // evaluate the Green-Lagrange strain
                    //                    const Matrix< DDRMat >& tLGStrain = tCMLeader->strain( CM_Function_Type::LAGRANGIAN );
                    //                    print( tLGStrain, "tLGStrain" );
                    //
                    //                    //                    // evaluate the Euler-Almansi strain
                    //                    //                    const Matrix< DDRMat >& tEAStrain = tCMLeader->strain( CM_Function_Type::EULERIAN );
                    //                    //                    print( tEAStrain, "tEAStrain" );

                    // evaluate the derivative of the deformation gradient analytical
                    Matrix< DDRMat > tdDGStraindu =
                            tCMLeader->dStraindDOF(
                                    tDofDerivative,
                                    CM_Function_Type::DEFORMATION_GRADIENT );
                    // print( tdDGStraindu, "tdDGStraindu" );

                    // evaluate the derivative of the deformation gradient by FD
                    Matrix< DDRMat > tdDGStrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdDGStrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_5,
                            CM_Function_Type::DEFORMATION_GRADIENT );
                    // print( tdDGStrainduFD, "tdDGStrainduFD" );

                    // check that analytical and FD match
                    bool tCheckDF = fem::check( tdDGStraindu, tdDGStrainduFD, tEpsilon );
                    REQUIRE( tCheckDF );

                    // evaluate the derivative of the Green-Lagrange strain analytical
                    Matrix< DDRMat > tdLGStraindu =
                            tCMLeader->dStraindDOF(
                                    tDofDerivative,
                                    CM_Function_Type::LAGRANGIAN );
                    // print( tdLGStraindu, "tdLGStraindu" );

                    // evaluate the derivative of the Green-Lagrange strain by FD
                    Matrix< DDRMat > tdLGStrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdLGStrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_5,
                            CM_Function_Type::LAGRANGIAN );
                    // print( tdLGStrainduFD, "tdLGStrainduFD" );

                    // check that analytical and FD match
                    bool tCheckLG = fem::check( tdLGStraindu, tdLGStrainduFD, tEpsilon );
                    REQUIRE( tCheckLG );

                    //                    // evaluate the derivative of the Euler-Almansi strain analytical
                    //                    Matrix< DDRMat > tdEAStraindu =
                    //                            tCMLeader->dStraindDOF(
                    //                                    tDofDerivative,
                    //                                    CM_Function_Type::EULERIAN );
                    //                    print( tdEAStraindu, "tdEAStraindu" );
                    //
                    //                    // evaluate the derivative of the Euler-Almansi  strain by FD
                    //                    Matrix< DDRMat > tdEAStrainduFD;
                    //                    tCMLeader->eval_derivative_FD(
                    //                            CM_Request_Type::STRAIN,
                    //                            tdEAStrainduFD,
                    //                            tDofDerivative,
                    //                            tPerturbation,
                    //                            tDofDerivative,    // dummy
                    //                            tNormal,           // dummy
                    //                            tJump,             // dummy
                    //                            FDScheme_Type::POINT_5,
                    //                            CM_Function_Type::EULERIAN );
                    //                    print( tdEAStrainduFD, "tdEAStrainduFD" );
                    //
                    //                    // check that analytical and FD match
                    //                    bool tCheckDF = fem::check( tdEAStraindu, tdEAStrainduFD, tEpsilon );
                    //                    REQUIRE( tCheckDF );

                    // test strain
                    //------------------------------------------------------------------------------

                    // loop over test dof type
                    for ( uint iTestDof = 0; iTestDof < tLeaderDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type

                        //                        // evaluate test strain for deformation gradient
                        //                        const Matrix< DDRMat >& tDGTestStrain =
                        //                                tCMLeader->testStrain( CM_Function_Type::DEFORMATION_GRADIENT );
                        //                        //print( tDGTestStrain, "tDGTestStrain" );
                        //
                        //                        // evaluate test strain for Green-Lagrange strain
                        //                        const Matrix< DDRMat >& tLGTestStrain =
                        //                                tCMLeader->testStrain( CM_Function_Type::LAGRANGIAN );
                        //                        //print( tLGTestStrain, "tLGTestStrain" );

                        //                        // evaluate test strain for Euler-Almansi strain
                        //                        const Matrix< DDRMat >& tEATestStrain =
                        //                                tCMLeader->testStrain( CM_Function_Type::EULERIAN );
                        //                        //print( tEATestStrain, "tEATestStrain" );

                        //                        // evaluate derivative of the test strain for deformation gradient by analytical
                        //                        Matrix< DDRMat > tdDGTestStraindu =
                        //                                tCMLeader->dTestStraindDOF(
                        //                                        tDofDerivative,
                        //                                        CM_Function_Type::DEFORMATION_GRADIENT );
                        //                        print( tdDGTestStraindu, "tdDGTestStraindu" );
                        //
                        //                        // evaluate derivative of the test strain for deformation gradient by FD
                        //                        Matrix< DDRMat > tdDGTestStrainduFD;
                        //                        tCMLeader->eval_derivative_FD(
                        //                                CM_Request_Type::TEST_STRAIN,
                        //                                tdDGTestStrainduFD,
                        //                                tDofDerivative,
                        //                                tPerturbation,
                        //                                tDofDerivative,    // dummy
                        //                                tNormal,           // dummy
                        //                                tJump,             // dummy
                        //                                FDScheme_Type::POINT_5,
                        //                                CM_Function_Type::DEFORMATION_GRADIENT );
                        //                        print( tdDGTestStrainduFD, "tdDGTestStrainduFD" );
                        //
                        //                        // check that analytical and FD match
                        //                        bool tCheckDGTestStrain = fem::check( tdDGTestStraindu, tdDGTestStrainduFD, tEpsilon );
                        //                        REQUIRE( tCheckDGTestStrain );

                        //                        // evaluate derivative of the test strain for Green-Lagrange strain by analytical
                        //                        Matrix< DDRMat > tdLGTestStraindu =
                        //                                tCMLeader->dTestStraindDOF(
                        //                                        tDofDerivative,
                        //                                        CM_Function_Type::LAGRANGIAN );
                        //                        print( tdLGTestStraindu, "tdLGTestStraindu" );
                        //
                        //                        // evaluate derivative of the test strain for Green-Lagrange strain by FD
                        //                        Matrix< DDRMat > tdLGTestStrainduFD;
                        //                        tCMLeader->eval_derivative_FD(
                        //                                CM_Request_Type::TEST_STRAIN,
                        //                                tdLGTestStrainduFD,
                        //                                tDofDerivative,
                        //                                tPerturbation,
                        //                                tDofDerivative,    // dummy
                        //                                tNormal,           // dummy
                        //                                tJump,             // dummy
                        //                                FDScheme_Type::POINT_5,
                        //                                CM_Function_Type::LAGRANGIAN );
                        //                        print( tdLGTestStrainduFD, "tdLGTestStrainduFD" );
                        //
                        //                        // check that analytical and FD match
                        //                        bool tCheckLGTestStrain = fem::check( tdLGTestStraindu, tdLGTestStrainduFD, tEpsilon );
                        //                        REQUIRE( tCheckLGTestStrain );

                        //                        // evaluate derivative of the test strain for Euler-Almansi strain by analytical
                        //                        Matrix< DDRMat > tdEATestStraindu =
                        //                                tCMLeader->dTestStraindDOF(
                        //                                        tDofDerivative,
                        //                                        CM_Function_Type::EULERIAN );
                        //                        print( tdEATestStraindu, "tdEATestStraindu" );
                        //
                        //                        // evaluate derivative of the test strain for Euler-Almansi strain by FD
                        //                        Matrix< DDRMat > tdEATestStrainduFD;
                        //                        tCMLeader->eval_derivative_FD(
                        //                                CM_Request_Type::TEST_STRAIN,
                        //                                tdEATestStrainduFD,
                        //                                tDofDerivative,
                        //                                tPerturbation,
                        //                                tDofDerivative,    // dummy
                        //                                tNormal,           // dummy
                        //                                tJump,             // dummy
                        //                                FDScheme_Type::POINT_5,
                        //                                CM_Function_Type::EULERIAN );
                        //                        print( tdEATestStrainduFD, "tdEATestStrainduFD" );
                        //
                        //                        // check that analytical and FD match
                        //                        bool tCheckEATestStrain = fem::check( tdEATestStraindu, tdEATestStrainduFD, tEpsilon );
                        //                        REQUIRE( tCheckEATestStrain );
                    }
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "CM_Struc_NonLinear_EigenStrain",
        "[CM_Struc_NonLinear_EigenStrain]" )
{
    // define an epsilon environment
    real tEpsilon = 5.0e-5;

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
    std::shared_ptr< fem::Property > tPropEigenstrain = std::make_shared< fem::Property >();

    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
    tCMLeader->set_dof_type_list( { tDispDofTypes } );
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
    tCMLeader->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tCMLeader->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeader->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // build global dof type list
    tCMLeader->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( iSpaceDim, 1, 10.0 );

        // create the eigenstrain
        Matrix< DDRMat > tEigenstrain;

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

                // fill vector for eigenstrain
                tEigenstrain.set_size( 4, 1, 0.5 );

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

                // fill vector for eigenstrain
                tEigenstrain.set_size( 9, 1, 0.5 );

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
        // set eigenstrain values
        tPropEigenstrain->set_parameters( { tEigenstrain } );

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
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_model_type( fem::Model_Type::PLANE_STRAIN );

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

            // create a Vector of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( 0.1 * tLeaderDOFHatVel );

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > > tDummyDv;
            Vector< Vector< mtk::Field_Type > >    tDummyField;
            Field_Interpolator_Manager             tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

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
                // reset CM evaluation flags
                tCMLeader->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMLeader->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // populate the requested leader dof type for CM
                Vector< Vector< MSI::Dof_Type > >
                        tRequestedLeaderGlobalDofTypes =
                                tCMLeader->get_global_dof_type_list();

                // populate the test leader dof type for CM
                Vector< Vector< MSI::Dof_Type > > tLeaderDofTypes =
                        tCMLeader->get_dof_type_list();

                // loop over requested dof type
                for ( uint iRequestedDof = 0; iRequestedDof < tRequestedLeaderGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    const Vector< MSI::Dof_Type >& tDofDerivative = tRequestedLeaderGlobalDofTypes( iRequestedDof );

                    // strain
                    //------------------------------------------------------------------------------
                    //                    // evaluate the deformation gradient
                    //                    const Matrix< DDRMat >& tDGStrain = tCMLeader->strain( CM_Function_Type::DEFORMATION_GRADIENT );
                    //                    print( tDGStrain, "tDGStrain" );
                    //
                    //                    // evaluate the Green-Lagrange strain
                    //                    const Matrix< DDRMat >& tLGStrain = tCMLeader->strain( CM_Function_Type::LAGRANGIAN );
                    //                    print( tLGStrain, "tLGStrain" );
                    //
                    //                    //                    // evaluate the Euler-Almansi strain
                    //                    //                    const Matrix< DDRMat >& tEAStrain = tCMLeader->strain( CM_Function_Type::EULERIAN );
                    //                    //                    print( tEAStrain, "tEAStrain" );

                    // evaluate the derivative of the deformation gradient analytical
                    Matrix< DDRMat > tdDGStraindu =
                            tCMLeader->dStraindDOF(
                                    tDofDerivative,
                                    CM_Function_Type::DEFORMATION_GRADIENT );
                    // print( tdDGStraindu, "tdDGStraindu" );

                    // evaluate the derivative of the deformation gradient by FD
                    Matrix< DDRMat > tdDGStrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdDGStrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_5,
                            CM_Function_Type::DEFORMATION_GRADIENT );
                    // print( tdDGStrainduFD, "tdDGStrainduFD" );

                    // check that analytical and FD match
                    bool tCheckDF = fem::check( tdDGStraindu, tdDGStrainduFD, tEpsilon );
                    REQUIRE( tCheckDF );

                    // evaluate the derivative of the Green-Lagrange strain analytical
                    Matrix< DDRMat > tdLGStraindu =
                            tCMLeader->dStraindDOF(
                                    tDofDerivative,
                                    CM_Function_Type::LAGRANGIAN );
                    // print( tdLGStraindu, "tdLGStraindu" );

                    // evaluate the derivative of the Green-Lagrange strain by FD
                    Matrix< DDRMat > tdLGStrainduFD;
                    tCMLeader->eval_derivative_FD(
                            CM_Request_Type::STRAIN,
                            tdLGStrainduFD,
                            tDofDerivative,
                            tPerturbation,
                            tDofDerivative,    // dummy
                            tNormal,           // dummy
                            tJump,             // dummy
                            FDScheme_Type::POINT_5,
                            CM_Function_Type::LAGRANGIAN );
                    // print( tdLGStrainduFD, "tdLGStrainduFD" );

                    // check that analytical and FD match
                    bool tCheckLG = fem::check( tdLGStraindu, tdLGStrainduFD, tEpsilon );
                    REQUIRE( tCheckLG );

                    //                    // evaluate the derivative of the Euler-Almansi strain analytical
                    //                    Matrix< DDRMat > tdEAStraindu =
                    //                            tCMLeader->dStraindDOF(
                    //                                    tDofDerivative,
                    //                                    CM_Function_Type::EULERIAN );
                    //                    print( tdEAStraindu, "tdEAStraindu" );
                    //
                    //                    // evaluate the derivative of the Euler-Almansi  strain by FD
                    //                    Matrix< DDRMat > tdEAStrainduFD;
                    //                    tCMLeader->eval_derivative_FD(
                    //                            CM_Request_Type::STRAIN,
                    //                            tdEAStrainduFD,
                    //                            tDofDerivative,
                    //                            tPerturbation,
                    //                            tDofDerivative,    // dummy
                    //                            tNormal,           // dummy
                    //                            tJump,             // dummy
                    //                            FDScheme_Type::POINT_5,
                    //                            CM_Function_Type::EULERIAN );
                    //                    print( tdEAStrainduFD, "tdEAStrainduFD" );
                    //
                    //                    // check that analytical and FD match
                    //                    bool tCheckDF = fem::check( tdEAStraindu, tdEAStrainduFD, tEpsilon );
                    //                    REQUIRE( tCheckDF );

                    // test strain
                    //------------------------------------------------------------------------------

                    // loop over test dof type
                    for ( uint iTestDof = 0; iTestDof < tLeaderDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type

                        //                        // evaluate test strain for deformation gradient
                        //                        const Matrix< DDRMat >& tDGTestStrain =
                        //                                tCMLeader->testStrain( CM_Function_Type::DEFORMATION_GRADIENT );
                        //                        //print( tDGTestStrain, "tDGTestStrain" );
                        //
                        //                        // evaluate test strain for Green-Lagrange strain
                        //                        const Matrix< DDRMat >& tLGTestStrain =
                        //                                tCMLeader->testStrain( CM_Function_Type::LAGRANGIAN );
                        //                        //print( tLGTestStrain, "tLGTestStrain" );

                        //                        // evaluate test strain for Euler-Almansi strain
                        //                        const Matrix< DDRMat >& tEATestStrain =
                        //                                tCMLeader->testStrain( CM_Function_Type::EULERIAN );
                        //                        //print( tEATestStrain, "tEATestStrain" );

                        //                        // evaluate derivative of the test strain for deformation gradient by analytical
                        //                        Matrix< DDRMat > tdDGTestStraindu =
                        //                                tCMLeader->dTestStraindDOF(
                        //                                        tDofDerivative,
                        //                                        CM_Function_Type::DEFORMATION_GRADIENT );
                        //                        print( tdDGTestStraindu, "tdDGTestStraindu" );
                        //
                        //                        // evaluate derivative of the test strain for deformation gradient by FD
                        //                        Matrix< DDRMat > tdDGTestStrainduFD;
                        //                        tCMLeader->eval_derivative_FD(
                        //                                CM_Request_Type::TEST_STRAIN,
                        //                                tdDGTestStrainduFD,
                        //                                tDofDerivative,
                        //                                tPerturbation,
                        //                                tDofDerivative,    // dummy
                        //                                tNormal,           // dummy
                        //                                tJump,             // dummy
                        //                                FDScheme_Type::POINT_5,
                        //                                CM_Function_Type::DEFORMATION_GRADIENT );
                        //                        print( tdDGTestStrainduFD, "tdDGTestStrainduFD" );
                        //
                        //                        // check that analytical and FD match
                        //                        bool tCheckDGTestStrain = fem::check( tdDGTestStraindu, tdDGTestStrainduFD, tEpsilon );
                        //                        REQUIRE( tCheckDGTestStrain );

                        //                        // evaluate derivative of the test strain for Green-Lagrange strain by analytical
                        //                        Matrix< DDRMat > tdLGTestStraindu =
                        //                                tCMLeader->dTestStraindDOF(
                        //                                        tDofDerivative,
                        //                                        CM_Function_Type::LAGRANGIAN );
                        //                        print( tdLGTestStraindu, "tdLGTestStraindu" );
                        //
                        //                        // evaluate derivative of the test strain for Green-Lagrange strain by FD
                        //                        Matrix< DDRMat > tdLGTestStrainduFD;
                        //                        tCMLeader->eval_derivative_FD(
                        //                                CM_Request_Type::TEST_STRAIN,
                        //                                tdLGTestStrainduFD,
                        //                                tDofDerivative,
                        //                                tPerturbation,
                        //                                tDofDerivative,    // dummy
                        //                                tNormal,           // dummy
                        //                                tJump,             // dummy
                        //                                FDScheme_Type::POINT_5,
                        //                                CM_Function_Type::LAGRANGIAN );
                        //                        print( tdLGTestStrainduFD, "tdLGTestStrainduFD" );
                        //
                        //                        // check that analytical and FD match
                        //                        bool tCheckLGTestStrain = fem::check( tdLGTestStraindu, tdLGTestStrainduFD, tEpsilon );
                        //                        REQUIRE( tCheckLGTestStrain );

                        //                        // evaluate derivative of the test strain for Euler-Almansi strain by analytical
                        //                        Matrix< DDRMat > tdEATestStraindu =
                        //                                tCMLeader->dTestStraindDOF(
                        //                                        tDofDerivative,
                        //                                        CM_Function_Type::EULERIAN );
                        //                        print( tdEATestStraindu, "tdEATestStraindu" );
                        //
                        //                        // evaluate derivative of the test strain for Euler-Almansi strain by FD
                        //                        Matrix< DDRMat > tdEATestStrainduFD;
                        //                        tCMLeader->eval_derivative_FD(
                        //                                CM_Request_Type::TEST_STRAIN,
                        //                                tdEATestStrainduFD,
                        //                                tDofDerivative,
                        //                                tPerturbation,
                        //                                tDofDerivative,    // dummy
                        //                                tNormal,           // dummy
                        //                                tJump,             // dummy
                        //                                FDScheme_Type::POINT_5,
                        //                                CM_Function_Type::EULERIAN );
                        //                        print( tdEATestStrainduFD, "tdEATestStrainduFD" );
                        //
                        //                        // check that analytical and FD match
                        //                        bool tCheckEATestStrain = fem::check( tdEATestStraindu, tdEATestStrainduFD, tEpsilon );
                        //                        REQUIRE( tCheckEATestStrain );
                    }
                }
            }
            // clean up
            tLeaderFIs.clear();
        }
    }
} /*END_TEST_CASE*/
