/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Diffusion_Test.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private public

#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"

#undef protected
#undef private

#include "cl_MTK_Enums.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"

#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Diffusion_UT.cpp"

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

void Test_IWG_Diffusion_Bulk(
        real& aConductivity,
        real& aLoad,
        real& aH2Pen,
        real& aH3Pen,
        real& aPhaseField,
        real& aEpsilon )
{
    // define an epsilon environment
    real tEpsilon = aEpsilon;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

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
    Vector< Vector< MSI::Dof_Type > > tTempDofTypes = { { MSI::Dof_Type::TEMP } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { { { aConductivity } } } );
    tPropConductivity->set_val_function( tConstValFunc_Diff );

    std::shared_ptr< fem::Property > tPropLoad = std::make_shared< fem::Property >();
    tPropLoad->set_parameters( { { { aLoad } } } );
    tPropLoad->set_val_function( tConstValFunc_Diff );

    std::shared_ptr< fem::Property > tPropH2 = std::make_shared< fem::Property >();
    tPropH2->set_parameters( { { { aH2Pen } } } );

    std::shared_ptr< fem::Property > tPropH3 = std::make_shared< fem::Property >();
    tPropH3->set_parameters( { { { aH3Pen } } } );

    std::shared_ptr< fem::Property > tPropPhaseField = std::make_shared< fem::Property >();
    tPropPhaseField->set_parameters( { { { aPhaseField } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiffusion =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMLeaderDiffusion->set_dof_type_list( { tTempDofTypes } );
    tCMLeaderDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMLeaderDiffusion->set_local_properties();

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
    tIWG->set_residual_dof_type( tTempDofTypes );
    tIWG->set_dof_type_list( tTempDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMLeaderDiffusion, "Diffusion" );
    tIWG->set_property( tPropLoad, "Load", mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropH2, "H2Penalty", mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropH3, "H3Penalty", mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropPhaseField, "PhaseField", mtk::Leader_Follower::LEADER );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

    // initialize error flag
    bool tErrorFlag = false;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        tCMLeaderDiffusion->set_space_dim( iSpaceDim );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = {
                    { 0.0, 0.0 },
                    { 0.9, 0.0 },
                    { 1.0, 1.1 },
                    { 0.0, 0.8 }
                };
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = {
                    { 0.0, 0.0, 0.0 },
                    { 0.9, 0.0, 0.0 },
                    { 1.0, 1.1, 0.0 },
                    { 0.0, 0.8, 0.0 },
                    { 0.2, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 0.9, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 }
                };
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // set space dimension to CM, SP
        tCMLeaderDiffusion->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
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
            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

            // set the coefficients xHat, tHat
            tGI.set_coeff( tXHat, tTHat );

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
            int tNumDofTemp = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDofHatTemp;
            fill_that( tLeaderDofHatTemp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tTempDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatTemp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tTempDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofTemp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = { { 0, tNumDofTemp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tTempDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofTemp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofTemp, tNumDofTemp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tTempDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >   tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager               tFIManager( tTempDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
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
                if ( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << '\n';

                    // set error flag to true
                    tErrorFlag = true;
                }
            }

            // clean up
            tLeaderFIs.clear();
        }

        // require error flag is false
        REQUIRE( !tErrorFlag );
    }
}

//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Diffusion_Bulk - Conductivity Only", "[moris],[fem],[IWG_Elasticity_Bulk]" )
{
    real tConductivity = 5.0;
    real tLoad         = 2.5;
    real tH2Pen        = 0.0;
    real tH3Pen        = 0.0;
    real tPhaseField   = 0.0;
    real tEpsilon      = 1e-6;

    Test_IWG_Diffusion_Bulk( tConductivity, tLoad, tH2Pen, tH3Pen, tPhaseField, tEpsilon );
}

//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Diffusion_Bulk - H2 Penalty Only", "[moris],[fem],[IWG_Elasticity_Bulk]" )
{
    real tConductivity = 0.0;
    real tLoad         = 0.0;
    real tH2Pen        = 1.0;
    real tH3Pen        = 0.0;
    real tPhaseField   = 0.0;
    real tEpsilon      = 2e-6;

    Test_IWG_Diffusion_Bulk( tConductivity, tLoad, tH2Pen, tH3Pen, tPhaseField, tEpsilon );
}

//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Diffusion_Bulk - H3 Penalty Only", "[moris],[fem],[IWG_Elasticity_Bulk]" )
{
    real tConductivity = 0.0;
    real tLoad         = 0.0;
    real tH2Pen        = 0.0;
    real tH3Pen        = 1.0;
    real tPhaseField   = 0.0;
    real tEpsilon      = 2e-5;

    Test_IWG_Diffusion_Bulk( tConductivity, tLoad, tH2Pen, tH3Pen, tPhaseField, tEpsilon );
}
//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Diffusion_Bulk - PhaseField Only", "[moris],[fem],[IWG_Elasticity_Bulk]" )
{
    real tConductivity = 0.0;
    real tLoad         = 0.0;
    real tH2Pen        = 0.0;
    real tH3Pen        = 0.0;
    real tPhaseField   = 1.0;
    real tEpsilon      = 2e-5;

    Test_IWG_Diffusion_Bulk( tConductivity, tLoad, tH2Pen, tH3Pen, tPhaseField, tEpsilon );
}

//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Diffusion_Bulk - All Properties", "[moris],[fem],[IWG_Elasticity_Bulk]" )
{
    real tConductivity = 1.0;
    real tLoad         = 0.0;
    real tH2Pen        = 1.0;
    real tH3Pen        = 1.0;
    real tPhaseField   = 1.0;
    real tEpsilon      = 2e-5;

    Test_IWG_Diffusion_Bulk( tConductivity, tLoad, tH2Pen, tH3Pen, tPhaseField, tEpsilon );
}

/*END_TEST_CASE*/
