/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Struc_Nonlinear_Dirichlet.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private public
// FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
// MTK/src
#include "cl_MTK_Enums.hpp"
// FEM//INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
// LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

void
Test_IWG_Struc_NL_Dirichlet(
        fem::Constitutive_Type aConstitutiveModel,
        IWG_Type               aNitscheFormulation )
{
    // define an epsilon environment
    real tEpsilon = 1E-5;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

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

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = tDispDofTypes;

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderEMod = std::make_shared< fem::Property >();
    tPropLeaderEMod->set_parameters( { { { 10.0 } } } );
    tPropLeaderEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderNu = std::make_shared< fem::Property >();
    tPropLeaderNu->set_parameters( { { { 0.3 } } } );
    tPropLeaderNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderSelect = std::make_shared< fem::Property >();
    tPropLeaderSelect->set_val_function( tMValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( aConstitutiveModel );
    tCMLeaderStrucLinIso->set_dof_type_list( { tDispDofTypes } );
    tCMLeaderStrucLinIso->set_property( tPropLeaderEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropLeaderNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { { { 100.0 } } } );
    tSPDirichletNitsche->set_property( tPropLeaderEMod, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPDirichletNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( aNitscheFormulation );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
    tIWG->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderSelect, "Select" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set selection matrix parameters
        tPropLeaderSelect->set_parameters( { { { static_cast< real >( iSpaceDim ) } } } );

        // Dirichlet BC
        Matrix< DDRMat > tImposedDisp( iSpaceDim, 1, 0.5 );
        tPropLeaderDirichlet->set_parameters( { tImposedDisp } );

        // set space dimension to CM, SP
        if ( iSpaceDim == 2 )
        {
            tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        }
        else
        {
            tCMLeaderStrucLinIso->set_model_type( fem::Model_Type::FULL );
        }
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // create an interpolation order
            mtk::Interpolation_Order tGIInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // space and time geometry interpolators
            //------------------------------------------------------------------------------
            // create a space geometry interpolation rule
            mtk::Interpolation_Rule tGIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tGIInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create a space time geometry interpolator
            Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

            // create time coeff tHat
            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

            Matrix< DDRMat > tXHat;
            fill_xhat_Elast( tXHat, iSpaceDim, iInterpOrder );

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
            int tNumDofDisp = tNumCoeff * iSpaceDim;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );

            tLeaderDofHatDisp = tLeaderDofHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp, tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

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
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}

//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_NL_Dirichlet_Symmetric_Nitsche_Saint_Venant",
        "[moris],[fem],[IWG_Struc_NL_Dirichlet_Symmetric_Nitsche_Saint_Venant]" )
{
    Test_IWG_Struc_NL_Dirichlet(
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF,
            IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_SE );

} /* END_TEST_CASE */

//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_NL_Dirichlet_Unsymmetric_Nitsche_Saint_Venant",
        "[moris],[fem],[IWG_Struc_NL_Dirichlet_Unsymmetric_Nitsche_Saint_Venant]" )
{
    Test_IWG_Struc_NL_Dirichlet(
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF,
            IWG_Type::STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_SE );

} /* END_TEST_CASE */

//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_NL_Dirichlet_Symmetric_Nitsche_Neo_Hookean_SE",
        "[moris],[fem],[IWG_Struc_NL_Dirichlet_Symmetric_Nitsche_Neo_Hookean_SE]" )
{
    Test_IWG_Struc_NL_Dirichlet(
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN,
            IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_SE );

} /* END_TEST_CASE */

//---------------------------------------------------------------------------------------------

TEST_CASE( "IWG_Struc_NL_Dirichlet_Unsymmetric_Nitsche_Neo_Hookean_SE",
        "[moris],[fem],[IWG_Struc_NL_Dirichlet_Unymmetric_Nitsche_Neo_Hookean_SE]" )
{
    Test_IWG_Struc_NL_Dirichlet(
            fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN,
            IWG_Type::STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_SE );

} /* END_TEST_CASE */

//---------------------------------------------------------------------------------------------

// FIXME: PF versions do not work - error in jacobian

// TEST_CASE( "IWG_Struc_NL_Dirichlet_Symmetric_Nitsche_Neo_Hookean_PF",
//         "[moris],[fem],[IWG_Struc_NL_Dirichlet_Symmetric_Nitsche_Neo_Hookean_PF]" )
//{
//     Test_IWG_Struc_NL_Dirichlet(
//             fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN,
//             IWG_Type::STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF );
//
// } /* END_TEST_CASE */

//---------------------------------------------------------------------------------------------

// FIXME: PF versions do not work - error in jacobian

// TEST_CASE( "IWG_Struc_NL_Dirichlet_Unsymmetric_Nitsche_Neo_Hookean_PF",
//         "[moris],[fem],[IWG_Struc_NL_Dirichlet_Unsymmetric_Nitsche_Neo_Hookean_PF]" )
//{
//     Test_IWG_Struc_NL_Dirichlet(
//             fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN,
//             IWG_Type::STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_PF );
//
// } /* END_TEST_CASE */
