/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Diffusion_Ghost.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
#include "cl_MTK_Enums.hpp"
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Ghost_UT.cpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

/**
 * @brief This UT tests the normal field ghost IWG for QUAD, HEX geometry type
 * for LINEAR, QUADRATIC and CUBIC interpolation order against an analytic solution
 */
TEST_CASE( "IWG_Diff_Ghost_Analytic", "[moris],[fem],[IWG_Diff_Ghost_Analytic]" )
{
    // define a tolerance
    real tEpsilon = 1E-13;

    //------------------------------------------------------------------------------
    // initialize the geometry/element

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
    Matrix< DDUMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< Vector< MSI::Dof_Type > > tDofTypes  = { { MSI::Dof_Type::TEMP } };

    //------------------------------------------------------------------------------
    // initialize the IWG
    
    // create a unit property for the unit stabilization parameter
    std::shared_ptr< fem::Property > tPropUnit = std::make_shared< fem::Property > ();
    tPropUnit->set_parameters( { {{ 1.0 }} } );
    tPropUnit->set_val_function( tConstValFunc );

    // define a unit ghost stabilization parameter
    fem::SP_Factory tSPFactory;
    std::shared_ptr< fem::Stabilization_Parameter > tSP = tSPFactory.create_SP( fem::Stabilization_Type::GHOST_DISPL );
    tSP->set_parameters( { {{ 1.0 }} });
    tSP->set_property( tPropUnit, "Material", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSP->set_cluster( tCluster );

    // define the IWG
    fem::IWG_Factory tIWGFactory;
    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::GHOST_NORMAL_FIELD );
    tIWG->set_stabilization_parameter( tSP, "GhostSP" );
    tIWG->set_residual_dof_type( tDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::FOLLOWER );

    //------------------------------------------------------------------------------
    // initialize set info

    // create a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
// for( uint iSpaceDim = 2; iSpaceDim < 3; iSpaceDim++ )
    {
        //------------------------------------------------------------------------------
        // space and time geometry interpolators

        // get geometry data for the unit element
        fill_X_hat( iSpaceDim, tXHat, tGeometryType );

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

        // loop on the interpolation orders
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // create and set normal
            Matrix< DDRMat > tNormal;
            fill_normal( tNormal, iSpaceDim, iInterpOrder );
            tIWG->set_normal( tNormal );

            //------------------------------------------------------------------------------
            // integration points

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

            // get integration weights
            Matrix< DDRMat > tIntegWeights;
            tIntegrator.get_weights( tIntegWeights );

            //------------------------------------------------------------------------------
            // field interpolators

            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // number of coefficients per element for the current interpolation order
            int tNumCoeff = (int) tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( 
                    tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatTemp;
            fill_u_hat_leader( tLeaderDOFHatTemp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator temperature
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatTemp );

            // fill random coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatTemp;
            fill_u_hat_follower( tFollowerDOFHatTemp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( tDofTypes.size() );

            // create the field interpolator temperature
            tFollowerFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatTemp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumCoeff - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumCoeff, 2 * tNumCoeff - 1 } };

            // assemble the two coefficient vectors into one
            Matrix< DDRMat > tGlobalUHat( 2 * tNumCoeff, 1, 0.0 );
            tGlobalUHat( { 0, tNumCoeff - 1 }, { 0, 0 } ) = tLeaderDOFHatTemp.matrix_data();
            tGlobalUHat( { tNumCoeff, 2 * tNumCoeff - 1 }, { 0, 0 } ) = tFollowerDOFHatTemp.matrix_data();
            Matrix< DDRMat > tTransGlobalUHat = trans( tGlobalUHat );

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = { { 0, tNumCoeff-1 }, { tNumCoeff, 2 * tNumCoeff - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumCoeff, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumCoeff, 2 * tNumCoeff, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > > tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tLeaderFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tFollowerFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tLeaderFIManager.mFI = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator = &tGI;
            tLeaderFIManager.mIGGeometryInterpolator = &tGI;
            tFollowerFIManager.mFI = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator = &tGI;
            tFollowerFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tLeaderFIManager;
            tIWG->mSet->mFollowerFIManager = &tFollowerFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tLeaderFIManager, mtk::Leader_Follower::LEADER );
            tIWG->set_field_interpolator_manager( &tFollowerFIManager, mtk::Leader_Follower::FOLLOWER );

            // set the interpolation order for the ghost sum
            tIWG->set_interpolation_order( iInterpOrder );

            // integrate the residual
            Matrix< DDRMat > tIntegrandFromResidual( 1, 1, 0.0 );
            Matrix< DDRMat > tIntegrandFromJacobian( 1, 1, 0.0 );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                //------------------------------------------------------------------------------
                // reset and set Gauss point

                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mFollowerFIManager->set_space_time( tParamPoint );

                // reset residual and jacobian
                tIWG->mSet->mResidual( 0 ).fill( 0.0 );
                tIWG->mSet->mJacobian.fill( 0.0 );

                //------------------------------------------------------------------------------
                // perform integration

                // compute residual and jacobian
                real tIntegWeight = std::pow( 0.5, (real) iSpaceDim + 1 ) * tIntegWeights( iGP );
                tIWG->compute_residual( tIntegWeight );
                tIWG->compute_jacobian( tIntegWeight );

                // add to integral sum
                tIntegrandFromResidual( { 0, 0 }, { 0, 0 } ) += tTransGlobalUHat * tIWG->mSet->mResidual( 0 );
                tIntegrandFromJacobian( { 0, 0 }, { 0, 0 } ) += tTransGlobalUHat * tIWG->mSet->mJacobian * tGlobalUHat;
            }

            // analytical reference value
            real tRefVal = get_ref_value( iSpaceDim, iInterpOrder );

            // get the check values
            bool tCheckResidual = ( std::abs( ( tIntegrandFromResidual( 0, 0 ) - tRefVal ) / tRefVal ) < tEpsilon );
            bool tCheckJacobian = ( std::abs( ( tIntegrandFromJacobian( 0, 0 ) - tRefVal ) / tRefVal ) < tEpsilon );

            // print information if something goes wrong
            if( !tCheckResidual || !tCheckJacobian )
            {
                std::cout << "\n-------------------------------------------" << std::endl;
                std::cout << "Residual or Jacobian check failed in " << iSpaceDim << "D and for p=" << iInterpOrder << std::endl;
                std::cout << "Reference value is: " << tRefVal << std::endl;
                std::cout << "Residual value is: " << tIntegrandFromResidual( 0, 0 ) << std::endl;
                std::cout << "Residual value (obtained through the Jacobian) is: " << tIntegrandFromJacobian( 0, 0 ) << std::endl;
                std::cout << "Absolute differences are: " << std::abs( tIntegrandFromResidual( 0, 0 ) - tRefVal ) << 
                        " and " << std::abs( tIntegrandFromResidual( 0, 0 ) - tRefVal ) << std::endl;
                std::cout << "Relative differences are: " << std::abs( ( tIntegrandFromResidual( 0, 0 ) - tRefVal ) / tRefVal ) << 
                        " and " << std::abs( ( tIntegrandFromJacobian( 0, 0 ) - tRefVal ) / tRefVal ) << std::endl;
                std::cout << "-------------------------------------------\n" << std::endl;
            }

            // check that the integral sums are as expected
            REQUIRE( tCheckResidual );
            REQUIRE( tCheckJacobian );

            // clean up before going to next Element type
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
}/* END_TEST_CASE */

