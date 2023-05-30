/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Struc_Linear_Axisymmetric_Virtual_Work_Ghost.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM//INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
//LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

// This UT tests the isotropic linear elasticity vitual work based ghost IWG
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order
TEST_CASE( "IWG_Elast_Axi_VWGhost", "[moris],[fem],[axi],[IWG_Elast_Axi_VWGhost]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderYoungModulus = std::make_shared< fem::Property > ();
    tPropLeaderYoungModulus->set_parameters( { {{ 1.0 }} } );
    tPropLeaderYoungModulus->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropLeaderPoissonRatio = std::make_shared< fem::Property > ();
    tPropLeaderPoissonRatio->set_parameters( { {{ 0.3 }} } );
    tPropLeaderPoissonRatio->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropFollowerYoungModulus = std::make_shared< fem::Property > ();
    tPropFollowerYoungModulus->set_parameters( { {{ 2.0 }} } );
    tPropFollowerYoungModulus->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropFollowerPoissonRatio = std::make_shared< fem::Property > ();
    tPropFollowerPoissonRatio->set_parameters( { {{ 0.3 }} } );
    tPropFollowerPoissonRatio->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tRotAxisRadVec = std::make_shared< fem::Property >();
    tRotAxisRadVec->set_parameters( { {{ -4.0, 0}, {0, -4.0}} } );
    tRotAxisRadVec->set_val_function( tAxisymRotAxisFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderElastLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderElastLinIso->set_model_type( fem::Model_Type::AXISYMMETRIC );
    tCMLeaderElastLinIso->set_dof_type_list( tDispDofTypes );
    tCMLeaderElastLinIso->set_property( tPropLeaderYoungModulus, "YoungsModulus" );
    tCMLeaderElastLinIso->set_property( tPropLeaderPoissonRatio, "PoissonRatio" );
    tCMLeaderElastLinIso->set_property( tRotAxisRadVec, "AxisymRotationAxis" );
    tCMLeaderElastLinIso->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMFollowerElastLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMFollowerElastLinIso->set_model_type( fem::Model_Type::AXISYMMETRIC );
    tCMFollowerElastLinIso->set_dof_type_list( tDispDofTypes );
    tCMFollowerElastLinIso->set_property( tPropFollowerYoungModulus, "YoungsModulus" );
    tCMFollowerElastLinIso->set_property( tPropFollowerPoissonRatio, "PoissonRatio" );
    tCMFollowerElastLinIso->set_property( tRotAxisRadVec, "AxisymRotationAxis" );
    tCMFollowerElastLinIso->set_local_properties();

    // define stabilization parameters
      fem::SP_Factory tSPFactory;

      std::shared_ptr< fem::Stabilization_Parameter > tSP1 =
              tSPFactory.create_SP( fem::Stabilization_Type::GHOST_VW );
      tSP1->set_parameters( {{{ 100.0 }} });

      // create a dummy fem cluster and set it to SP
      fem::Cluster * tCluster = new fem::Cluster();
      tSP1->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_VW_GHOST );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( { tDispDofTypes }, mtk::Leader_Follower::LEADER );
    tIWG->set_dof_type_list( { tDispDofTypes }, mtk::Leader_Follower::FOLLOWER );
    tIWG->set_constitutive_model( tCMLeaderElastLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMFollowerElastLinIso,  "ElastLinIso", mtk::Leader_Follower::FOLLOWER );
    tIWG->set_stabilization_parameter( tSP1, "GhostVW" );
    tIWG->set_property( tRotAxisRadVec, "Thickness" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set follower dof type map
    tIWG->mSet->mFollowerDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mFollowerDofTypeMap ( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 3; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set space dim for CM
        tCMLeaderElastLinIso->set_space_dim( iSpaceDim );
        tCMFollowerElastLinIso->set_space_dim( iSpaceDim );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch( iSpaceDim )
        {
            case 2 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD only." );
                break;
            }
        }

        // FIXME higher order
        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 2; iInterpOrder++ )
        {
            // set interpolation order for ghost
            tIWG->set_interpolation_order( iInterpOrder );

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
            Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

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
            int tNumDofDisp  = tNumCoeff * iSpaceDim;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatDisp;
            fill_uhat_Elast( tLeaderDOFHatDisp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDisp );

            // fill random coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatDisp;
            fill_uhat_Elast( tFollowerDOFHatDisp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tFollowerFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDispDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, 2 * tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofDisp - 1 },
                    { tNumDofDisp, 2 * tNumDofDisp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofDisp, 2 * tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDispDofTypes;
            tIWG->mRequestedFollowerGlobalDofTypes  = tDispDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tLeaderFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tFollowerFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );

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

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mFollowerFIManager->set_space_time( tParamPoint );

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
                if( !tCheckJacobian )
                {
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<" iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }
}/* END_TEST_CASE */

