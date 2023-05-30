/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#define protected public
#define private public
// FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp"
#undef protected
#undef private
// LINALG/src
#include "op_equal_equal.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
// FEM//INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

using namespace moris;
using namespace fem;

inline void
tWallDistanceValFunc(
        moris::Matrix< moris::DDRMat >&                aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    moris::fem::Field_Interpolator* tFIWallDist =
            aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::L2 );

    aPropMatrix = aParameters( 0 ) * tFIWallDist->val();
}

inline void
tWallDistanceDerFunc(
        moris::Matrix< moris::DDRMat >&                aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    moris::fem::Field_Interpolator* tFIWallDist =
            aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::L2 );

    aPropMatrix = aParameters( 0 ) * tFIWallDist->N();
}

TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Bulk", "[IWG_Spalart_Allmaras_Turbulence_Bulk]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

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

    // dof type list
    moris::Cell< MSI::Dof_Type > tVelDofTypes = { MSI::Dof_Type::VX };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVisDofTypes = { { MSI::Dof_Type::VISCOSITY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes    = { tVelDofTypes, tVisDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropWallDistance = std::make_shared< fem::Property >();
    tPropWallDistance->set_parameters( { { { 1.0 } } } );
    tPropWallDistance->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    // tPropViscosity->set_parameters( { {{ 2.0 }} } );
    tPropViscosity->set_val_function( tConstValFunc );
    tPropViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );
    // tPropViscosity->set_dof_type_list( { tVisDofTypes } );
    // tPropViscosity->set_val_function( tVISCOSITYFIValFunc );
    // tPropViscosity->set_dof_derivative_functions( { tVISCOSITYFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMLeaderSATurbulence->set_dof_type_list( tDofTypes );
    tCMLeaderSATurbulence->set_property( tPropWallDistance, "WallDistance" );
    tCMLeaderSATurbulence->set_property( tPropViscosity, "KinViscosity" );
    tCMLeaderSATurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPSUPG =
            tSPFactory.create_SP( fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
    tSPSUPG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tSPSUPG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence" );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPSUPG->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
    tIWG->set_residual_dof_type( tVisDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence" );
    tIWG->set_stabilization_parameter( tSPSUPG, "SUPG" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

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

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };

                // set velocity dof types
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

                // set viscosity property parameters
                tPropViscosity->set_parameters( { { { 2.0 } }, { { 0.0 }, { 0.0 } } } );
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
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };

                // set viscosity property parameters
                tPropViscosity->set_parameters( { { { 2.0 } }, { { 0.0 }, { 0.0 }, { 0.0 } } } );
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        Matrix< DDRMat > tGravity( iSpaceDim, 1, 10.0 );

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

        // set space dimension to CM, SP
        tCMLeaderSATurbulence->set_space_dim( iSpaceDim );
        tSPSUPG->set_space_dim( iSpaceDim );

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
            int tNumDofVel = tNumCoeff * iSpaceDim;
            int tNumDofVis = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // loop over different state variable configurations
            for ( uint iCu = 0; iCu < 4; ++iCu )
            {
                for ( uint iCp = 0; iCp < 3; ++iCp )
                {
                    // fill coefficients for leader FI
                    Matrix< DDRMat > tLeaderDOFHatVel;
                    fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder, iCu );
                    Matrix< DDRMat > tLeaderDOFHatVis;
                    fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder, iCp );

                    // create a cell of field interpolators for IWG
                    Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

                    // create the field interpolator velocity
                    tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

                    // create the field interpolator viscosity
                    tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVis );

                    // set size and fill the set residual assembly map
                    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofVis - 1 } };

                    // set size and fill the set jacobian assembly map
                    Matrix< DDSMat > tJacAssembly = {
                        { 0, tNumDofVel - 1 },
                        { tNumDofVel, tNumDofVel + tNumDofVis - 1 }
                    };
                    tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
                    tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

                    // set size and init the set residual and jacobian
                    tIWG->mSet->mResidual.resize( 1 );
                    tIWG->mSet->mResidual( 0 ).set_size( tNumDofVel + tNumDofVis, 1, 0.0 );
                    tIWG->mSet->mJacobian.set_size( tNumDofVel + tNumDofVis, tNumDofVel + tNumDofVis, 0.0 );

                    // build global dof type list
                    tIWG->get_global_dof_type_list();

                    // populate the requested leader dof type
                    tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

                    // create a field interpolator manager
                    moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
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

                    // loop over integration points
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
                            std::cout << "Case: Geometry " << iSpaceDim
                                      << " Order " << iInterpOrder
                                      << " VelConf " << iCu
                                      << " PresConf " << iCp
                                      << " iGP " << iGP << std::endl;
                        }

                        // require check is true
                        REQUIRE( tCheckJacobian );
                    }

                    // clean up
                    tLeaderFIs.clear();
                }
            }
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Bulk_Small_Wall_Distance",
        "[IWG_Spalart_Allmaras_Turbulence_Bulk_Small_Wall_Distance]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-5;

    // define turbulent viscosity scaling (should be similar to wall distance, can be positive or negative)
    real tViscosityScaling = 1e-5;

    // define wall distance scaling (needs to be positive and small)
    real tWallDistScaling = 1e-5;

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

    // dof type list
    moris::Cell< MSI::Dof_Type >                tVelDofTypes      = { MSI::Dof_Type::VX };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVisDofTypes      = { { MSI::Dof_Type::VISCOSITY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tWallDistDofTypes = { { MSI::Dof_Type::L2 } };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = { tVelDofTypes, tVisDofTypes( 0 ), tWallDistDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropWallDistance = std::make_shared< fem::Property >();
    tPropWallDistance->set_parameters( { { { 1.0 } } } );
    tPropWallDistance->set_dof_type_list( tWallDistDofTypes );
    tPropWallDistance->set_val_function( tWallDistanceValFunc );
    tPropWallDistance->set_dof_derivative_functions( { tWallDistanceDerFunc } );

    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_val_function( tConstValFunc );
    // tPropViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMLeaderSATurbulence->set_dof_type_list( tDofTypes );
    tCMLeaderSATurbulence->set_property( tPropWallDistance, "WallDistance" );
    tCMLeaderSATurbulence->set_property( tPropViscosity, "KinViscosity" );
    tCMLeaderSATurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPSUPG =
            tSPFactory.create_SP( fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
    tSPSUPG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tSPSUPG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence" );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPSUPG->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
    tIWG->set_residual_dof_type( tVisDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence" );
    tIWG->set_stabilization_parameter( tSPSUPG, "SUPG" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::L2 ) )        = 2;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::L2 ) )        = 2;

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

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };

                // set velocity dof types
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

                // set viscosity property parameters
                tPropViscosity->set_parameters( { { { 2.0 } }, { { 0.0 }, { 0.0 } } } );
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
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };

                // set viscosity property parameters
                tPropViscosity->set_parameters( { { { 2.0 } }, { { 0.0 }, { 0.0 }, { 0.0 } } } );
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        Matrix< DDRMat > tGravity( iSpaceDim, 1, 10.0 );

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

        // set space dimension to CM, SP
        tCMLeaderSATurbulence->set_space_dim( iSpaceDim );
        tSPSUPG->set_space_dim( iSpaceDim );

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
            int tNumDofVel      = tNumCoeff * iSpaceDim;
            int tNumDofVis      = tNumCoeff;
            int tNumDofWallDist = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // loop over different state variable configurations
            for ( uint iCu = 0; iCu < 4; ++iCu )
            {
                for ( uint iCp = 0; iCp < 3; ++iCp )
                {
                    // fill coefficients for leader FI
                    Matrix< DDRMat > tLeaderDOFHatVel;
                    fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );

                    Matrix< DDRMat > tLeaderDOFHatVis;
                    fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );
                    tLeaderDOFHatVis = tViscosityScaling * tLeaderDOFHatVis;

                    Matrix< DDRMat > tLeaderDOFHatWallDist;
                    fill_phat( tLeaderDOFHatWallDist, iSpaceDim, iInterpOrder );
                    tLeaderDOFHatWallDist = tWallDistScaling * tLeaderDOFHatWallDist;

                    // create a cell of field interpolators for IWG
                    Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

                    // create the field interpolator velocity
                    tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

                    // create the field interpolator viscosity
                    tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVis );

                    // create the field interpolator wall distance
                    tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tWallDistDofTypes( 0 ) );
                    tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatWallDist );

                    // set size and fill the set residual assembly map
                    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofVis - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofVel + tNumDofVis, tNumDofVel + tNumDofVis + tNumDofWallDist - 1 } };

                    // set size and fill the set jacobian assembly map
                    Matrix< DDSMat > tJacAssembly = {
                        { 0, tNumDofVel - 1 },
                        { tNumDofVel, tNumDofVel + tNumDofVis - 1 },
                        { tNumDofVel + tNumDofVis, tNumDofVel + tNumDofVis + tNumDofWallDist - 1 }
                    };

                    tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
                    tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
                    tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;

                    // set size and init the set residual and jacobian
                    tIWG->mSet->mResidual.resize( 1 );
                    tIWG->mSet->mResidual( 0 ).set_size( tNumDofVel + tNumDofVis + tNumDofWallDist, 1, 0.0 );
                    tIWG->mSet->mJacobian.set_size( tNumDofVel + tNumDofVis + tNumDofWallDist, tNumDofVel + tNumDofVis + tNumDofWallDist, 0.0 );

                    // build global dof type list
                    tIWG->get_global_dof_type_list();

                    // populate the requested leader dof type
                    tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

                    // create a field interpolator manager
                    moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
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

                    // loop over integration points
                    uint tNumGPs = tIntegPoints.n_cols();
                    for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
                    {
                        // reset IWG evaluation flags
                        tIWG->reset_eval_flags();

                        // create evaluation point xi, tau
                        Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                        // set integration point
                        tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                        //-------------------------------------------------------------------------------------------------------------------------
                        // check components of constitutive model - derivatives with respect to velocity dofs
                        {
                            real tErrorLimit = 1e-5;

                            moris::fem::IWG_Spalart_Allmaras_Turbulence_Bulk* tSAIWG =
                                    dynamic_cast< moris::fem::IWG_Spalart_Allmaras_Turbulence_Bulk* >( tIWG.get() );

                            const Matrix< DDRMat >& tJacDestructionCoefficient = tCMLeaderSATurbulence->dwalldestructioncoeffdu( tVelDofTypes );
                            const Matrix< DDRMat >& tJacProductionCoefficient  = tCMLeaderSATurbulence->dproductioncoeffdu( tVelDofTypes );
                            const Matrix< DDRMat >& tJacDiffusionCoefficient   = tCMLeaderSATurbulence->ddiffusioncoeffdu( tVelDofTypes );
                            const Matrix< DDRMat >& tJacDivFlux                = tCMLeaderSATurbulence->ddivfluxdu( tVelDofTypes );
                            const Matrix< DDRMat >& tJacSUPG                   = tSPSUPG->dSPdLeaderDOF( tVelDofTypes );
                            const Matrix< DDRMat >& tJacDestructionTerm        = tCMLeaderSATurbulence->dwalldestructiontermdu( tVelDofTypes );
                            const Matrix< DDRMat >& tJacProductionTerm         = tCMLeaderSATurbulence->dproductiontermdu( tVelDofTypes );

                            Matrix< DDRMat > tJacStrongForm, tResStrongForm;
                            tSAIWG->compute_jacobian_strong_form( tVelDofTypes, tJacStrongForm );

                            Matrix< DDRMat > tJacDestructionCoefficientFD( 1, tLeaderDOFHatVel.numel(), 0.0 );
                            Matrix< DDRMat > tJacProductionCoefficientFD( 1, tLeaderDOFHatVel.numel(), 0.0 );
                            Matrix< DDRMat > tJacDiffusionCoefficientFD( 1, tLeaderDOFHatVel.numel(), 0.0 );
                            Matrix< DDRMat > tJacDivFluxFD( 1, tLeaderDOFHatVel.numel(), 0.0 );
                            Matrix< DDRMat > tJacSUPGFD( 1, tLeaderDOFHatVel.numel(), 0.0 );
                            Matrix< DDRMat > tJacStrongFormFD( 1, tLeaderDOFHatVel.numel(), 0.0 );
                            Matrix< DDRMat > tJacDestructionTermFD( 1, tLeaderDOFHatVel.numel(), 0.0 );
                            Matrix< DDRMat > tJacProductionTermFD( 1, tLeaderDOFHatVel.numel(), 0.0 );

                            uint tCount = 0;

                            for ( uint iDim = 0; iDim < tLeaderDOFHatVel.n_cols(); ++iDim )
                            {
                                for ( uint iCoef = 0; iCoef < tLeaderDOFHatVel.n_rows(); ++iCoef )
                                {
                                    tLeaderDOFHatVel( iCoef, iDim ) += tPerturbation;
                                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );
                                    tLeaderFIs( 0 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();
                                    tJacDestructionCoefficientFD( tCount ) = tCMLeaderSATurbulence->wall_destruction_coefficient()( 0 );
                                    tJacProductionCoefficientFD( tCount )  = tCMLeaderSATurbulence->production_coefficient()( 0 );
                                    tJacDestructionTermFD( tCount )        = tCMLeaderSATurbulence->wall_destruction_term()( 0 );
                                    tJacProductionTermFD( tCount )         = tCMLeaderSATurbulence->production_term()( 0 );
                                    tJacDiffusionCoefficientFD( tCount )   = tCMLeaderSATurbulence->diffusion_coefficient()( 0 );
                                    tJacDivFluxFD( tCount )                = tCMLeaderSATurbulence->divflux()( 0 );
                                    tJacSUPGFD( tCount )                   = tSPSUPG->val()( 0 );
                                    tSAIWG->compute_residual_strong_form( tResStrongForm );
                                    tJacStrongFormFD( tCount ) = tResStrongForm( 0 );

                                    tLeaderDOFHatVel( iCoef, iDim ) -= 2.0 * tPerturbation;
                                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );
                                    tLeaderFIs( 0 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();
                                    tJacDestructionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->wall_destruction_coefficient()( 0 );
                                    tJacProductionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->production_coefficient()( 0 );
                                    tJacDestructionTermFD( tCount ) -= tCMLeaderSATurbulence->wall_destruction_term()( 0 );
                                    tJacProductionTermFD( tCount ) -= tCMLeaderSATurbulence->production_term()( 0 );
                                    tJacDiffusionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->diffusion_coefficient()( 0 );
                                    tJacDivFluxFD( tCount ) -= tCMLeaderSATurbulence->divflux()( 0 );
                                    tJacSUPGFD( tCount ) -= tSPSUPG->val()( 0 );
                                    tSAIWG->compute_residual_strong_form( tResStrongForm );
                                    tJacStrongFormFD( tCount ) -= tResStrongForm( 0 );

                                    tJacDestructionCoefficientFD( tCount ) = tJacDestructionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacProductionCoefficientFD( tCount )  = tJacProductionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacDestructionTermFD( tCount )        = tJacDestructionTermFD( tCount ) / 2.0 / tPerturbation;
                                    tJacProductionTermFD( tCount )         = tJacProductionTermFD( tCount ) / 2.0 / tPerturbation;
                                    tJacDiffusionCoefficientFD( tCount )   = tJacDiffusionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacDivFluxFD( tCount )                = tJacDivFluxFD( tCount ) / 2.0 / tPerturbation;
                                    tJacSUPGFD( tCount )                   = tJacSUPGFD( tCount ) / 2.0 / tPerturbation;
                                    tJacStrongFormFD( tCount )             = tJacStrongFormFD( tCount ) / 2.0 / tPerturbation;

                                    tLeaderDOFHatVel( iCoef, iDim ) += tPerturbation;
                                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );
                                    tLeaderFIs( 0 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();

                                    tCount++;
                                }
                            }

                            real tRelError;
                            for ( uint iI = 0; iI < tCount; ++iI )
                            {
                                tRelError = std::abs( tJacDestructionCoefficient( iI ) - tJacDestructionCoefficientFD( iI ) ) / ( std::abs( tJacDestructionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDestructionCoefficient velocity dof: ana = "
                                              << tJacDestructionCoefficient( iI ) << " FD = "
                                              << tJacDestructionCoefficientFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacDestructionTerm( iI ) - tJacDestructionTermFD( iI ) ) / ( std::abs( tJacDestructionTermFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDestructionTerm velocity dof: ana = "
                                              << tJacDestructionTerm( iI ) << " FD = "
                                              << tJacDestructionTermFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacProductionCoefficient( iI ) - tJacProductionCoefficientFD( iI ) ) / ( std::abs( tJacProductionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacProductionCoefficient velocity dof: ana = "
                                              << tJacProductionCoefficient( iI ) << " FD = "
                                              << tJacProductionCoefficientFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacProductionTerm( iI ) - tJacProductionTermFD( iI ) ) / ( std::abs( tJacProductionTermFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacProductionTerm velocity dof: ana = "
                                              << tJacProductionTerm( iI ) << " FD = "
                                              << tJacProductionTermFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacDiffusionCoefficient( iI ) - tJacDiffusionCoefficientFD( iI ) ) / ( std::abs( tJacDiffusionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDiffusionCoefficientFD velocity dof:"
                                              << " ana = " << tJacDiffusionCoefficient( iI )
                                              << " FD  = " << tJacDiffusionCoefficientFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }

                                tRelError = std::abs( tJacDivFlux( iI ) - tJacDivFluxFD( iI ) ) / ( std::abs( tJacDivFluxFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDivFluxFD velocity dof:"
                                              << " ana = " << tJacDivFlux( iI )
                                              << " FD  = " << tJacDivFluxFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }

                                tRelError = std::abs( tJacSUPG( iI ) - tJacSUPGFD( iI ) ) / ( std::abs( tJacSUPGFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacSUPGFD velocity dof:"
                                              << " ana = " << tJacSUPG( iI )
                                              << " FD  = " << tJacSUPGFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }

                                tRelError = std::abs( tJacStrongForm( iI ) - tJacStrongFormFD( iI ) ) / ( std::abs( tJacStrongFormFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacStrongForm velocity dof:"
                                              << " ana = " << tJacStrongForm( iI )
                                              << " FD  = " << tJacStrongFormFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }
                            }
                        }

                        //-------------------------------------------------------------------------------------------------------------------------
                        // check components of constitutive model - derivatives with respect to viscosity dofs
                        {
                            real tErrorLimit = 1e-5;

                            moris::fem::IWG_Spalart_Allmaras_Turbulence_Bulk* tSAIWG =
                                    dynamic_cast< moris::fem::IWG_Spalart_Allmaras_Turbulence_Bulk* >( tIWG.get() );

                            const Matrix< DDRMat >& tJacDestructionCoefficient = tCMLeaderSATurbulence->dwalldestructioncoeffdu( tVisDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacProductionCoefficient  = tCMLeaderSATurbulence->dproductioncoeffdu( tVisDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacDiffusionCoefficient   = tCMLeaderSATurbulence->ddiffusioncoeffdu( tVisDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacDivFlux                = tCMLeaderSATurbulence->ddivfluxdu( tVisDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacSUPG                   = tSPSUPG->dSPdLeaderDOF( tVisDofTypes( 0 ) );

                            Matrix< DDRMat > tJacStrongForm, tResStrongForm;
                            tSAIWG->compute_jacobian_strong_form( tVisDofTypes( 0 ), tJacStrongForm );

                            Matrix< DDRMat > tJacDestructionCoefficientFD( 1, tLeaderDOFHatVis.numel(), 0.0 );
                            Matrix< DDRMat > tJacProductionCoefficientFD( 1, tLeaderDOFHatVis.numel(), 0.0 );
                            Matrix< DDRMat > tJacDiffusionCoefficientFD( 1, tLeaderDOFHatVis.numel(), 0.0 );
                            Matrix< DDRMat > tJacDivFluxFD( 1, tLeaderDOFHatVis.numel(), 0.0 );
                            Matrix< DDRMat > tJacSUPGFD( 1, tLeaderDOFHatVis.numel(), 0.0 );
                            Matrix< DDRMat > tJacStrongFormFD( 1, tLeaderDOFHatVel.numel(), 0.0 );

                            uint tCount = 0;

                            for ( uint iDim = 0; iDim < tLeaderDOFHatVis.n_cols(); ++iDim )
                            {
                                for ( uint iCoef = 0; iCoef < tLeaderDOFHatVis.n_rows(); ++iCoef )
                                {
                                    tLeaderDOFHatVis( iCoef, iDim ) += tPerturbation;
                                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVis );
                                    tLeaderFIs( 1 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();
                                    tJacDestructionCoefficientFD( tCount ) = tCMLeaderSATurbulence->wall_destruction_coefficient()( 0 );
                                    tJacProductionCoefficientFD( tCount )  = tCMLeaderSATurbulence->production_coefficient()( 0 );
                                    tJacDiffusionCoefficientFD( tCount )   = tCMLeaderSATurbulence->diffusion_coefficient()( 0 );
                                    tJacDivFluxFD( tCount )                = tCMLeaderSATurbulence->divflux()( 0 );
                                    tJacSUPGFD( tCount )                   = tSPSUPG->val()( 0 );
                                    tSAIWG->compute_residual_strong_form( tResStrongForm );
                                    tJacStrongFormFD( tCount ) = tResStrongForm( 0 );

                                    tLeaderDOFHatVis( iCoef, iDim ) -= 2.0 * tPerturbation;
                                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVis );
                                    tLeaderFIs( 1 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();
                                    tJacDestructionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->wall_destruction_coefficient()( 0 );
                                    tJacProductionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->production_coefficient()( 0 );
                                    tJacDiffusionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->diffusion_coefficient()( 0 );
                                    tJacDivFluxFD( tCount ) -= tCMLeaderSATurbulence->divflux()( 0 );
                                    tJacSUPGFD( tCount ) -= tSPSUPG->val()( 0 );
                                    tSAIWG->compute_residual_strong_form( tResStrongForm );
                                    tJacStrongFormFD( tCount ) -= tResStrongForm( 0 );

                                    tJacDestructionCoefficientFD( tCount ) = tJacDestructionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacProductionCoefficientFD( tCount )  = tJacProductionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacDiffusionCoefficientFD( tCount )   = tJacDiffusionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacDivFluxFD( tCount )                = tJacDivFluxFD( tCount ) / 2.0 / tPerturbation;
                                    tJacSUPGFD( tCount )                   = tJacSUPGFD( tCount ) / 2.0 / tPerturbation;
                                    tJacStrongFormFD( tCount )             = tJacStrongFormFD( tCount ) / 2.0 / tPerturbation;

                                    tLeaderDOFHatVis( iCoef, iDim ) += tPerturbation;
                                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVis );
                                    tLeaderFIs( 1 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();

                                    tCount++;
                                }
                            }

                            real tRelError;
                            for ( uint iI = 0; iI < tCount; ++iI )
                            {
                                tRelError = std::abs( tJacDestructionCoefficient( iI ) - tJacDestructionCoefficientFD( iI ) ) / ( std::abs( tJacDestructionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDestructionCoefficient viscosity dof: ana = "
                                              << tJacDestructionCoefficient( iI ) << " FD = "
                                              << tJacDestructionCoefficientFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacProductionCoefficient( iI ) - tJacProductionCoefficientFD( iI ) ) / ( std::abs( tJacProductionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacProductionCoefficient viscosity dof: ana = "
                                              << tJacProductionCoefficient( iI ) << " FD = "
                                              << tJacProductionCoefficientFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacDiffusionCoefficient( iI ) - tJacDiffusionCoefficientFD( iI ) ) / ( std::abs( tJacDiffusionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDiffusionCoefficientFD viscosity dof:"
                                              << " ana = " << tJacDiffusionCoefficient( iI )
                                              << " FD  = " << tJacDiffusionCoefficientFD( iI )
                                              << " rel error = " << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacDivFlux( iI ) - tJacDivFluxFD( iI ) ) / ( std::abs( tJacDivFluxFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDivFluxFD viscosity dof:"
                                              << " ana = " << tJacDivFlux( iI )
                                              << " FD  = " << tJacDivFluxFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }

                                tRelError = std::abs( tJacSUPG( iI ) - tJacSUPGFD( iI ) ) / ( std::abs( tJacSUPGFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacSUPGFD viscosity dof:"
                                              << " ana = " << tJacSUPG( iI )
                                              << " FD  = " << tJacSUPGFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }

                                tRelError = std::abs( tJacStrongForm( iI ) - tJacStrongFormFD( iI ) ) / ( std::abs( tJacStrongFormFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacStrongForm viscosity dof:"
                                              << " ana = " << tJacStrongForm( iI )
                                              << " FD  = " << tJacStrongFormFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }
                            }
                        }

                        //-------------------------------------------------------------------------------------------------------------------------
                        // check components of constitutive model - derivatives with respect to wall distance dofs
                        {
                            real tErrorLimit = 1e-5;

                            moris::fem::IWG_Spalart_Allmaras_Turbulence_Bulk* tSAIWG =
                                    dynamic_cast< moris::fem::IWG_Spalart_Allmaras_Turbulence_Bulk* >( tIWG.get() );

                            const Matrix< DDRMat >& tJacDestructionCoefficient = tCMLeaderSATurbulence->dwalldestructioncoeffdu( tWallDistDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacProductionCoefficient  = tCMLeaderSATurbulence->dproductioncoeffdu( tWallDistDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacDiffusionCoefficient   = tCMLeaderSATurbulence->ddiffusioncoeffdu( tWallDistDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacDestructionTerm        = tCMLeaderSATurbulence->dwalldestructiontermdu( tWallDistDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacProductionTerm         = tCMLeaderSATurbulence->dproductiontermdu( tWallDistDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacDivFlux                = tCMLeaderSATurbulence->ddivfluxdu( tWallDistDofTypes( 0 ) );
                            const Matrix< DDRMat >& tJacSUPG                   = tSPSUPG->dSPdLeaderDOF( tWallDistDofTypes( 0 ) );

                            Matrix< DDRMat > tJacStrongForm, tResStrongForm;
                            tSAIWG->compute_jacobian_strong_form( tWallDistDofTypes( 0 ), tJacStrongForm );

                            Matrix< DDRMat > tJacDestructionCoefficientFD( 1, tLeaderDOFHatWallDist.numel(), 0.0 );
                            Matrix< DDRMat > tJacProductionCoefficientFD( 1, tLeaderDOFHatWallDist.numel(), 0.0 );
                            Matrix< DDRMat > tJacDiffusionCoefficientFD( 1, tLeaderDOFHatWallDist.numel(), 0.0 );
                            Matrix< DDRMat > tJacDestructionTermFD( 1, tLeaderDOFHatWallDist.numel(), 0.0 );
                            Matrix< DDRMat > tJacProductionTermFD( 1, tLeaderDOFHatWallDist.numel(), 0.0 );
                            Matrix< DDRMat > tJacDivFluxFD( 1, tLeaderDOFHatWallDist.numel(), 0.0 );
                            Matrix< DDRMat > tJacSUPGFD( 1, tLeaderDOFHatWallDist.numel(), 0.0 );
                            Matrix< DDRMat > tJacStrongFormFD( 1, tLeaderDOFHatWallDist.numel(), 0.0 );

                            uint tCount = 0;

                            for ( uint iDim = 0; iDim < tLeaderDOFHatWallDist.n_cols(); ++iDim )
                            {
                                for ( uint iCoef = 0; iCoef < tLeaderDOFHatWallDist.n_rows(); ++iCoef )
                                {
                                    tLeaderDOFHatWallDist( iCoef, iDim ) += tPerturbation;
                                    tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatWallDist );
                                    tLeaderFIs( 2 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();
                                    tJacDestructionCoefficientFD( tCount ) = tCMLeaderSATurbulence->wall_destruction_coefficient()( 0 );
                                    tJacProductionCoefficientFD( tCount )  = tCMLeaderSATurbulence->production_coefficient()( 0 );
                                    tJacDiffusionCoefficientFD( tCount )   = tCMLeaderSATurbulence->diffusion_coefficient()( 0 );
                                    tJacDestructionTermFD( tCount )        = tCMLeaderSATurbulence->wall_destruction_term()( 0 );
                                    tJacProductionTermFD( tCount )         = tCMLeaderSATurbulence->production_term()( 0 );
                                    tJacDivFluxFD( tCount )                = tCMLeaderSATurbulence->divflux()( 0 );
                                    tJacSUPGFD( tCount )                   = tSPSUPG->val()( 0 );
                                    tSAIWG->compute_residual_strong_form( tResStrongForm );
                                    tJacStrongFormFD( tCount ) = tResStrongForm( 0 );

                                    tLeaderDOFHatWallDist( iCoef, iDim ) -= 2.0 * tPerturbation;
                                    tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatWallDist );
                                    tLeaderFIs( 2 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();
                                    tJacDestructionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->wall_destruction_coefficient()( 0 );
                                    tJacProductionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->production_coefficient()( 0 );
                                    tJacDiffusionCoefficientFD( tCount ) -= tCMLeaderSATurbulence->diffusion_coefficient()( 0 );
                                    tJacDestructionTermFD( tCount ) -= tCMLeaderSATurbulence->wall_destruction_term()( 0 );
                                    tJacProductionTermFD( tCount ) -= tCMLeaderSATurbulence->production_term()( 0 );
                                    tJacDivFluxFD( tCount ) -= tCMLeaderSATurbulence->divflux()( 0 );
                                    tJacSUPGFD( tCount ) -= tSPSUPG->val()( 0 );
                                    tSAIWG->compute_residual_strong_form( tResStrongForm );
                                    tJacStrongFormFD( tCount ) -= tResStrongForm( 0 );

                                    tJacDestructionCoefficientFD( tCount ) = tJacDestructionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacProductionCoefficientFD( tCount )  = tJacProductionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacDiffusionCoefficientFD( tCount )   = tJacDiffusionCoefficientFD( tCount ) / 2.0 / tPerturbation;
                                    tJacDestructionTermFD( tCount )        = tJacDestructionTermFD( tCount ) / 2.0 / tPerturbation;
                                    tJacProductionTermFD( tCount )         = tJacProductionTermFD( tCount ) / 2.0 / tPerturbation;
                                    tJacDivFluxFD( tCount )                = tJacDivFluxFD( tCount ) / 2.0 / tPerturbation;
                                    tJacSUPGFD( tCount )                   = tJacSUPGFD( tCount ) / 2.0 / tPerturbation;
                                    tJacStrongFormFD( tCount )             = tJacStrongFormFD( tCount ) / 2.0 / tPerturbation;

                                    tLeaderDOFHatWallDist( iCoef, iDim ) += tPerturbation;
                                    tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatWallDist );
                                    tLeaderFIs( 2 )->reset_eval_flags();
                                    tCMLeaderSATurbulence->reset_eval_flags();
                                    tSPSUPG->reset_eval_flags();

                                    tCount++;
                                }
                            }

                            real tRelError;
                            for ( uint iI = 0; iI < tCount; ++iI )
                            {
                                tRelError = std::abs( tJacDestructionCoefficient( iI ) - tJacDestructionCoefficientFD( iI ) ) / ( std::abs( tJacDestructionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDestructionCoefficient wall distance dof: ana = "
                                              << tJacDestructionCoefficient( iI ) << " FD = "
                                              << tJacDestructionCoefficientFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacProductionCoefficient( iI ) - tJacProductionCoefficientFD( iI ) ) / ( std::abs( tJacProductionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacProductionCoefficient wall distance dof: ana = "
                                              << tJacProductionCoefficient( iI ) << " FD = "
                                              << tJacProductionCoefficientFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacDiffusionCoefficient( iI ) - tJacDiffusionCoefficientFD( iI ) ) / ( std::abs( tJacDiffusionCoefficientFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDiffusionCoefficientFD wall distance dof:"
                                              << " ana = " << tJacDiffusionCoefficient( iI )
                                              << " FD  = " << tJacDiffusionCoefficientFD( iI )
                                              << " rel error = " << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacProductionTerm( iI ) - tJacProductionTermFD( iI ) ) / ( std::abs( tJacProductionTermFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacProductionTerm wall distance dof: ana = "
                                              << tJacProductionTerm( iI ) << " FD = "
                                              << tJacProductionTermFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacDestructionTerm( iI ) - tJacDestructionTermFD( iI ) ) / ( std::abs( tJacDestructionTermFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDestructionTerm wall distance dof: ana = "
                                              << tJacDestructionTerm( iI ) << " FD = "
                                              << tJacDestructionTermFD( iI ) << " rel error = "
                                              << tRelError
                                              << std::endl;
                                }

                                tRelError = std::abs( tJacDivFlux( iI ) - tJacDivFluxFD( iI ) ) / ( std::abs( tJacDivFluxFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacDivFluxFD wall distance dof:"
                                              << " ana = " << tJacDivFlux( iI )
                                              << " FD  = " << tJacDivFluxFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }

                                tRelError = std::abs( tJacSUPG( iI ) - tJacSUPGFD( iI ) ) / ( std::abs( tJacSUPGFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacSUPGFD wall distance dof:"
                                              << " ana = " << tJacSUPG( iI )
                                              << " FD  = " << tJacSUPGFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }

                                tRelError = std::abs( tJacStrongForm( iI ) - tJacStrongFormFD( iI ) ) / ( std::abs( tJacStrongFormFD( iI ) ) + MORIS_REAL_EPS );
                                if ( tRelError > tErrorLimit )
                                {
                                    std::cout << "tJacStrongForm wall distance dof:"
                                              << " ana = " << tJacStrongForm( iI )
                                              << " FD  = " << tJacStrongFormFD( iI )
                                              << " rel error = " << tRelError << std::endl;
                                }
                            }
                        }

                        //-------------------------------------------------------------------------------------------------------------------------

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
                                true,
                                false );

                        // print for debug
                        if ( !tCheckJacobian )
                        {
                            std::cout << "Case: Geometry " << iSpaceDim
                                      << " Order " << iInterpOrder
                                      << " VelConf " << iCu
                                      << " PresConf " << iCp
                                      << " iGP " << iGP << std::endl;
                        }

                        // require check is true
                        REQUIRE( tCheckJacobian );
                    }

                    // clean up
                    tLeaderFIs.clear();
                }
            }
        }
    }
} /*END_TEST_CASE*/

TEST_CASE( "IWG_Spalart_Allmaras_Turbulence_Bulk_Negative",
        "[IWG_Spalart_Allmaras_Turbulence_Bulk_Negative]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

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

    // dof type list
    moris::Cell< MSI::Dof_Type > tVelDofTypes = { MSI::Dof_Type::VX };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVisDofTypes = { { MSI::Dof_Type::VISCOSITY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes    = { tVelDofTypes, tVisDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropWallDistance = std::make_shared< fem::Property >();
    tPropWallDistance->set_parameters( { { { 1.0 } } } );
    tPropWallDistance->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_val_function( tConstValFunc );
    tPropViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMLeaderSATurbulence->set_dof_type_list( tDofTypes );
    tCMLeaderSATurbulence->set_property( tPropWallDistance, "WallDistance" );
    tCMLeaderSATurbulence->set_property( tPropViscosity, "KinViscosity" );
    tCMLeaderSATurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPSUPG =
            tSPFactory.create_SP( fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
    tSPSUPG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tSPSUPG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence" );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPSUPG->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::SPALART_ALLMARAS_TURBULENCE_BULK );
    tIWG->set_residual_dof_type( tVisDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMLeaderSATurbulence, "SpalartAllmarasTurbulence" );
    tIWG->set_stabilization_parameter( tSPSUPG, "SUPG" );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

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

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };

                // set velocity dof types
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

                // set viscosity property parameters
                tPropViscosity->set_parameters( { { { 2.0 } }, { { 0.0 }, { 0.0 } } } );
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
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };

                // set viscosity property parameters
                tPropViscosity->set_parameters( { { { 2.0 } }, { { 0.0 }, { 0.0 }, { 0.0 } } } );
                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        Matrix< DDRMat > tGravity( iSpaceDim, 1, 10.0 );

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

        // set space dimension to CM, SP
        tCMLeaderSATurbulence->set_space_dim( iSpaceDim );
        tSPSUPG->set_space_dim( iSpaceDim );

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
            int tNumDofVel = tNumCoeff * iSpaceDim;
            int tNumDofVis = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // loop over different state variable configurations
            for ( uint iCu = 0; iCu < 4; ++iCu )
            {
                for ( uint iCp = 0; iCp < 3; ++iCp )
                {
                    // fill coefficients for leader FI
                    Matrix< DDRMat > tLeaderDOFHatVel;
                    fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tLeaderDOFHatVis;
                    fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );
                    tLeaderDOFHatVis = -1.0 * tLeaderDOFHatVis;

                    // create a cell of field interpolators for IWG
                    Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

                    // create the field interpolator velocity
                    tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

                    // create the field interpolator viscosity
                    tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVis );

                    // set size and fill the set residual assembly map
                    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofVis - 1 } };

                    // set size and fill the set jacobian assembly map
                    Matrix< DDSMat > tJacAssembly = {
                        { 0, tNumDofVel - 1 },
                        { tNumDofVel, tNumDofVel + tNumDofVis - 1 }
                    };
                    tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
                    tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

                    // set size and init the set residual and jacobian
                    tIWG->mSet->mResidual.resize( 1 );
                    tIWG->mSet->mResidual( 0 ).set_size( tNumDofVel + tNumDofVis, 1, 0.0 );
                    tIWG->mSet->mJacobian.set_size( tNumDofVel + tNumDofVis, tNumDofVel + tNumDofVis, 0.0 );

                    // build global dof type list
                    tIWG->get_global_dof_type_list();

                    // populate the requested leader dof type
                    tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

                    // create a field interpolator manager
                    moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
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

                    // loop over integration points
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
                            std::cout << "Case: Geometry " << iSpaceDim
                                      << " Order " << iInterpOrder
                                      << " VelConf " << iCu
                                      << " PresConf " << iCp
                                      << " iGP " << iGP << std::endl;
                        }

                        // require check is true
                        REQUIRE( tCheckJacobian );
                    }

                    // clean up
                    tLeaderFIs.clear();
                }
            }
        }
    }
} /*END_TEST_CASE*/
