/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Incompressible_NS_Pressure_Bulk.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

//MTK/src
#include "cl_MTK_Enums.hpp"
//LINALG/src
#include "op_equal_equal.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
//FEM//INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_MTK_Integrator.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Incompressible_NS_Pressure_Bulk_No_Turbulence",
        "[IWG_Incompressible_NS_Pressure_Bulk_No_Turbulence]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

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
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tTEMPDofTypes = { MSI::Dof_Type::TEMP };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tTEMPDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunc );
    //tPropViscosity->set_dof_type_list( { tPDofTypes } );
    //tPropViscosity->set_val_function( tPFIValFunc );
    //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 2.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //tPropDensity->set_dof_type_list( { tPDofTypes } );
    //tPropDensity->set_val_function( tPFIValFunc );
    //tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropGravity = std::make_shared< fem::Property >();
    tPropGravity->set_val_function( tConstValFunc );
    //tPropGravity->set_dof_type_list( { tPDofTypes } );
    //tPropGravity->set_val_function( tPFIValFunc );
    //tPropGravity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropThermalExp = std::make_shared< fem::Property >();
    tPropThermalExp->set_parameters( { {{ 23.0 }} } );
    tPropThermalExp->set_val_function( tConstValFunc );
    //tPropThermalExp->set_dof_type_list( { tTEMPDofTypes } );
    //tPropThermalExp->set_val_function( tTEMPFIValFunc );
    //tPropThermalExp->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropRefTemp = std::make_shared< fem::Property >();
    tPropRefTemp->set_parameters( { {{ 15.0 }} } );
    tPropRefTemp->set_val_function( tConstValFunc );
    //tPropRefTemp->set_dof_type_list( { tTEMPDofTypes } );
    //tPropRefTemp->set_val_function( tTEMPFIValFunc );
    //tPropRefTemp->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderIncFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMLeaderIncFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) } );
    tCMLeaderIncFluid->set_property( tPropViscosity, "Viscosity" );
    tCMLeaderIncFluid->set_property( tPropDensity, "Density" );
    tCMLeaderIncFluid->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow =
            tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
    tSPIncFlow->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) }, mtk::Leader_Follower::LEADER );
    tSPIncFlow->set_property( tPropDensity, "Density", mtk::Leader_Follower::LEADER );
    tSPIncFlow->set_property( tPropViscosity, "Viscosity", mtk::Leader_Follower::LEADER );
    tSPIncFlow->set_parameters( { {{ 36.0 }} } );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPIncFlow->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
    tIWG->set_residual_dof_type( tPDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropGravity, "Gravity" );
    tIWG->set_property( tPropThermalExp, "ThermalExpansion" );
    tIWG->set_property( tPropRefTemp, "ReferenceTemp" );
    tIWG->set_constitutive_model( tCMLeaderIncFluid, "IncompressibleFluid" );
    tIWG->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 2;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set geometry inputs
        //------------------------------------------------------------------------------
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

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };
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

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };
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
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tPropGravity->set_parameters( { tGravity } );
        tCMLeaderIncFluid->set_space_dim( iSpaceDim );
        tSPIncFlow->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder   = tIntegrationOrders( iSpaceDim - 2 );

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
            int tNumDofVel  = tNumCoeff * iSpaceDim;
            int tNumDofP    = tNumCoeff;
            int tNumDofTEMP = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_phat( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator temperature
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTEMPDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatTEMP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofP - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofP + tNumDofTEMP,
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    tNumDofVel + tNumDofP + tNumDofTEMP,
                    tNumDofVel + tNumDofP + tNumDofTEMP,
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop over integartion points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
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
                if( !tCheckJacobian )
                {
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "IWG_Incompressible_NS_Pressure_Bulk_With_Turbulence",
        "[IWG_Incompressible_NS_Pressure_Bulk_With_Turbulence]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

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
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tVisDofTypes  = { MSI::Dof_Type::VISCOSITY };
    moris::Cell< MSI::Dof_Type > tTEMPDofTypes = { MSI::Dof_Type::TEMP };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tTEMPDofTypes, tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_val_function( tConstValFunc );
    tPropViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );
    //tPropViscosity->set_dof_type_list( { tPDofTypes } );
    //tPropViscosity->set_val_function( tPFIValFunc );
    //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropKinViscosity = std::make_shared< fem::Property >();
    tPropKinViscosity->set_val_function( tConstValFunc );
    tPropKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 2.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //tPropDensity->set_dof_type_list( { tPDofTypes } );
    //tPropDensity->set_val_function( tPFIValFunc );
    //tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropGravity = std::make_shared< fem::Property >();
    tPropGravity->set_val_function( tConstValFunc );
    //tPropGravity->set_dof_type_list( { tPDofTypes } );
    //tPropGravity->set_val_function( tPFIValFunc );
    //tPropGravity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropThermalExp = std::make_shared< fem::Property >();
    tPropThermalExp->set_parameters( { {{ 23.0 }} } );
    tPropThermalExp->set_val_function( tConstValFunc );
    //tPropThermalExp->set_dof_type_list( { tTEMPDofTypes } );
    //tPropThermalExp->set_val_function( tTEMPFIValFunc );
    //tPropThermalExp->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropRefTemp = std::make_shared< fem::Property >();
    tPropRefTemp->set_parameters( { {{ 15.0 }} } );
    tPropRefTemp->set_val_function( tConstValFunc );
    //tPropRefTemp->set_dof_type_list( { tTEMPDofTypes } );
    //tPropRefTemp->set_val_function( tTEMPFIValFunc );
    //tPropRefTemp->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropMassSource = std::make_shared< fem::Property >();
    tPropMassSource->set_parameters( { {{ 3.5 }} } );
    tPropMassSource->set_val_function( tConstValFunc );
    //tPropMassSource->set_dof_type_list( { tTEMPDofTypes } );
    //tPropMassSource->set_val_function( tTEMPFIValFunc );
    //tPropMassSource->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropBodyLoad = std::make_shared< fem::Property >();
    tPropBodyLoad->set_val_function( tConstValFunc );
    //tPropMassSource->set_dof_type_list( { tTEMPDofTypes } );
    //tPropMassSource->set_val_function( tTEMPFIValFunc );
    //tPropMassSource->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMLeaderTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMLeaderTurbulence->set_property( tPropViscosity, "Viscosity" );
    tCMLeaderTurbulence->set_property( tPropKinViscosity, "KinViscosity" );
    tCMLeaderTurbulence->set_property( tPropDensity, "Density" );
    tCMLeaderTurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow =
            tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
    tSPIncFlow->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) }, mtk::Leader_Follower::LEADER );
    tSPIncFlow->set_property( tPropDensity, "Density", mtk::Leader_Follower::LEADER );
    tSPIncFlow->set_property( tPropViscosity, "Viscosity", mtk::Leader_Follower::LEADER );
    tSPIncFlow->set_parameters( { {{ 36.0 }} } );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPIncFlow->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
    tIWG->set_residual_dof_type( tPDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropGravity,        "Gravity" );
    tIWG->set_property( tPropThermalExp,     "ThermalExpansion" );
    tIWG->set_property( tPropRefTemp,        "ReferenceTemp" );
    tIWG->set_property( tPropMassSource,     "MassSource" );
    tIWG->set_property( tPropBodyLoad,       "Load" );

    tIWG->set_constitutive_model( tCMLeaderTurbulence, "IncompressibleFluid" );
    tIWG->set_stabilization_parameter( tSPIncFlow,     "IncompressibleFlow" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 2;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 3;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 2;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 3;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch( iSpaceDim )
        {
            case 2 :
            {
                // set momentum source term
                tPropBodyLoad->set_parameters( { { { 0.5 },{ 0.35 } } } );

                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                        { 1.0, 0.0 },
                        { 1.0, 1.0 },
                        { 0.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }};

                // set viscosity property parameters
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0}} } );

                break;
            }
            case 3 :
            {
                // set momentum source term
                tPropBodyLoad->set_parameters( { { { 0.5 },{ 0.35 },{ 0.42 } } } );

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

                // set velocity dof types
                tVelDofTypes = {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }};

                // set viscosity property parameters
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );
                tPropViscosity->set_parameters( { {{ 1.0 }}, {{0.0},{0.0},{0.0}} } );

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
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimension to CM, SP
        tPropGravity->set_parameters( { tGravity } );
        tCMLeaderTurbulence->set_space_dim( iSpaceDim );
        tSPIncFlow->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder   = tIntegrationOrders( iSpaceDim - 2 );

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
            int tNumDofVel  = tNumCoeff * iSpaceDim;
            int tNumDofP    = tNumCoeff;
            int tNumDofTEMP = tNumCoeff;
            int tNumDofVis  = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_phat( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatVis;
            fill_phat( tLeaderDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // create the field interpolator temperature
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTEMPDofTypes );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatTEMP );

            // create the field interpolator viscosity
            tLeaderFIs( 3 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tLeaderFIs( 3 )->set_coeff( tLeaderDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel-1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofP - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel - 1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofVel + tNumDofP, tNumDofVel + tNumDofP + tNumDofTEMP - 1 },
                    { tNumDofVel + tNumDofP + tNumDofTEMP, tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis,
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis,
                    tNumDofVel + tNumDofP + tNumDofTEMP + tNumDofVis,
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum gen::PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop over integartion points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
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
                if( !tCheckJacobian )
                {
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/

