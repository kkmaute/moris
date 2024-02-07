/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Compressible_NS_Dirichlet_Nitsche.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private
//LINALG/src
#include "op_equal_equal.hpp"
//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_MM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Compressible_UT.cpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Compressible_NS_Dirichlet_Nitsche_Pressure_Primitive",
        "[IWG_Compressible_NS_Dirichlet_Nitsche_Pressure_Primitive]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-6;
    real tEpsilonCubic = 3.0E-6;

    // define a perturbation relative size
    real tPerturbation = 2.0E-4;
    real tPerturbationCubic = 5.0E-4;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

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
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    Vector< MSI::Dof_Type > tPressureDof = { MSI::Dof_Type::P };
    Vector< MSI::Dof_Type > tVelocityDof = { MSI::Dof_Type::VX };
    Vector< MSI::Dof_Type > tTempDof     = { MSI::Dof_Type::TEMP };

    Vector< Vector< MSI::Dof_Type > > tDofTypes         = { tPressureDof, tVelocityDof, tTempDof };
    Vector< Vector< MSI::Dof_Type > > tResidualDofTypes = tDofTypes;

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties

    // dynamic viscosity
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 11.9 }} } );
    tPropViscosity->set_val_function( tConstValFunc );

    // isochoric heat capacity
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 5.7 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );

    // specific gas constant
    std::shared_ptr< fem::Property > tPropGasConstant = std::make_shared< fem::Property >();
    tPropGasConstant->set_parameters( { {{ 2.8 }} } );
    tPropGasConstant->set_val_function( tConstValFunc );

    // thermal conductivity
    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 3.7 }} } );
    tPropConductivity->set_val_function( tConstValFunc );

    //------------------------------------------------------------------------------
    // prescribed values

    // prescribe pressure
    std::shared_ptr< fem::Property > tPropPrescPres = std::make_shared< fem::Property >();
    tPropPrescPres->set_parameters( { {{ 5.1 }} } );
    tPropPrescPres->set_val_function( tConstValFunc );

    // prescribed velocity
    std::shared_ptr< fem::Property > tPropPrescVel = std::make_shared< fem::Property >();
    tPropPrescVel->set_val_function( tConstValFunc );

    // selection matrix for velocity
    std::shared_ptr< fem::Property > tPropVelSelect = std::make_shared< fem::Property >();
    tPropVelSelect->set_val_function( tConstValFunc );

    // prescribed temperature
    std::shared_ptr< fem::Property > tPropPrescTemp = std::make_shared< fem::Property >();
    tPropPrescTemp->set_parameters( { {{ 3.8 }} } );
    tPropPrescTemp->set_val_function( tConstValFunc );

    // upwinding weight
    std::shared_ptr< fem::Property > tPropUpwind = std::make_shared< fem::Property >();
    tPropUpwind->set_parameters( { {{ 1.3 }} } );
    tPropUpwind->set_val_function( tConstValFunc );

    // define material model and assign properties
    fem::MM_Factory tMMFactory;

    std::shared_ptr< fem::Material_Model > tMMFluid =
            tMMFactory.create_MM( fem::Material_Type::PERFECT_GAS );
    tMMFluid->set_dof_type_list( {tPressureDof, tTempDof } );
    tMMFluid->set_property( tPropHeatCapacity, "IsochoricHeatCapacity" );
    tMMFluid->set_property( tPropGasConstant,  "SpecificGasConstant" );

    // define constitutive model and assign properties
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_COMPRESSIBLE_NEWTONIAN );
    tCMLeaderFluid->set_dof_type_list( {tPressureDof, tVelocityDof, tTempDof } );
    tCMLeaderFluid->set_property( tPropViscosity,    "DynamicViscosity" );
    tCMLeaderFluid->set_property( tPropConductivity, "ThermalConductivity" );
    tCMLeaderFluid->set_material_model( tMMFluid, "ThermodynamicMaterialModel" );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::COMPRESSIBLE_DIRICHLET_NITSCHE );
    tSPNitsche->set_property( tPropViscosity, "DynamicViscosity", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_property( tPropConductivity, "ThermalConductivity", mtk::Leader_Follower::LEADER );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    // loop over two different configurations for the IWG
    for( uint iIWG = 0; iIWG < 2; iIWG++ )
    {
        // output for debugging
        // std::cout << "----------------------------------------------------------------------------------\n" << std::flush;
        // std::cout << "Performing Tests for IWG ( 0 : All fields prescribed, 1: Some fields prescribed + Upwinding ) : " << iIWG << "\n" << std::flush;

        // create IWG (pointer)
        std::shared_ptr< fem::IWG > tIWG =
                tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tIWG->set_residual_dof_type( tResidualDofTypes );
        tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
        tIWG->set_property( tPropViscosity,    "DynamicViscosity" );
        tIWG->set_property( tPropConductivity, "ThermalConductivity" );
        tIWG->set_material_model( tMMFluid, "FluidMM" );
        tIWG->set_constitutive_model( tCMLeaderFluid, "FluidCM" );
        tIWG->set_stabilization_parameter( tSPNitsche, "NitschePenaltyParameter" );

        // test different configurations
        switch ( iIWG )
        {
            case 0 :
            {
                tIWG->set_property( tPropPrescPres, "PrescribedDof1" );
                tIWG->set_property( tPropPrescVel,  "PrescribedVelocity" );
                tIWG->set_property( tPropPrescTemp, "PrescribedDof3" );
                break;
            }
            case 1 :
            {
                tIWG->set_property( tPropPrescPres, "PrescribedDof1" );
                tIWG->set_property( tPropPrescVel,  "PrescribedVelocity" );
                tIWG->set_property( tPropVelSelect, "SelectVelocity" );
                //tIWG->set_property( tPropPrescTemp, "PrescribedDof3" );
                tIWG->set_property( tPropUpwind,    "PressureUpwind" );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "UT - IWG_Compressible_NS_Boundary_Pressure_Primitive - iIWG out of bounds, loop incorrect" );
            }
        }

        //------------------------------------------------------------------------------
        // set a fem set pointer

        MSI::Equation_Set * tSet = new fem::Set();
        static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
        tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

        // set size for the set EqnObjDofTypeList
        tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

        // set size and populate the set dof type map
        tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )     = 0;
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
        tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

        // set size and populate the set leader dof type map
        tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )     = 0;
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
        tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

        // loop on the space dimension
        for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
        {
            // output for debugging
            // std::cout << "-------------------------------------------------------------------\n" << std::flush;
            // std::cout << "Performing Tests For Number of Spatial dimensions: " << iSpaceDim << "\n" << std::flush;
            // std::cout << "-------------------------------------------------------------------\n\n" << std::flush;

            // prescribed velocity
            Matrix< DDRMat > tPrescVel( iSpaceDim, 1, 2.7 );

            // select matrix for velocity
            Matrix< DDRMat > tSelectMat;
            eye( iSpaceDim, iSpaceDim, tSelectMat );
            tSelectMat( 1, 1 ) = 0.0;

            // weights for nitsche SP
            Matrix< DDRMat > tOnes( iSpaceDim + 2, 1, 1.0 );
            Vector< Matrix< DDRMat > > tSPWeights( 1, tOnes );

            // create normal for IWG
            Matrix< DDRMat > tNormal( iSpaceDim, 1, 3.8 );

            // switch on space dimension
            switch( iSpaceDim )
            {
                case 2 :
                {
                    // set geometry type
                    tGeometryType = mtk::Geometry_Type::QUAD;

                    // set velocity dof types
                    tVelocityDof = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

                    // set prescribed velocity
                    tPrescVel( 1 ) = 5.6;

                    // set normal
                    tNormal( 1 ) = -2.6;

                    break;
                }
                case 3 :
                {
                    // set geometry type
                    tGeometryType = mtk::Geometry_Type::HEX;

                    // set velocity dof types
                    tVelocityDof = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };

                    // set prescribed velocity
                    tPrescVel( 1 ) = 5.6;
                    tPrescVel( 2 ) = 7.4;

                    // set normal
                    tNormal( 1 ) = -2.6;
                    tNormal( 2 ) =  9.4;

                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "Unit Test Error: QUAD or HEX only." );
                    break;
                }
            }

            // assign normal to IWG
            tNormal = tNormal / norm( tNormal );
            tIWG->set_normal( tNormal );

            // set space dimension to CM and SP
            tMMFluid->set_space_dim( iSpaceDim );
            tCMLeaderFluid->set_space_dim( iSpaceDim );
            tSPNitsche->set_space_dim( iSpaceDim );

            // set prescribed Values
            tPropPrescVel->set_parameters( { tPrescVel } );
            tPropVelSelect->set_parameters( { tSelectMat } );
            tSPNitsche->set_parameters( tSPWeights );

            // loop on the interpolation order
            for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
            {
                // tune finite differencing for cubic shape functions
                if ( iInterpOrder == 3 )
                {
                    tEpsilon = tEpsilonCubic;
                    tPerturbation = tPerturbationCubic;
                }

                // output for debugging
                // std::cout << "-------------------------------------------------------------------\n" << std::flush;
                // std::cout << "-------------------------------------------------------------------\n" << std::flush;
                // std::cout << "Performing Tests For Interpolation Order:" << iInterpOrder << "\n\n" << std::flush;

                // create an interpolation order
                mtk::Interpolation_Order tGIInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

                //------------------------------------------------------------------------------
                // space and time geometry interpolators
                // create a space geometry interpolation rule
                mtk::Interpolation_Rule tGIRule( tGeometryType,
                        mtk::Interpolation_Type::LAGRANGE,
                        tGIInterpolationOrder,
                        mtk::Interpolation_Type::LAGRANGE,
                        mtk::Interpolation_Order::LINEAR );

                // create a space time geometry interpolator
                Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

                // create time coeff tHat
                Matrix< DDRMat > tTHat = {{ 0.26 }, { 0.87 }};

                Matrix< DDRMat > tXHat;
                fill_xhat( tXHat, iSpaceDim, iInterpOrder );

                // set the coefficients xHat, tHat
                tGI.set_coeff( tXHat, tTHat );

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

                //------------------------------------------------------------------------------
                // field interpolators
                // create an interpolation order
                mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

                // number of dof for interpolation order
                uint tNumCoeff = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

                // get number of dof per type
                int tNumDofPres = tNumCoeff;
                int tNumDofVel  = tNumCoeff * iSpaceDim;
                int tNumDofTemp = tNumCoeff;
                int tTotalNumDof = tNumDofPres + tNumDofVel + tNumDofTemp;

                //create a space time interpolation rule
                mtk::Interpolation_Rule tFIRule (
                        tGeometryType,
                        mtk::Interpolation_Type::LAGRANGE,
                        tInterpolationOrder,
                        mtk::Interpolation_Type::LAGRANGE,
                        mtk::Interpolation_Order::LINEAR );

                // fill coefficients for leader FI
                Matrix< DDRMat > tLeaderDOFHatP;
                fill_RhoHat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );
                Matrix< DDRMat > tLeaderDOFHatVel;
                fill_UHat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
                Matrix< DDRMat > tLeaderDOFHatTemp;
                fill_TempHat( tLeaderDOFHatTemp, iSpaceDim, iInterpOrder );

                // create a cell of field interpolators for IWG
                Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

                // create the field interpolator density
                tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPressureDof );
                tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatP );

                // create the field interpolator velocity
                tLeaderFIs( 1 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDof );
                tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVel );

                // create the field interpolator pressure
                tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDof );
                tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatTemp );

                // set size and fill the set residual assembly map
                tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
                tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofPres - 1 } };
                tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofPres, tNumDofPres + tNumDofVel - 1 } };
                tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofPres + tNumDofVel, tTotalNumDof - 1 } };

                // set size and fill the set jacobian assembly map
                Matrix< DDSMat > tJacAssembly = {
                        { 0, tNumDofPres - 1 },
                        { tNumDofPres, tNumDofPres + tNumDofVel - 1 },
                        { tNumDofPres + tNumDofVel, tTotalNumDof - 1 } };
                tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
                tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
                tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
                tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;

                // set size and init the set residual and jacobian
                tIWG->mSet->mResidual.resize( 1 );
                tIWG->mSet->mResidual( 0 ).set_size(
                        tTotalNumDof,
                        1,
                        0.0 );
                tIWG->mSet->mJacobian.set_size(
                        tTotalNumDof,
                        tTotalNumDof,
                        0.0 );

                // build global dof type list
                tIWG->get_global_dof_type_list();

                // populate the requested leader dof type
                tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

                // create a field interpolator manager
                Vector< Vector< enum gen::PDV_Type > > tDummyDv;
                Vector< Vector< enum mtk::Field_Type > > tDummyField;
                Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

                // populate the field interpolator manager
                tFIManager.mFI = tLeaderFIs;
                tFIManager.mIPGeometryInterpolator = &tGI;
                tFIManager.mIGGeometryInterpolator = &tGI;

                // set the interpolator manager to the set
                tIWG->mSet->mLeaderFIManager = &tFIManager;

                // set IWG field interpolator manager
                tIWG->set_field_interpolator_manager( &tFIManager );

                // loop over integration points
                uint tNumGPs = tIntegPoints.n_cols();
                for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
                {
                    // output for debugging
                    // std::cout << "-------------------------------------------------------------------\n" << std::flush;
                    // std::cout << "Looping over Gauss points. Current GP-#: " << iGP << "\n\n" << std::flush;

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
                    tIWG->compute_jacobian( 1.0 );

                    // check evaluation of the jacobian by FD
                    //------------------------------------------------------------------------------
                    // reset jacobian
                    tIWG->mSet->mJacobian.fill( 0.0 );

                    // init the jacobian for IWG and FD evaluation
                    Matrix< DDRMat > tJacobian = tIWG->mSet->get_jacobian();
                    Matrix< DDRMat > tJacobianTest;
                    Matrix< DDRMat > tJacobianFD;

                    // check jacobian by FD
                    bool tCheckJacobian = tIWG->check_jacobian_multi_residual(
                            tPerturbation,
                            tEpsilon,
                            1.0,
                            tJacobianTest,
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

            } // loop over interpolation orders

        } // loop over spatial dimensions

    } // loop over different IWG configurations

}/*END_TEST_CASE*/

