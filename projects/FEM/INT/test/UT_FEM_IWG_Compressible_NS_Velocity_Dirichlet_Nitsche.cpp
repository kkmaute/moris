/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.cpp
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
#include "cl_FEM_Cluster.hpp"
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
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Compressible_UT.cpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

using namespace moris;
using namespace fem;

void tMValFunc_Select_2D(
        moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =
          { { 1.0, 0.0 },
            { 0.0, 0.0 } };
}

void tMValFunc_Select_3D(
        moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =
          { { 1.0, 0.0, 0.0 },
            { 0.0, 0.0, 0.0 },
            { 0.0, 0.0, 1.0 } };
}

TEST_CASE( "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche_Ideal",
        "[IWG_Compressible_NS_Velocity_Dirichlet_Nitsche_Ideal]" )
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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDensityDof  = { { MSI::Dof_Type::RHO } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelocityDof = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tTempDof     = { { MSI::Dof_Type::TEMP } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes    = { tDensityDof( 0 ), tVelocityDof( 0 ), tTempDof( 0 ) };

    moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > tResDofTypes = { tDensityDof, tVelocityDof, tTempDof };

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

    // prescribed velocity
    std::shared_ptr< fem::Property > tPropVelocity = std::make_shared< fem::Property >();
    tPropVelocity->set_val_function( tConstValFunc );

    // define constitutive model and assign properties
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_COMPRESSIBLE_IDEAL );
    tCMLeaderFluid->set_dof_type_list( {tDensityDof( 0 ), tVelocityDof( 0 ), tTempDof( 0 ) } );
    tCMLeaderFluid->set_property( tPropViscosity,    "DynamicViscosity" );
    tCMLeaderFluid->set_property( tPropHeatCapacity, "IsochoricHeatCapacity" );
    tCMLeaderFluid->set_property( tPropGasConstant,  "SpecificGasConstant" );
    tCMLeaderFluid->set_property( tPropConductivity, "ThermalConductivity" );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::COMPRESSIBLE_VELOCITY_DIRICHLET_NITSCHE );
    tSPNitsche->set_property( tPropViscosity, "DynamicViscosity", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_parameters( { {{ 1.0 }}, {{ 1.0 }} } );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // loop over different Nitsche types
    for( uint iNitsche = 0; iNitsche < 2; iNitsche++ )
    {
        // output for debugging
//        std::cout << "-------------------------------------------------------------------\n" << std::flush;
//        std::cout << "Symmetric/Unsymmetric Nitsche (0:Symmetric, 1:Unsymmetric) : " << iNitsche << "\n" << std::flush;

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        // create the two different
        std::shared_ptr< fem::IWG > tIWG;
        if ( iNitsche == 0 )
            tIWG = tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE );
        else // iNitsche == 1
            tIWG = tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );

        tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
        tIWG->set_property( tPropVelocity, "PrescribedValue" );
        tIWG->set_constitutive_model( tCMLeaderFluid, "Fluid" );
        tIWG->set_stabilization_parameter( tSPNitsche, "NitschePenaltyParameter" );

        // loop over different residual dof Types
        for( uint iResDofType = 1; iResDofType < 3; iResDofType++ )
        {
            // output for debugging
//            std::cout << "-------------------------------------------------------------------\n" << std::flush;
//            std::cout << "Residual DoF type (0:density, 1:velocity, 2:Temperature) : " << iResDofType << "\n" << std::flush;

            // set residual dof type
            tIWG->set_residual_dof_type( tResDofTypes( iResDofType ) );
            tSPNitsche->set_dof_type_list( { tDofTypes( iResDofType ) }, mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            // set a fem set pointer

            MSI::Equation_Set * tSet = new fem::Set();
            static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::RHO ) )   = 0;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

            // set size and populate the set leader dof type map
            tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::RHO ) )   = 0;
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

            // loop on the space dimension
            for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
            {
                // output for debugging
//                std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                std::cout << "Performing Tests For Number of Spatial dimensions: " << iSpaceDim << "\n" << std::flush;

                // create normal for IWG
                Matrix< DDRMat > tNormal( iSpaceDim, 1, 3.8 );

                // prescribed velocity
                Matrix< DDRMat > tVelocity( iSpaceDim, 1, 2.7 );

                // switch on space dimension
                switch( iSpaceDim )
                {
                    case 2 :
                    {
                        // set geometry type
                        tGeometryType = mtk::Geometry_Type::QUAD;

                        // set velocity dof types
                        tVelocityDof = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };

                        // set normal
                        tNormal( 1 ) = -2.6;

                        // set prescribed velocity
                        tVelocity( 1 ) = 5.6;

                        break;
                    }
                    case 3 :
                    {
                        // set geometry type
                        tGeometryType = mtk::Geometry_Type::HEX;

                        // set velocity dof types
                        tVelocityDof = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };

                        // set normal
                        tNormal( 1 ) = -2.6;
                        tNormal( 2 ) =  9.4;

                        // set prescribed velocity
                        tVelocity( 1 ) = 5.6;
                        tVelocity( 2 ) = 7.4;

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

                // set space dimension to CM
                tPropVelocity->set_parameters( { tVelocity } );
                tCMLeaderFluid->set_space_dim( iSpaceDim );
                tSPNitsche->set_space_dim( iSpaceDim );

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
//                    std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                    std::cout << "Performing Tests For Interpolation Order:" << iInterpOrder << "\n\n" << std::flush;

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
                    int tNumDofRho  = tNumCoeff;
                    int tNumDofVel  = tNumCoeff * iSpaceDim;
                    int tNumDofTemp = tNumCoeff;
                    int tTotalNumDof = tNumDofRho + tNumDofVel + tNumDofTemp;

                    //create a space time interpolation rule
                    mtk::Interpolation_Rule tFIRule (
                            tGeometryType,
                            mtk::Interpolation_Type::LAGRANGE,
                            tInterpolationOrder,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    // fill coefficients for leader FI
                    Matrix< DDRMat > tLeaderDOFHatRho;
                    fill_RhoHat( tLeaderDOFHatRho, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tLeaderDOFHatVel;
                    fill_UHat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tLeaderDOFHatTemp;
                    fill_TempHat( tLeaderDOFHatTemp, iSpaceDim, iInterpOrder );

                    // create a cell of field interpolators for IWG
                    Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

                    // create the field interpolator density
                    tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tDensityDof( 0 ));
                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatRho );

                    // create the field interpolator velocity
                    tLeaderFIs( 1 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDof( 0 ) );
                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVel );

                    // create the field interpolator pressure
                    tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDof( 0 ) );
                    tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatTemp );

                    // set size and fill the set residual assembly map
                    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofRho - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofRho, tNumDofRho + tNumDofVel - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofRho + tNumDofVel, tTotalNumDof - 1 } };

                    // set size and fill the set jacobian assembly map
                    Matrix< DDSMat > tJacAssembly = {
                            { 0, tNumDofRho - 1 },
                            { tNumDofRho, tNumDofRho + tNumDofVel - 1 },
                            { tNumDofRho + tNumDofVel, tTotalNumDof - 1 } };
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
                    moris::Cell< moris::Cell< enum ge::PDV_Type > > tDummyDv;
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

                    // set the interpolator manager to the set
                    tCMLeaderFluid->mSet->mLeaderFIManager = &tFIManager;

                    // set IWG field interpolator manager
                    tCMLeaderFluid->set_field_interpolator_manager( &tFIManager );

                    // loop over integration points
                    uint tNumGPs = tIntegPoints.n_cols();
                    for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
                    {
                        // output for debugging
//                        std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                        std::cout << "Looping over Gauss points. Current GP-#: " << iGP << "\n\n" << std::flush;

                        // reset CM evaluation flags
                        tCMLeaderFluid->reset_eval_flags();

                        // reset IWG evaluation flags
                        tIWG->reset_eval_flags();

                        // create evaluation point xi, tau
                        Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                        // set integration point
                        tCMLeaderFluid->mSet->mLeaderFIManager->set_space_time( tParamPoint );
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
                            std::cout<<"Case: Residual DoF type " << iResDofType << " Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
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
}/*END_TEST_CASE*/

//------------------------------------------------------------------------------

TEST_CASE( "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche_VdW",
        "[IWG_Compressible_NS_Velocity_Dirichlet_Nitsche_VdW]" )
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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDensityDof  = { { MSI::Dof_Type::RHO } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelocityDof = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tTempDof     = { { MSI::Dof_Type::TEMP } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes    = { tDensityDof( 0 ), tVelocityDof( 0 ), tTempDof( 0 ) };

    moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > tResDofTypes = { tDensityDof, tVelocityDof, tTempDof };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties

    // isochoric heat capacity
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 5.7 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );

    // specific gas constant
    std::shared_ptr< fem::Property > tPropGasConstant = std::make_shared< fem::Property >();
    tPropGasConstant->set_parameters( { {{ 2.8 }} } );
    tPropGasConstant->set_val_function( tConstValFunc );

    // dynamic viscosity
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 11.9 }} } );
    tPropViscosity->set_val_function( tConstValFunc );

    // thermal conductivity
    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 3.7 }} } );
    tPropConductivity->set_val_function( tConstValFunc );

    // Capillarity Coefficient
    std::shared_ptr< fem::Property > tPropCapillarity = std::make_shared< fem::Property >();
    tPropCapillarity->set_parameters( { {{ 2.1 }} } );
    tPropCapillarity->set_val_function( tConstValFunc );

    // First VdW constant
    std::shared_ptr< fem::Property > tPropFirstVdWconst = std::make_shared< fem::Property >();
    tPropFirstVdWconst->set_parameters( { {{ 9.1 }} } );
    tPropFirstVdWconst->set_val_function( tConstValFunc );

    // Second VdW constant
    std::shared_ptr< fem::Property > tPropSecondVdWconst = std::make_shared< fem::Property >();
    tPropSecondVdWconst->set_parameters( { {{ 6.9 }} } );
    tPropSecondVdWconst->set_val_function( tConstValFunc );

    // prescribed velocity
    std::shared_ptr< fem::Property > tPropVelocity = std::make_shared< fem::Property >();
    tPropVelocity->set_val_function( tConstValFunc );

    // define constitutive model and assign properties
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_COMPRESSIBLE_VDW );
    tCMLeaderFluid->set_dof_type_list( {tDensityDof( 0 ), tVelocityDof( 0 ), tTempDof( 0 ) } );
    tCMLeaderFluid->set_property( tPropHeatCapacity,   "IsochoricHeatCapacity" );
    tCMLeaderFluid->set_property( tPropGasConstant,    "SpecificGasConstant" );
    tCMLeaderFluid->set_property( tPropViscosity,      "DynamicViscosity" );
    tCMLeaderFluid->set_property( tPropConductivity,   "ThermalConductivity" );
    tCMLeaderFluid->set_property( tPropCapillarity,    "CapillarityCoefficient" );
    tCMLeaderFluid->set_property( tPropFirstVdWconst,  "FirstVdWconstant" );
    tCMLeaderFluid->set_property( tPropSecondVdWconst, "SecondVdWconstant" );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::COMPRESSIBLE_VELOCITY_DIRICHLET_NITSCHE );
    tSPNitsche->set_property( tPropViscosity, "DynamicViscosity", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_parameters( { {{ 1.0 }}, {{ 1.0 }} } );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // loop over different Nitsche types
    for( uint iNitsche = 0; iNitsche < 2; iNitsche++ )
    {
        // output for debugging
//        std::cout << "-------------------------------------------------------------------\n" << std::flush;
//        std::cout << "Symmetric/Unsymmetric Nitsche (0:Symmetric, 1:Unsymmetric) : " << iNitsche << "\n" << std::flush;

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        // create the two different
        std::shared_ptr< fem::IWG > tIWG;
        if ( iNitsche == 0 )
            tIWG = tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE );
        else // iNitsche == 1
            tIWG = tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );

        tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
        tIWG->set_property( tPropVelocity, "PrescribedValue" );
        tIWG->set_constitutive_model( tCMLeaderFluid, "Fluid" );
        tIWG->set_stabilization_parameter( tSPNitsche, "NitschePenaltyParameter" );

        // loop over different residual dof Types
        for( uint iResDofType = 1; iResDofType < 3; iResDofType++ )
        {
            // output for debugging
//            std::cout << "-------------------------------------------------------------------\n" << std::flush;
//            std::cout << "Residual DoF type (0:density, 1:velocity, 2:Temperature) : " << iResDofType << "\n" << std::flush;

            // set residual dof type
            tIWG->set_residual_dof_type( tResDofTypes( iResDofType ) );
            tSPNitsche->set_dof_type_list( { tDofTypes( iResDofType ) }, mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            // set a fem set pointer

            MSI::Equation_Set * tSet = new fem::Set();
            static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::RHO ) )   = 0;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

            // set size and populate the set leader dof type map
            tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::RHO ) )   = 0;
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

            // loop on the space dimension
            for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
            {
                // output for debugging
//                std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                std::cout << "Performing Tests For Number of Spatial dimensions: " << iSpaceDim << "\n" << std::flush;

                // create normal for IWG
                Matrix< DDRMat > tNormal( iSpaceDim, 1, 3.8 );

                // prescribed velocity
                Matrix< DDRMat > tVelocity( iSpaceDim, 1, 2.7 );

                // switch on space dimension
                switch( iSpaceDim )
                {
                    case 2 :
                    {
                        // set geometry type
                        tGeometryType = mtk::Geometry_Type::QUAD;

                        // set velocity dof types
                        tVelocityDof = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY  } };

                        // set normal
                        tNormal( 1 ) = -2.6;

                        // set prescribed velocity
                        tVelocity( 1 ) = 5.6;

                        break;
                    }
                    case 3 :
                    {
                        // set geometry type
                        tGeometryType = mtk::Geometry_Type::HEX;

                        // set velocity dof types
                        tVelocityDof = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };

                        // set normal
                        tNormal( 1 ) = -2.6;
                        tNormal( 2 ) =  9.4;

                        // set prescribed velocity
                        tVelocity( 1 ) = 5.6;
                        tVelocity( 2 ) = 7.4;

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

                // set space dimension to CM
                tPropVelocity->set_parameters( { tVelocity } );
                tCMLeaderFluid->set_space_dim( iSpaceDim );
                tSPNitsche->set_space_dim( iSpaceDim );

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
//                    std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                    std::cout << "Performing Tests For Interpolation Order:" << iInterpOrder << "\n\n" << std::flush;

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
                    int tNumDofRho  = tNumCoeff;
                    int tNumDofVel  = tNumCoeff * iSpaceDim;
                    int tNumDofTemp = tNumCoeff;
                    int tTotalNumDof = tNumDofRho + tNumDofVel + tNumDofTemp;

                    //create a space time interpolation rule
                    mtk::Interpolation_Rule tFIRule (
                            tGeometryType,
                            mtk::Interpolation_Type::LAGRANGE,
                            tInterpolationOrder,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    // fill coefficients for leader FI
                    Matrix< DDRMat > tLeaderDOFHatRho;
                    fill_RhoHat( tLeaderDOFHatRho, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tLeaderDOFHatVel;
                    fill_UHat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tLeaderDOFHatTemp;
                    fill_TempHat( tLeaderDOFHatTemp, iSpaceDim, iInterpOrder );

                    // create a cell of field interpolators for IWG
                    Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

                    // create the field interpolator density
                    tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tDensityDof( 0 ) );
                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatRho );

                    // create the field interpolator velocity
                    tLeaderFIs( 1 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDof( 0 ) );
                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVel );

                    // create the field interpolator pressure
                    tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDof( 0 ) );
                    tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatTemp );

                    // set size and fill the set residual assembly map
                    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofRho - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofRho, tNumDofRho + tNumDofVel - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofRho + tNumDofVel, tTotalNumDof - 1 } };

                    // set size and fill the set jacobian assembly map
                    Matrix< DDSMat > tJacAssembly = {
                            { 0, tNumDofRho - 1 },
                            { tNumDofRho, tNumDofRho + tNumDofVel - 1 },
                            { tNumDofRho + tNumDofVel, tTotalNumDof - 1 } };
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
                    moris::Cell< moris::Cell< enum ge::PDV_Type > > tDummyDv;
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

                    // set the interpolator manager to the set
                    tCMLeaderFluid->mSet->mLeaderFIManager = &tFIManager;

                    // set IWG field interpolator manager
                    tCMLeaderFluid->set_field_interpolator_manager( &tFIManager );

                    // loop over integration points
                    uint tNumGPs = tIntegPoints.n_cols();
                    for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
                    {
                        // output for debugging
//                        std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                        std::cout << "Looping over Gauss points. Current GP-#: " << iGP << "\n\n" << std::flush;

                        // reset CM evaluation flags
                        tCMLeaderFluid->reset_eval_flags();

                        // reset IWG evaluation flags
                        tIWG->reset_eval_flags();

                        // create evaluation point xi, tau
                        Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                        // set integration point
                        tCMLeaderFluid->mSet->mLeaderFIManager->set_space_time( tParamPoint );
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
                            std::cout<<"Case: Residual DoF type " << iResDofType << " Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
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
}/*END_TEST_CASE*/

//------------------------------------------------------------------------------

TEST_CASE( "IWG_Compressible_NS_Velocity_Dirichlet_Nitsche_Ideal_Select",
        "[IWG_Compressible_NS_Velocity_Dirichlet_Nitsche_Ideal_Select]" )
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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDensityDof  = { { MSI::Dof_Type::RHO } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelocityDof = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tTempDof     = { { MSI::Dof_Type::TEMP } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes    = { tDensityDof( 0 ), tVelocityDof( 0 ), tTempDof( 0 ) };

    moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > tResDofTypes = { tDensityDof, tVelocityDof, tTempDof };

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

    // prescribed velocity
    std::shared_ptr< fem::Property > tPropVelocity = std::make_shared< fem::Property >();
    tPropVelocity->set_val_function( tConstValFunc );

    // select matrix
    std::shared_ptr< fem::Property > tPropSelect = std::make_shared< fem::Property >();
    tPropSelect->set_parameters( { {{ 1.0 }} } );
    tPropSelect->set_val_function( tMValFunc_Select_2D );

    // define constitutive model and assign properties
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_COMPRESSIBLE_IDEAL );
    tCMLeaderFluid->set_dof_type_list( {tDensityDof( 0 ), tVelocityDof( 0 ), tTempDof( 0 ) } );
    tCMLeaderFluid->set_property( tPropViscosity,    "DynamicViscosity" );
    tCMLeaderFluid->set_property( tPropHeatCapacity, "IsochoricHeatCapacity" );
    tCMLeaderFluid->set_property( tPropGasConstant,  "SpecificGasConstant" );
    tCMLeaderFluid->set_property( tPropConductivity, "ThermalConductivity" );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::COMPRESSIBLE_VELOCITY_DIRICHLET_NITSCHE );
    tSPNitsche->set_property( tPropViscosity, "DynamicViscosity", mtk::Leader_Follower::LEADER );
    tSPNitsche->set_parameters( { {{ 1.0 }}, {{ 1.0 }} } );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // loop over different Nitsche types
    for( uint iNitsche = 0; iNitsche < 2; iNitsche++ )
    {
        // output for debugging
//        std::cout << "-------------------------------------------------------------------\n" << std::flush;
//        std::cout << "Symmetric/Unsymmetric Nitsche (0:Symmetric, 1:Unsymmetric) : " << iNitsche << "\n" << std::flush;

        // define the IWGs
        fem::IWG_Factory tIWGFactory;

        // create the two different
        std::shared_ptr< fem::IWG > tIWG;
        if ( iNitsche == 0 )
            tIWG = tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE );
        else // iNitsche == 1
            tIWG = tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE );

        tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
        tIWG->set_property( tPropVelocity, "PrescribedValue" );
        tIWG->set_constitutive_model( tCMLeaderFluid, "Fluid" );
        tIWG->set_stabilization_parameter( tSPNitsche, "NitschePenaltyParameter" );

        // loop over different residual dof Types
        for( uint iResDofType = 1; iResDofType < 3; iResDofType++ )
        {
            // output for debugging
//            std::cout << "-------------------------------------------------------------------\n" << std::flush;
//            std::cout << "Residual DoF type (0:density, 1:velocity, 2:Temperature) : " << iResDofType << "\n" << std::flush;

            // set residual dof type
            tIWG->set_residual_dof_type( tResDofTypes( iResDofType ) );
            tSPNitsche->set_dof_type_list( { tDofTypes( iResDofType ) }, mtk::Leader_Follower::LEADER );

            //------------------------------------------------------------------------------
            // set a fem set pointer

            MSI::Equation_Set * tSet = new fem::Set();
            static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::SIDESET );
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::RHO ) )   = 0;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

            // set size and populate the set leader dof type map
            tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::RHO ) )   = 0;
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

            // loop on the space dimension
            for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
            {
                // output for debugging
//                std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                std::cout << "Performing Tests For Number of Spatial dimensions: " << iSpaceDim << "\n" << std::flush;

                // create normal for IWG
                Matrix< DDRMat > tNormal( iSpaceDim, 1, 3.8 );

                // prescribed velocity
                Matrix< DDRMat > tVelocity( iSpaceDim, 1, 2.7 );

                // switch on space dimension
                switch( iSpaceDim )
                {
                    case 2 :
                    {
                        // set geometry type
                        tGeometryType = mtk::Geometry_Type::QUAD;

                        // set velocity dof types
                        tVelocityDof = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };

                        // set normal
                        tNormal( 1 ) = -2.6;

                        // set prescribed velocity
                        tVelocity( 1 ) = 5.6;

                        // use correct select matrix
                        tPropSelect->set_val_function( tMValFunc_Select_2D );
                        tIWG->set_property( tPropSelect, "Select" );

                        break;
                    }
                    case 3 :
                    {
                        // set geometry type
                        tGeometryType = mtk::Geometry_Type::HEX;

                        // set velocity dof types
                        tVelocityDof = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };

                        // set normal
                        tNormal( 1 ) = -2.6;
                        tNormal( 2 ) =  9.4;

                        // set prescribed velocity
                        tVelocity( 1 ) = 5.6;
                        tVelocity( 2 ) = 7.4;

                        // use correct select matrix
                        tPropSelect->set_val_function( tMValFunc_Select_3D );
                        tIWG->set_property( tPropSelect, "Select" );

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

                // set space dimension to CM
                tPropVelocity->set_parameters( { tVelocity } );
                tCMLeaderFluid->set_space_dim( iSpaceDim );
                tSPNitsche->set_space_dim( iSpaceDim );

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
//                    std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                    std::cout << "Performing Tests For Interpolation Order:" << iInterpOrder << "\n\n" << std::flush;

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
                    int tNumDofRho  = tNumCoeff;
                    int tNumDofVel  = tNumCoeff * iSpaceDim;
                    int tNumDofTemp = tNumCoeff;
                    int tTotalNumDof = tNumDofRho + tNumDofVel + tNumDofTemp;

                    //create a space time interpolation rule
                    mtk::Interpolation_Rule tFIRule (
                            tGeometryType,
                            mtk::Interpolation_Type::LAGRANGE,
                            tInterpolationOrder,
                            mtk::Interpolation_Type::LAGRANGE,
                            mtk::Interpolation_Order::LINEAR );

                    // fill coefficients for leader FI
                    Matrix< DDRMat > tLeaderDOFHatRho;
                    fill_RhoHat( tLeaderDOFHatRho, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tLeaderDOFHatVel;
                    fill_UHat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tLeaderDOFHatTemp;
                    fill_TempHat( tLeaderDOFHatTemp, iSpaceDim, iInterpOrder );

                    // create a cell of field interpolators for IWG
                    Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

                    // create the field interpolator density
                    tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tDensityDof( 0 ) );
                    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatRho );

                    // create the field interpolator velocity
                    tLeaderFIs( 1 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDof( 0 ) );
                    tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatVel );

                    // create the field interpolator pressure
                    tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDof( 0 ) );
                    tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatTemp );

                    // set size and fill the set residual assembly map
                    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
                    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofRho - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofRho, tNumDofRho + tNumDofVel - 1 } };
                    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofRho + tNumDofVel, tTotalNumDof - 1 } };

                    // set size and fill the set jacobian assembly map
                    Matrix< DDSMat > tJacAssembly = {
                            { 0, tNumDofRho - 1 },
                            { tNumDofRho, tNumDofRho + tNumDofVel - 1 },
                            { tNumDofRho + tNumDofVel, tTotalNumDof - 1 } };
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
                    moris::Cell< moris::Cell< enum ge::PDV_Type > > tDummyDv;
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

                    // set the interpolator manager to the set
                    tCMLeaderFluid->mSet->mLeaderFIManager = &tFIManager;

                    // set IWG field interpolator manager
                    tCMLeaderFluid->set_field_interpolator_manager( &tFIManager );

                    // loop over integration points
                    uint tNumGPs = tIntegPoints.n_cols();
                    for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
                    {
                        // output for debugging
//                        std::cout << "-------------------------------------------------------------------\n" << std::flush;
//                        std::cout << "Looping over Gauss points. Current GP-#: " << iGP << "\n\n" << std::flush;

                        // reset CM evaluation flags
                        tCMLeaderFluid->reset_eval_flags();

                        // reset IWG evaluation flags
                        tIWG->reset_eval_flags();

                        // create evaluation point xi, tau
                        Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                        // set integration point
                        tCMLeaderFluid->mSet->mLeaderFIManager->set_space_time( tParamPoint );
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
                            std::cout<<"Case: Residual DoF type " << iResDofType << " Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
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
}/*END_TEST_CASE*/

