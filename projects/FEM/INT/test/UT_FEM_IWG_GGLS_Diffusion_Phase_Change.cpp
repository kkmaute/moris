/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_GGLS_Diffusion_Phase_Change.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private
//MTK/src
#include "cl_MTK_Enums.hpp"
//LINALG/src
#include "op_equal_equal.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"

void tConstValFunction_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tGeoValFunction_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
}

void tFIValFunction_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) + 0.1 * aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tFIDerFunction_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = 0.1 * aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

void tFIDer0Function_UTIWGGGLSDIFFBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = 0.0 * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

namespace moris
{
    namespace fem
    {

        moris::Vector<bool> test_IWG_Diffusion_Phase_Change_GGLS(
                Matrix< DDRMat > aXHat,
                Matrix< DDRMat > aTHat,
                mtk::Interpolation_Rule aGIRule,
                mtk::Interpolation_Rule aFIRule,
                Matrix< DDRMat > aDOFHat,
                Matrix< DDRMat > aParamPoint,
                uint aNumDOFs,
                uint aSpatialDim = 2 )
        {
            // initialize cell of checks
            moris::Vector<bool> tChecks( 1, false );

            // define an epsilon environment
            real tEpsilonRel = 1.0E-6;

            // define a perturbation relative size
            real tPerturbation = 1.0E-6;

            // create the properties ------------------------------------------------------------------- //

            std::shared_ptr< fem::Property > tPropLeaderConductivity = std::make_shared< fem::Property > ();
            tPropLeaderConductivity->set_parameters( {{{ 1.1 }}, {{ 1.1 }}} );
            //    tPropLeaderConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            //    tPropLeaderConductivity->set_val_function( tFIValFunction_UTIWGDIFFBULK );
            //    tPropLeaderConductivity->set_dof_derivative_functions( { tFIDerFunction_UTIWGGGLSDIFFBULK } );
            tPropLeaderConductivity->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            std::shared_ptr< fem::Property > tPropLeaderDensity = std::make_shared< fem::Property > ();
            tPropLeaderDensity->set_parameters( { {{ 1.2 }} } );
            tPropLeaderDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropLeaderDensity->set_val_function( tFIValFunction_UTIWGGGLSDIFFBULK );
            tPropLeaderDensity->set_dof_derivative_functions( { tFIDerFunction_UTIWGGGLSDIFFBULK } );
//            tPropLeaderDensity->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );
//            tPropLeaderDensity->set_dof_derivative_functions( { tFIDer0Function_UTIWGGGLSDIFFBULK } );

            std::shared_ptr< fem::Property > tPropLeaderHeatCapacity = std::make_shared< fem::Property > ();
            tPropLeaderHeatCapacity->set_parameters( { {{ 0.3 }} } );
            tPropLeaderHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropLeaderHeatCapacity->set_val_function( tFIValFunction_UTIWGGGLSDIFFBULK );
            tPropLeaderHeatCapacity->set_dof_derivative_functions( { tFIDerFunction_UTIWGGGLSDIFFBULK } );
//            tPropLeaderHeatCapacity->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );
//            tPropLeaderHeatCapacity->set_dof_derivative_functions( { tFIDer0Function_UTIWGGGLSDIFFBULK } );

            // latent heat
            std::shared_ptr< fem::Property > tPropLeaderLatentHeat = std::make_shared< fem::Property >();
            tPropLeaderLatentHeat->set_parameters( {{{ 100.0}}} );
            //tPropLeaderLatentHeat->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropLeaderLatentHeat->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            // phase change temp
            std::shared_ptr< fem::Property > tPropLeaderTmelt = std::make_shared< fem::Property >();
            tPropLeaderTmelt->set_parameters( {{{ 5.2 }}} );
            //tPropLeaderTupper->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropLeaderTmelt->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            // phase change constant
            std::shared_ptr< fem::Property > tPropLeaderPCconst = std::make_shared< fem::Property >();
            tPropLeaderPCconst->set_parameters( {{{ 2.2 }}} );
            //tPropLeaderPCconst->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropLeaderPCconst->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            // phase state function type
            std::shared_ptr< fem::Property > tPropLeaderPCfunction = std::make_shared< fem::Property >();
            tPropLeaderPCfunction->set_parameters( {{{ 2 }}} );
            //tPropLeaderPCfunction->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropLeaderPCfunction->set_val_function( tConstValFunction_UTIWGGGLSDIFFBULK );

            // temperature load
            std::shared_ptr< fem::Property > tPropLeaderBodyLoad = nullptr;

            // define constitutive model ---------------------------------------------------------------- //
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiffLinIsoPC = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO_PC );
            tCMLeaderDiffLinIsoPC->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMLeaderDiffLinIsoPC->set_property( tPropLeaderConductivity, "Conductivity" );
            tCMLeaderDiffLinIsoPC->set_property( tPropLeaderDensity     , "Density" );
            tCMLeaderDiffLinIsoPC->set_property( tPropLeaderHeatCapacity, "HeatCapacity" );
            tCMLeaderDiffLinIsoPC->set_property( tPropLeaderLatentHeat  , "LatentHeat" );
            tCMLeaderDiffLinIsoPC->set_property( tPropLeaderTmelt       , "PCTemp" );
            tCMLeaderDiffLinIsoPC->set_property( tPropLeaderPCfunction  , "PhaseStateFunction" );
            tCMLeaderDiffLinIsoPC->set_property( tPropLeaderPCconst     , "PhaseChangeConst" );
            tCMLeaderDiffLinIsoPC->set_space_dim( 3 );
            tCMLeaderDiffLinIsoPC->set_local_properties();

            // define stabilization parameter ----------------------------------------------------------- //
            fem::SP_Factory tSPFactory;

            std::shared_ptr< fem::Stabilization_Parameter > tSPGGLSParam = tSPFactory.create_SP( fem::Stabilization_Type::GGLS_DIFFUSION );
            tSPGGLSParam->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::LEADER );
            tSPGGLSParam->set_parameters( { {{ 1.0 }} });
            tSPGGLSParam->set_property( tPropLeaderConductivity, "Conductivity", mtk::Leader_Follower::LEADER );
            tSPGGLSParam->set_property( tPropLeaderDensity, "Density", mtk::Leader_Follower::LEADER );
            tSPGGLSParam->set_property( tPropLeaderHeatCapacity, "HeatCapacity", mtk::Leader_Follower::LEADER );

            // create a dummy fem cluster and set it to SP
            fem::Cluster * tCluster = new fem::Cluster();
            tSPGGLSParam->set_cluster( tCluster );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
            tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
            tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::LEADER );
            tIWG->set_constitutive_model( tCMLeaderDiffLinIsoPC, "Diffusion", mtk::Leader_Follower::LEADER );
            tIWG->set_stabilization_parameter( tSPGGLSParam, "GGLSParam");
            tIWG->set_property( tPropLeaderBodyLoad, "Load", mtk::Leader_Follower::LEADER );

            // space and time geometry interpolators
            //------------------------------------------------------------------------------

            // create a space time geometry interpolator
            Geometry_Interpolator tGI( aGIRule );

            // set the coefficients xHat, tHat
            tGI.set_coeff( aXHat, aTHat );

            // set the evaluation point
            tGI.set_space_time( aParamPoint );

            // field interpolators
            //------------------------------------------------------------------------------

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFIs( 1 );

            // create the field interpolator
            tFIs( 0 ) = new Field_Interpolator( 1, aFIRule, &tGI, { MSI::Dof_Type::TEMP } );

            // set the coefficients uHat
            tFIs( 0 )->set_coeff( aDOFHat );

            //set the evaluation point xi, tau
            tFIs( 0 )->set_space_time( aParamPoint );

            // set a fem set pointer
            MSI::Equation_Set * tSet = new fem::Set();
            static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            // set size and populate the set leader dof type map
            tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            int aInt = (aNumDOFs-1);

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 1 );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, aInt } };

            // set size and fill the set jacobian assembly map
            tIWG->mSet->mJacDofAssemblyMap.resize( 1 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, aInt } };

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( aNumDOFs, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( aNumDOFs, aNumDOFs, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};

            // create a field interpolator manager
            moris::Vector< moris::Vector< enum MSI::Dof_Type > > tDummy;
            Field_Interpolator_Manager tFIManager( tDummy, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // check evaluation of the residual for IWG
            //------------------------------------------------------------------------------
            // evaluate the residual
            tIWG->compute_residual( 1.0 );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tJacobian;
            Matrix< DDRMat > tJacobianFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                    tEpsilonRel,
                    1.0,
                    tJacobian,
                    tJacobianFD );

            // require check is true
            //REQUIRE( tCheckJacobian );
            tChecks(0) = tCheckJacobian;

            // debug
            // moris::Matrix<DDRMat> test1 = tJacobianFD-tJacobian;
            // real tMax = test1.max();
            // print( tJacobian,   "tJacobian" );
            // print( tJacobianFD, "tJacobianFD" );
            //print( test1, "JacobianDifference" );
            // std::cout << "Maximum difference = " << tMax << " \n" << std::flush;

            return tChecks;

        } // end test function

        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "IWG_GGLS_Diffusion_Phase_Change_HEX8", "[moris],[fem],[IWG_GGLS_Diffusion_Phase_Change_HEX8]" )
        {
            //create a quad4 space element
            Matrix< DDRMat > tXHat = {
                    { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 }};

            //create a line time element
            Matrix< DDRMat > tTHat( 2, 1 );
            tTHat( 0 ) = 1.0e-3;
            tTHat( 1 ) = 1.1e-3;

            //create a space geometry interpolation rule
            mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create an interpolation rule
            mtk::Interpolation_Rule tIPRule (
                    mtk::Geometry_Type::HEX,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0 = {
                    {4.7},{4.8},{4.9},{5.3},{5.2},{5.3},{5.4},{4.5},
                    {4.9},{5.2},{4.9},{5.3},{5.7},{5.1},{4.8},{4.8}};
            uint tNumDOFs = 16;
            Matrix< DDRMat > tParametricPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.4 }};

            // run test
            moris::Vector<bool> tChecks = test_IWG_Diffusion_Phase_Change_GGLS(
                    tXHat,
                    tTHat,
                    tGeomInterpRule,
                    tIPRule,
                    tUHat0,
                    tParametricPoint,
                    tNumDOFs,
                    3);

            // checks
            bool tCheckJacobian = tChecks(0);
            REQUIRE( tCheckJacobian );

        } // end TEST_CASE

        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "IWG_GGLS_Diffusion_Phase_Change_HEX27", "[moris],[fem],[IWG_GGLS_Diffusion_Phase_Change_HEX27]" )
        {
            //create a HEX27 space element
            Matrix< DDRMat > tXHat = {
                    { 0.0, 0.0, 0.0}, { 4.0, 0.0, 0.0}, { 4.0, 1.0, 0.0}, { 0.0, 1.0, 0.0},
                    { 0.0, 0.0, 3.0}, { 4.0, 0.0, 3.0}, { 4.0, 1.0, 3.0}, { 0.0, 1.0, 3.0},
                    { 2.0, 0.0, 0.0}, { 4.0, 0.5, 0.0}, { 2.0, 1.0, 0.0}, { 0.0, 0.5, 0.0},
                    { 0.0, 0.0, 1.5}, { 4.0, 0.0, 1.5}, { 4.0, 1.0, 1.5}, { 0.0, 1.0, 1.5},
                    { 2.0, 0.0, 3.0}, { 4.0, 0.5, 3.0}, { 2.0, 1.0, 3.0}, { 0.0, 0.5, 3.0},
                    { 2.0, 0.5, 1.5}, { 2.0, 0.5, 0.0}, { 2.0, 0.5, 3.0},
                    { 2.0, 0.0, 1.5}, { 4.0, 0.5, 1.5}, { 2.0, 1.0, 1.5}, { 0.0, 0.5, 1.5}};

            //create a line time element
            Matrix< DDRMat > tTHat( 3, 1 );
            tTHat( 0 ) = 1.00e-3;
            tTHat( 2 ) = 1.05e-3;
            tTHat( 1 ) = 1.10e-3;

            //create a space geometry interpolation rule
            mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC );

            // create an interpolation rule
            mtk::Interpolation_Rule tIPRule (
                    mtk::Geometry_Type::HEX,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC );

            // set coefficients for field interpolators
            Matrix< DDRMat > tDOFHat = {
                    {5.1},{5.2},{5.3},{5.4},{5.5},{4.6},{4.7},{5.4},{4.9},{5.1},{5.2},{5.3},{5.4},{5.5},{4.6},{4.7},{4.8},{5.1},{5.3},{5.2},{5.3},{5.4},{4.5},{4.6},{4.7},{4.8},{4.9},
                    {5.1},{5.2},{5.3},{5.3},{5.5},{4.8},{4.7},{4.8},{4.9},{5.3},{5.2},{5.3},{5.4},{4.5},{5.2},{4.7},{4.8},{4.9},{5.1},{5.4},{5.3},{4.6},{5.5},{4.9},{4.7},{4.8},{4.9},
                    {5.4},{5.2},{5.3},{5.4},{5.2},{5.2},{4.7},{5.1},{5.3},{5.1},{5.1},{5.1},{5.2},{5.3},{5.2},{4.7},{4.6},{4.5},{5.6},{5.2},{5.3},{5.4},{5.3},{4.5},{5.3},{5.2},{5.4}};

            uint tNumDOFs = 27 * 3;
            Matrix< DDRMat > tParametricPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.4 }};

            // run test
            moris::Vector<bool> tChecks = test_IWG_Diffusion_Phase_Change_GGLS(
                    tXHat,
                    tTHat,
                    tGeomInterpRule,
                    tIPRule,
                    tDOFHat,
                    tParametricPoint,
                    tNumDOFs,
                    3);

            // checks
            bool tCheckJacobian = tChecks(0);
            REQUIRE( tCheckJacobian );

        } // end TEST_CASE

    } // end namespace fem
} // end namespace moris

