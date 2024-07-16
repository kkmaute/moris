/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Diffusion_Bulk.cpp
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
// MTK/src
#include "cl_MTK_Enums.hpp"
// LINALG/src
#include "op_equal_equal.hpp"
// FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"

using namespace moris;
using namespace fem;

inline void
tConstValFunction_UTIWGDIFFBULK(
        moris::Matrix< moris::DDRMat >&           aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*   aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

inline void
tGeoValFunction_UTIWGDIFFBULK(
        moris::Matrix< moris::DDRMat >&           aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*   aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
}

inline void
tFIValFunction_UTIWGDIFFBULK(
        moris::Matrix< moris::DDRMat >&           aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*   aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

inline void
tFIDerFunction_UTIWGDIFFBULK(
        moris::Matrix< moris::DDRMat >&           aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*   aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

inline void
tFIValDvFunction_UTIWGDIFFBULK(
        moris::Matrix< moris::DDRMat >&           aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*   aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::gen::PDV_Type::DENSITY )->val();
}

inline void
tFIDerDvFunction_UTIWGDIFFBULK(
        moris::Matrix< moris::DDRMat >&           aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*   aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::gen::PDV_Type::DENSITY )->N();
}

inline Vector< bool >
test_IWG_Diffusion_Bulk(
        Matrix< DDRMat >        aXHat,
        Matrix< DDRMat >        aTHat,
        mtk::Interpolation_Rule aGIRule,
        mtk::Interpolation_Rule aFIRule,
        Matrix< DDRMat >        aDOFHat,
        Matrix< DDRMat >        aParamPoint,
        uint                    aNumDOFs,
        uint                    aSpatialDim = 2 )
{
    // initialize cell of checks
    Vector< bool > tChecks( 1, false );

    // define an epsilon environment
    real tEpsilon = 1.0E-4;

    // define a perturbation relative size
    real tPerturbation = 1.0E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderConductivity = std::make_shared< fem::Property >();
    tPropLeaderConductivity->set_parameters( { { { 1.2 } } } );
    tPropLeaderConductivity->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
    tPropLeaderConductivity->set_val_function( tFIValFunction_UTIWGDIFFBULK );
    tPropLeaderConductivity->set_dof_derivative_functions( { tFIDerFunction_UTIWGDIFFBULK } );
    //    tPropLeaderConductivity->set_val_function( tConstValFunction_UTIWGDIFFBULK );

    std::shared_ptr< fem::Property > tPropLeaderTempLoad = std::make_shared< fem::Property >();
    tPropLeaderTempLoad->set_parameters( { { { 2.0 } } } );
    tPropLeaderTempLoad->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
    tPropLeaderTempLoad->set_val_function( tFIValFunction_UTIWGDIFFBULK );
    tPropLeaderTempLoad->set_dof_derivative_functions( { tFIDerFunction_UTIWGDIFFBULK } );
    //    tPropLeaderTempLoad->set_val_function( tConstValFunction_UTIWGDIFFBULK );

    std::shared_ptr< fem::Property > tPropLeaderDensity = std::make_shared< fem::Property >();
    tPropLeaderDensity->set_parameters( { { { 3.0 } } } );
    tPropLeaderDensity->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
    tPropLeaderDensity->set_val_function( tFIValFunction_UTIWGDIFFBULK );
    tPropLeaderDensity->set_dof_derivative_functions( { tFIDerFunction_UTIWGDIFFBULK } );
    //    tPropLeaderDensity->set_val_function( tConstValFunction_UTIWGDIFFBULK );

    std::shared_ptr< fem::Property > tPropLeaderHeatCapacity = std::make_shared< fem::Property >();
    tPropLeaderHeatCapacity->set_parameters( { { { 4.0 } } } );
    tPropLeaderHeatCapacity->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
    tPropLeaderHeatCapacity->set_val_function( tFIValFunction_UTIWGDIFFBULK );
    tPropLeaderHeatCapacity->set_dof_derivative_functions( { tFIDerFunction_UTIWGDIFFBULK } );
    //    tPropLeaderHeatCapacity->set_val_function( tConstValFunction_UTIWGDIFFBULK );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiffLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMLeaderDiffLinIso->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
    tCMLeaderDiffLinIso->set_property( tPropLeaderConductivity, "Conductivity" );
    tCMLeaderDiffLinIso->set_property( tPropLeaderDensity, "Density" );
    tCMLeaderDiffLinIso->set_property( tPropLeaderHeatCapacity, "HeatCapacity" );
    tCMLeaderDiffLinIso->set_space_dim( 3 );
    tCMLeaderDiffLinIso->set_local_properties();

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    tIWG->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMLeaderDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropLeaderTempLoad, "Load", mtk::Leader_Follower::LEADER );

    // space and time geometry interpolators
    //------------------------------------------------------------------------------

    // create a space time geometry interpolator
    Geometry_Interpolator tGI( aGIRule );

    // create time coeff tHat
    //    Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

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

    // set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( aParamPoint );

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    int aInt = ( aNumDOFs - 1 );

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
    tIWG->mRequestedLeaderGlobalDofTypes = { { MSI::Dof_Type::TEMP } };

    // create a field interpolator manager
    Vector< Vector< enum MSI::Dof_Type > > tDummy;
    Field_Interpolator_Manager             tFIManager( tDummy, tSet );

    // populate the field interpolator manager
    tFIManager.mFI                     = tFIs;
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
            tEpsilon,
            1.0,
            tJacobian,
            tJacobianFD );

    // require check is true
    // REQUIRE( tCheckJacobian );
    tChecks( 0 ) = tCheckJacobian;

    // debug
    // moris::Matrix<DDRMat> test1 = tJacobianFD-tJacobian;
    // real tMax = test1.max();
    // print( tJacobian,   "tJacobian" );
    // print( tJacobianFD, "tJacobianFD" );
    // print( test1, "JacobianDifference" );
    // std::cout << "Maximum difference = " << tMax << " \n" << std::flush;

    return tChecks;

}    // end test function

// ------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------- //
TEST_CASE( "IWG_Diffusion_Bulk_HEX8", "[moris],[fem],[IWG_Diffusion_Bulk_HEX8]" )
{
    // create a quad4 space element
    Matrix< DDRMat > tXHat = {
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 1.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
        { 1.0, 0.0, 1.0 },
        { 1.0, 1.0, 1.0 },
        { 0.0, 1.0, 1.0 }
    };

    // create a line time element
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 1.0e-3;
    tTHat( 1 ) = 1.1e-3;

    // create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create an interpolation rule
    mtk::Interpolation_Rule tIPRule(
            mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // set coefficients for field interpolators
    Matrix< DDRMat > tUHat0 = {
        { 3.9 },    //
        { 4.4 },
        { 4.9 },
        { 4.2 },
        { 4.9 },
        { 5.4 },
        { 5.9 },
        { 6.0 },
        { 5.9 },
        { 4.4 },
        { 3.9 },
        { 2.2 },
        { 4.9 },
        { 6.4 },
        { 4.9 },
        { 7.0 }
    };
    uint             tNumDOFs         = 16;
    Matrix< DDRMat > tParametricPoint = { { 0.35 }, { -0.25 }, { 0.75 }, { 0.4 } };

    // run test
    Vector< bool > tChecks = test_IWG_Diffusion_Bulk(
            tXHat,
            tTHat,
            tGeomInterpRule,
            tIPRule,
            tUHat0,
            tParametricPoint,
            tNumDOFs,
            3 );

    // checks
    bool tCheckJacobian = tChecks( 0 );
    REQUIRE( tCheckJacobian );

}    // end TEST_CASE

// ------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------- //
TEST_CASE( "IWG_Diffusion_Bulk_HEX27", "[moris],[fem],[IWG_Diffusion_Bulk_HEX27]" )
{
    // create a quad4 space element
    Matrix< DDRMat > tXHat = {
        { 0.0, 0.0, 0.0 },    //
        { 1.0, 0.0, 0.0 },
        { 1.0, 1.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
        { 1.0, 0.0, 1.0 },
        { 1.0, 1.0, 1.0 },
        { 0.0, 1.0, 1.0 },
        { 0.5, 0.0, 0.0 },
        { 1.0, 0.5, 0.0 },
        { 0.5, 1.0, 0.0 },
        { 0.0, 0.5, 0.0 },
        { 0.0, 0.0, 0.5 },
        { 1.0, 0.0, 0.5 },
        { 1.0, 1.0, 0.5 },
        { 0.0, 1.0, 0.5 },
        { 0.5, 0.0, 1.0 },
        { 1.0, 0.5, 1.0 },
        { 0.5, 1.0, 1.0 },
        { 0.0, 0.5, 1.0 },
        { 0.5, 0.5, 0.5 },
        { 0.5, 0.5, 0.0 },
        { 0.5, 0.5, 1.0 },
        { 0.5, 0.0, 0.5 },
        { 1.0, 0.5, 0.5 },
        { 0.5, 1.0, 0.5 },
        { 0.0, 0.5, 0.5 }
    };

    // create a line time element
    Matrix< DDRMat > tTHat( 3, 1 );
    tTHat( 0 ) = 1.00e-3;
    tTHat( 2 ) = 1.05e-3;
    tTHat( 1 ) = 1.10e-3;

    // create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC );

    // create an interpolation rule
    mtk::Interpolation_Rule tIPRule(
            mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC );

    // set coefficients for field interpolators
    Matrix< DDRMat > tDOFHat = {
        { 4.1 },    //
        { 4.2 },
        { 4.3 },
        { 4.4 },
        { 4.5 },
        { 4.6 },
        { 4.7 },
        { 4.4 },
        { 4.9 },
        { 4.1 },
        { 4.2 },
        { 4.3 },
        { 4.4 },
        { 4.5 },
        { 4.6 },
        { 4.7 },
        { 4.8 },
        { 4.1 },
        { 4.3 },
        { 4.2 },
        { 4.3 },
        { 4.4 },
        { 4.5 },
        { 4.6 },
        { 4.7 },
        { 4.8 },
        { 4.9 },
        { 5.1 },
        { 5.2 },
        { 5.3 },
        { 5.3 },
        { 5.5 },
        { 5.6 },
        { 5.7 },
        { 5.8 },
        { 5.9 },
        { 5.3 },
        { 5.2 },
        { 5.3 },
        { 5.4 },
        { 5.5 },
        { 5.2 },
        { 5.7 },
        { 5.8 },
        { 5.9 },
        { 5.1 },
        { 5.4 },
        { 5.3 },
        { 5.6 },
        { 5.5 },
        { 5.9 },
        { 5.7 },
        { 5.8 },
        { 5.9 },
        { 6.4 },
        { 6.2 },
        { 6.3 },
        { 6.4 },
        { 6.2 },
        { 6.6 },
        { 6.7 },
        { 6.1 },
        { 6.9 },
        { 6.1 },
        { 6.1 },
        { 6.3 },
        { 6.4 },
        { 6.9 },
        { 6.8 },
        { 6.7 },
        { 6.6 },
        { 6.5 },
        { 6.6 },
        { 6.2 },
        { 6.3 },
        { 6.4 },
        { 6.5 },
        { 6.6 },
        { 6.7 },
        { 6.8 },
        { 6.9 }
    };

    uint             tNumDOFs         = 27 * 3;
    Matrix< DDRMat > tParametricPoint = { { 0.35 }, { -0.25 }, { 0.75 }, { 0.4 } };

    // run test
    Vector< bool > tChecks = test_IWG_Diffusion_Bulk(
            tXHat,
            tTHat,
            tGeomInterpRule,
            tIPRule,
            tDOFHat,
            tParametricPoint,
            tNumDOFs,
            3 );

    // checks
    bool tCheckJacobian = tChecks( 0 );
    REQUIRE( tCheckJacobian );

}    // end TEST_CASE

// ------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------- //
TEST_CASE( "IWG_Diffusion_Bulk_Geo_Prop", "[moris],[fem],[IWG_Diff_Bulk_Geo_Prop]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderConductivity = std::make_shared< fem::Property >();
    tPropLeaderConductivity->set_parameters( { { { 1.0 } } } );
    tPropLeaderConductivity->set_val_function( tGeoValFunction_UTIWGDIFFBULK );

    std::shared_ptr< fem::Property > tPropLeaderDensity = std::make_shared< fem::Property >();
    tPropLeaderDensity->set_parameters( { { { 0.0 } } } );
    tPropLeaderDensity->set_val_function( tConstValFunction_UTIWGDIFFBULK );

    std::shared_ptr< fem::Property > tPropLeaderHeatCapacity = std::make_shared< fem::Property >();
    tPropLeaderHeatCapacity->set_parameters( { { { 0.0 } } } );
    tPropLeaderHeatCapacity->set_val_function( tConstValFunction_UTIWGDIFFBULK );

    std::shared_ptr< fem::Property > tPropLeaderTempLoad = std::make_shared< fem::Property >();
    tPropLeaderTempLoad->set_parameters( { { { 1.0 } } } );
    tPropLeaderTempLoad->set_val_function( tGeoValFunction_UTIWGDIFFBULK );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiffLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMLeaderDiffLinIso->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
    tCMLeaderDiffLinIso->set_property( tPropLeaderConductivity, "Conductivity" );
    tCMLeaderDiffLinIso->set_property( tPropLeaderDensity, "Density" );
    tCMLeaderDiffLinIso->set_property( tPropLeaderHeatCapacity, "HeatCapacity" );
    tCMLeaderDiffLinIso->set_space_dim( 3 );
    tCMLeaderDiffLinIso->set_local_properties();

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    tIWG->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMLeaderDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropLeaderTempLoad, "Load", mtk::Leader_Follower::LEADER );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = { { 0.35 }, { -0.25 }, { 0.75 }, { 0.0 } };

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    mtk::Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI( tGIRule );

    // create space coeff xHat
    Matrix< DDRMat > tXHat = { { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 1.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
        { 1.0, 0.0, 1.0 },
        { 1.0, 1.0, 1.0 },
        { 0.0, 1.0, 1.0 } };

    // create time coeff tHat
    Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

    // set the coefficients xHat, tHat
    tGI.set_coeff( tXHat, tTHat );

    // set the evaluation point
    tGI.set_space_time( tParamPoint );

    // field interpolators
    //------------------------------------------------------------------------------
    // create a space time interpolation rule
    mtk::Interpolation_Rule tFIRule( mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create random coefficients
    arma::Mat< double > tMatrix;
    tMatrix.randu( 8, 1 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat = 10.0 * tMatrix;

    // create a cell of field interpolators for IWG
    Vector< Field_Interpolator* > tFIs( 1 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );

    // set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );

    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );

    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    tIWG->mSet->mResDofAssemblyMap.resize( 1 );
    tIWG->mSet->mJacDofAssemblyMap.resize( 1 );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 7 } };
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 7 } };

    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size( 8, 1, 0.0 );
    tIWG->mSet->mJacobian.set_size( 8, 8, 0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    tIWG->mRequestedLeaderGlobalDofTypes = { { MSI::Dof_Type::TEMP } };

    Vector< Vector< enum MSI::Dof_Type > > tDummy;
    Field_Interpolator_Manager             tFIManager( tDummy, tSet );

    tFIManager.mFI                     = tFIs;
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
    Matrix< DDRMat > tJacobians;
    Matrix< DDRMat > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
            tEpsilon,
            1.0,
            tJacobians,
            tJacobiansFD );
    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    tFIs.clear();

} /* END_TEST_CASE */

TEST_CASE( "IWG_Diffusion_Bulk_Dv_Prop", "[moris],[fem],[IWG_Diff_Bulk_Dv_Prop]" )
{
    //    // define an epsilon environment
    //    real tEpsilon = 1E-6;
    //
    //    // define aperturbation relative size
    //    real tPerturbation = 1E-6;
    //
    //    // create the properties
    //    std::shared_ptr< fem::Property > tPropLeaderConductivity = std::make_shared< fem::Property > ();
    //    tPropLeaderConductivity->set_parameters( { {{ 1.0 }} } );
    //    tPropLeaderConductivity->set_dv_type_list( {{ gen::PDV_Type::DENSITY0 }} );
    //    tPropLeaderConductivity->set_val_function( tFIValDvFunction_UTIWGDIFFBULK );
    //    tPropLeaderConductivity->set_dv_derivative_functions( { tFIDerDvFunction_UTIWGDIFFBULK } );
    //
    //    std::shared_ptr< fem::Property > tPropLeaderTempLoad = std::make_shared< fem::Property > ();
    //    tPropLeaderTempLoad->set_parameters( { {{ 1.0 }} } );
    //    tPropLeaderTempLoad->set_val_function( tGeoValFunction_UTIWGDIFFBULK );
    //
    //    // define constitutive models
    //    fem::CM_Factory tCMFactory;
    //
    //    std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    //    tCMLeaderDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tCMLeaderDiffLinIso->set_property( tPropLeaderConductivity, "Conductivity" );
    //    tCMLeaderDiffLinIso->set_space_dim( 3 );
    //
    //    // define the IWGs
    //    fem::IWG_Factory tIWGFactory;
    //
    //    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_BULK );
    //    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    //    tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Leader_Follower::LEADER );
    //    tIWG->set_constitutive_model( tCMLeaderDiffLinIso, "Diffusion", mtk::Leader_Follower::LEADER );
    //    tIWG->set_property( tPropLeaderTempLoad, "Load", mtk::Leader_Follower::LEADER );
    //
    //    // create evaluation point xi, tau
    //    //------------------------------------------------------------------------------
    //    Matrix< DDRMat > tParamPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.0 }};
    //
    //    // space and time geometry interpolators
    //    //------------------------------------------------------------------------------
    //    // create a space geometry interpolation rule
    //    mtk::Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
    //                                mtk::Interpolation_Type::LAGRANGE,
    //                                mtk::Interpolation_Order::LINEAR,
    //                                mtk::Interpolation_Type::LAGRANGE,
    //                                mtk::Interpolation_Order::LINEAR );
    //
    //    // create a space time geometry interpolator
    //    Geometry_Interpolator tGI( tGIRule );
    //
    //    // create space coeff xHat
    //    arma::Mat< double > tXMatrix;
    //    tXMatrix.randu( 8, 3 );
    //    Matrix< DDRMat > tXHat;
    //    tXHat = 10.0 * tXMatrix;
    //
    //    // create time coeff tHat
    //    Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};
    //
    //    // set the coefficients xHat, tHat
    //    tGI.set_coeff( tXHat, tTHat );
    //
    //    // set the evaluation point
    //    tGI.set_space_time( tParamPoint );
    //
    //    // field interpolators
    //    //------------------------------------------------------------------------------
    //    //create a space time interpolation rule
    //    mtk::Interpolation_Rule tFIRule ( mtk::Geometry_Type::HEX,
    //                                 mtk::Interpolation_Type::LAGRANGE,
    //                                 mtk::Interpolation_Order::LINEAR,
    //                                 mtk::Interpolation_Type::CONSTANT,
    //                                 mtk::Interpolation_Order::CONSTANT );
    //
    //    // create random coefficients
    //    arma::Mat< double > tMatrix;
    //    tMatrix.randu( 8, 1 );
    //    Matrix< DDRMat > tDOFHat;
    //    tDOFHat = 10.0 * tMatrix;
    //
    //    // create a cell of field interpolators for IWG
    //    Vector< Field_Interpolator* > tFIs( 1 );
    //
    //    // create the field interpolator
    //    tFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );
    //
    //    // set the coefficients uHat
    //    tFIs( 0 )->set_coeff( tDOFHat );
    //
    //    //set the evaluation point xi, tau
    //    tFIs( 0 )->set_space_time( tParamPoint );
    //
    //    // create random coefficients
    //    arma::Mat< double > tDvMatrix;
    //    tDvMatrix.randu( 8, 1 );
    //    Matrix< DDRMat > tDvHat;
    //    tDvHat = 10.0 * tDvMatrix;
    //
    //    // create a cell of dv field interpolators for IWG
    //    Vector< Field_Interpolator* > tDvFIs( 1 );
    //
    //    // create the field interpolator
    //    tDvFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { gen::PDV_Type::DENSITY0 } );
    //
    //    // set the coefficients
    //    tDvFIs( 0 )->set_coeff( tDvHat );
    //
    //    //set the evaluation point xi, tau
    //    tDvFIs( 0 )->set_space_time( tParamPoint );
    //
    //    // create a fem set
    //    MSI::Equation_Set * tSet = new fem::Set();
    //
    //    tIWG->set_set_pointer(static_cast<fem::Set*>(tSet));
    //
    //    tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );
    //    tIWG->mSet->mUniqueDvTypeList.resize( 5, gen::PDV_Type::END_ENUM );
    //
    //    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    //    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >(MSI::Dof_Type::TEMP) ) = 0;
    //
    //    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    //    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >(MSI::Dof_Type::TEMP) ) = 0;
    //
    //    tIWG->mSet->mUniqueDvTypeMap.set_size( static_cast< int >( gen::PDV_Type::END_ENUM ) + 1, 1, -1 );
    //    tIWG->mSet->mUniqueDvTypeMap( static_cast< int >( gen::PDV_Type::DENSITY0 ) ) = 0;
    //
    //    tIWG->mSet->mLeaderDvTypeMap.set_size( static_cast< int >( gen::PDV_Type::END_ENUM ) + 1, 1, -1 );
    //    tIWG->mSet->mLeaderDvTypeMap( static_cast< int >( gen::PDV_Type::DENSITY0 ) ) = 0;
    //
    //    tIWG->mSet->mResDofAssemblyMap.resize( 1 );
    //    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 7 } };
    //
    //    tIWG->mSet->mDvAssemblyMap.resize( 1 );
    //    tIWG->mSet->mDvAssemblyMap( 0 ) = { { 0, 7 } };
    //
    //    tIWG->mSet->mResidual.resize( 1 );
    //    tIWG->mSet->mResidual( 0 ).set_size( 8, 1, 0.0 );
    //
    //    // build global dof type list
    //    tIWG->get_global_dof_type_list();
    //    tIWG->get_global_dv_type_list();
    //
    //    tIWG->mRequestedLeaderGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};
    //
    //    Vector< Vector< enum MSI::Dof_Type > > tDummy;
    //    Field_Interpolator_Manager tFIManager( tDummy, tSet );
    //
    //    tFIManager.mFI = tFIs;
    //    tFIManager.mDvFI = tDvFIs;
    //    tFIManager.mIPGeometryInterpolator = &tGI;
    //    tFIManager.mIGGeometryInterpolator = &tGI;
    //
    //    // set the interpolator manager to the set
    //    tIWG->mSet->mLeaderFIManager = &tFIManager;
    //
    //    // set IWG field interpolator manager
    //    tIWG->set_field_interpolator_manager( &tFIManager );
    //
    //    // check evaluation of drdpdv  by FD
    //    //------------------------------------------------------------------------------
    //    // init the jacobian for IWG and FD evaluation
    //    Vector< Matrix< DDRMat > > tdRdpMatFD;
    //    Vector< Matrix< DDRMat > > tdRdpGeoFD;
    //    Vector< Matrix< DDSMat > > tIsActive = { {{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 }},
    //                                           {{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 }},
    //                                           {{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 },{ 1 }} };
    //
    //    // check jacobian by FD
    //    tIWG->compute_dRdp_FD_material( 1.0,
    //                                    tPerturbation,
    //                                    tdRdpMatFD );
    //
    //    tIWG->compute_dRdp_FD_geometry( 1.0,
    //                                    tPerturbation,
    //                                    tIsActive,
    //                                    tdRdpGeoFD );
    //
    //    // print for debug
    //    //print( tdrdpdvMatFD( 0 ), "tdrdpdvMatFD" );
    //    //print( tdrdpdvGeoFD( 0 ), "tdrdpdvGeoFD" );
    //
    //    // clean up
    //    tFIs.clear();
    //    tDvFIs.clear();

} /* END_TEST_CASE */
