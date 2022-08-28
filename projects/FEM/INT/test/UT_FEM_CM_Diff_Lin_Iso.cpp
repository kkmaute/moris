/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_CM_Diff_Lin_Iso.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp" //FEM/INT/src
#include "fn_FEM_Check.hpp"

#define protected public
#define private   public
#include "cl_FEM_Constitutive_Model.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#undef protected
#undef private

void tValFunctionCM_Diff_Lin_Iso(
        moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 )
                 + aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tConstValFunctionCM_Diff_Lin_Iso(
        moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tDerFunctionCM_Diff_Lin_Iso(
        moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

using namespace moris;
using namespace fem;

moris::Cell<bool> test_diffusion_constitutive_model(
        Matrix< DDRMat > aXHat,
        Matrix< DDRMat > aTHat,
        mtk::Interpolation_Rule aGeomInterpRule,
        mtk::Interpolation_Rule aIPRule,
        Matrix< DDRMat > aUHat0,
        Matrix< DDRMat > aParametricPoint,
        uint aSpatialDim = 2)
{
    // real for check
    real tEpsilonRel = 2.0E-6;

    // initialize cell of checks
    moris::Cell<bool> tChecks( 7, false );

    // create the properties --------------------------------------------------------------------- //
    std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property >();
    tPropMasterConductivity->set_parameters( {{{ 1.1}}, {{1.1 }}} );
    //            tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //            tPropMasterConductivity->set_val_function( tValFunctionCM_Diff_Lin_Iso );
    //            tPropMasterConductivity->set_dof_derivative_functions( { tDerFunctionCM_Diff_Lin_Iso } );
    tPropMasterConductivity->set_val_function( tConstValFunctionCM_Diff_Lin_Iso );

    std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
    tPropMasterDensity->set_parameters( {{{ 1.2 }}} );
    tPropMasterDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropMasterDensity->set_val_function( tValFunctionCM_Diff_Lin_Iso );
    tPropMasterDensity->set_dof_derivative_functions( { tDerFunctionCM_Diff_Lin_Iso } );
    //            tPropMasterDensity->set_val_function( tConstValFunctionCM_Diff_Lin_Iso );

    std::shared_ptr< fem::Property > tPropMasterHeatCapacity = std::make_shared< fem::Property >();
    tPropMasterHeatCapacity->set_parameters( {{{ 1.3 }}} );
    tPropMasterHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropMasterHeatCapacity->set_val_function( tValFunctionCM_Diff_Lin_Iso );
    tPropMasterHeatCapacity->set_dof_derivative_functions( { tDerFunctionCM_Diff_Lin_Iso } );
    //            tPropMasterHeatCapacity->set_val_function( tConstValFunctionCM_Diff_Lin_Iso );

    // define constitutive models ---------------------------------------------------------------- //
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMMasterDiffLinIso->set_property( tPropMasterConductivity, "Conductivity" );
    tCMMasterDiffLinIso->set_property( tPropMasterDensity, "Density" );
    tCMMasterDiffLinIso->set_property( tPropMasterHeatCapacity, "HeatCapacity" );
    tCMMasterDiffLinIso->set_space_dim( aSpatialDim );
    tCMMasterDiffLinIso->set_local_properties();

    //create a space and a time geometry interpolator
    Geometry_Interpolator tGI( aGeomInterpRule );

    //set the coefficients xHat, tHat
    tGI.set_coeff( aXHat, aTHat );
    tGI.set_space_time(aParametricPoint);

    // create a TEMP field interpolator
    Cell< Field_Interpolator* > tFIs( 1, nullptr );
    tFIs( 0 ) = new Field_Interpolator ( 1, aIPRule, & tGI, { MSI::Dof_Type::TEMP } );

    // set coefficients for field interpolators
    tFIs( 0 )->set_coeff( aUHat0 );
    tFIs( 0 )->set_space_time(aParametricPoint);

    // create a fem set
    //MSI::Equation_Set * tSet = new fem::Set();
    fem::Set tSet;

    // set size for the set EqnObjDofTypeList
    tSet.mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tSet.mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSet.mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // set size and populate the set master dof type map
    tSet.mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSet.mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // create a field interpolator manager
    Cell< Cell < MSI::Dof_Type > > tDofTypes = {{ MSI::Dof_Type::TEMP }};
    Field_Interpolator_Manager tFIManager( tDofTypes, &tSet );
    tFIManager.mFI = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;

    // set field interpolator manager
    tCMMasterDiffLinIso->set_field_interpolator_manager( &tFIManager );

    // check flux-------------------------------------------------------------------
    //------------------------------------------------------------------------------
    // evaluate the constitutive model flux

    Matrix< DDRMat > tFlux = tCMMasterDiffLinIso->flux();
    //print( tFlux, "tFlux");

    // evaluate the constitutive model flux derivative
    Matrix< DDRMat > tdFluxdDOF = tCMMasterDiffLinIso->dFluxdDOF( { MSI::Dof_Type::TEMP } );
    //print( tdFluxdDOF, "tdFluxdDOF");

    // evaluate the constitutive model stress derivative by FD
    Matrix< DDRMat > tdFluxdDOF_FD;
    tCMMasterDiffLinIso->eval_dFluxdDOF_FD( { MSI::Dof_Type::TEMP }, tdFluxdDOF_FD, 1E-6 );

    //check stress derivative
    bool tCheckdStress = fem::check( tdFluxdDOF, tdFluxdDOF_FD, tEpsilonRel );

    // set bool in check list
    tChecks( 0 ) = tCheckdStress;

//    // debug
//    moris::print(tdFluxdDOF, "tdFluxdDOF");
//    moris::print(tdFluxdDOF_FD, "tdFluxdDOF_FD");

    // check strain-----------------------------------------------------------------
    //------------------------------------------------------------------------------
    // evaluate the constitutive model strain
    Matrix< DDRMat > tStrain = tCMMasterDiffLinIso->strain();
    //print( tStrain, "tStrain");

    // evaluate the constitutive model strain derivative
    Matrix< DDRMat > tdStraindDOF = tCMMasterDiffLinIso->dStraindDOF( { MSI::Dof_Type::TEMP } );
    //print( tdStraindDOF, "tdStraindDOF" );

    // evaluate the constitutive model strain derivative by FD
    Matrix< DDRMat > tdStraindDOF_FD;
    tCMMasterDiffLinIso->eval_dStraindDOF_FD( { MSI::Dof_Type::TEMP }, tdStraindDOF_FD, 1E-6 );

    //check strain derivative
    bool tCheckdStrain = fem::check( tdStraindDOF, tdStraindDOF_FD, tEpsilonRel );

    //set bool in check list
    tChecks( 1 ) = tCheckdStrain;

//    // debug
//    moris::print(tdStraindDOF, "tdStraindDOF");
//    moris::print(tdStraindDOF_FD, "tdStraindDOF_FD");

    // check constitutive matrix----------------------------------------------------
    //------------------------------------------------------------------------------
    // evaluate the constitutive model constitutive matrix
    Matrix< DDRMat > tConst = tCMMasterDiffLinIso->constitutive();
    //print( tConst, "tConst");

    // check traction---------------------------------------------------------------
    //------------------------------------------------------------------------------
    // define a normal
    Matrix< DDRMat > tNormal;
    if (aSpatialDim == 2)
        tNormal = {{1.0},{0.0}};
    else
        tNormal = {{0.0},{1.0},{0.0}};

    // evaluate the constitutive model traction
    Matrix< DDRMat > tTraction = tCMMasterDiffLinIso->traction( tNormal );
    //print( tTraction, "tTraction");

    // evaluate the constitutive model traction derivative
    Matrix< DDRMat > tdTractiondDOF = tCMMasterDiffLinIso->dTractiondDOF( { MSI::Dof_Type::TEMP }, tNormal);
    //print( tdTractiondDOF, "tdTractiondDOF" );

    // check test traction----------------------------------------------------------
    //------------------------------------------------------------------------------
    // evaluate the constitutive model test traction
    Matrix< DDRMat > tTestTraction = tCMMasterDiffLinIso->testTraction( tNormal, { MSI::Dof_Type::TEMP } );
    //print( tTestTraction, "tTestTraction");

    // evaluate the constitutive model test traction derivative
    Matrix< DDRMat > tdTestTractiondDOF = tCMMasterDiffLinIso->dTestTractiondDOF( { MSI::Dof_Type::TEMP }, tNormal, { { 0.0 } }, { MSI::Dof_Type::TEMP } );
    //print( tdTestTractiondDOF, "tdTestTractiondDOF" );

    // check test strain------------------------------------------------------------
    //------------------------------------------------------------------------------
    // evaluate the constitutive model test strain
    Matrix< DDRMat > tTestStrain = tCMMasterDiffLinIso->testStrain();
    //print( tTestStrain, "tTestStrain");

    // check Energy ------------------------------------------------------------------
    //------------------------------------------------------------------------------

    // evaluate the constitutive model flux derivative
    Matrix< DDRMat > tdEnergydDOF = tCMMasterDiffLinIso->dEnergydDOF( { MSI::Dof_Type::TEMP } );
    //print( tdFluxdDOF, "tdFluxdDOF");

    // evaluate the constitutive model stress derivative by FD
    Matrix< DDRMat > tdEnergydDOF_FD;
    tCMMasterDiffLinIso->eval_dEnergydDOF_FD( { MSI::Dof_Type::TEMP }, tdEnergydDOF_FD, 1E-6 );

    //check stress derivative
    bool tCheckEnergy = fem::check( tdEnergydDOF, tdEnergydDOF_FD, tEpsilonRel );

    // set bool in check list
    tChecks( 2 ) = tCheckEnergy;

    // check EnergyDot ------------------------------------------------------------------
    //------------------------------------------------------------------------------

    // evaluate the constitutive model flux derivative
    Matrix< DDRMat > tdEnergyDotdDOF = tCMMasterDiffLinIso->dEnergyDotdDOF( { MSI::Dof_Type::TEMP } );
    //print( tdFluxdDOF, "tdFluxdDOF");

    // evaluate the constitutive model stress derivative by FD
    Matrix< DDRMat > tdEnergyDotdDOF_FD;
    tCMMasterDiffLinIso->eval_dEnergyDotdDOF_FD( { MSI::Dof_Type::TEMP }, tdEnergyDotdDOF_FD, 1E-6 );

    //check stress derivative
    bool tCheckEnergyDot = fem::check( tdEnergyDotdDOF, tdEnergyDotdDOF_FD, tEpsilonRel );

    // set bool in check list
    tChecks( 3 ) = tCheckEnergyDot;

//    // debug
//    moris::print(tdEnergyDotdDOF, "tdEnergyDotdDOF");
//    moris::print(tdEnergyDotdDOF_FD, "tdEnergyDotdDOF_FD");

    // check gradEnergy -----------------------------------------------------------------
    //------------------------------------------------------------------------------

    // evaluate the constitutive model flux derivative
    Matrix< DDRMat > tdGradEnergydDOF = tCMMasterDiffLinIso->dGradEnergydDOF( { MSI::Dof_Type::TEMP } );
    //print( tdFluxdDOF, "tdFluxdDOF");

    // evaluate the constitutive model stress derivative by FD
    Matrix< DDRMat > tdGradEnergydDOF_FD;
    tCMMasterDiffLinIso->eval_dGradEnergydDOF_FD( { MSI::Dof_Type::TEMP }, tdGradEnergydDOF_FD, 1E-6 );

    //check stress derivative
    bool tCheckGradH = fem::check( tdGradEnergydDOF, tdGradEnergydDOF_FD, tEpsilonRel );

    // set bool to check list
    tChecks( 4 ) = tCheckGradH;

//    // debug
//    moris::print(tdGradEnergydDOF, "tdGradEnergydDOF");
//    moris::print(tdGradEnergydDOF_FD, "tdGradEnergydDOF_FD");

    // check gradEnergyDot --------------------------------------------------------------
    //------------------------------------------------------------------------------

    // evaluate the constitutive model flux derivative
    Matrix< DDRMat > tdGradEnergyDotdDOF = tCMMasterDiffLinIso->dGradEnergyDotdDOF( { MSI::Dof_Type::TEMP } );
    //print( tdFluxdDOF, "tdFluxdDOF");

    // evaluate the constitutive model stress derivative by FD
    Matrix< DDRMat > tdGradEnergyDotdDOF_FD;
    tCMMasterDiffLinIso->eval_dGradEnergyDotdDOF_FD( { MSI::Dof_Type::TEMP }, tdGradEnergyDotdDOF_FD, 1E-6 );

    //check stress derivative
    bool tCheckGradEnergyDot = fem::check( tdGradEnergyDotdDOF, tdGradEnergyDotdDOF_FD, tEpsilonRel );

    // set bool in check list
    tChecks( 5 ) = tCheckGradEnergyDot;

//    // debug
//    moris::print(tdGradEnergyDotdDOF, "tdGradEnergyDotdDOF");
//    moris::print(tdGradEnergyDotdDOF_FD, "tdGradEnergyDotdDOF_FD");

    // check graddivflux -----------------------------------------------------------
    //------------------------------------------------------------------------------

    // evaluate the constitutive model strain derivative
    Matrix< DDRMat > tdGradDivFluxdDOF = tCMMasterDiffLinIso->dGradDivFluxdDOF( { MSI::Dof_Type::TEMP } );
    //print( tdStraindDOF, "tdStraindDOF" );

    // evaluate the constitutive model strain derivative by FD
    Matrix< DDRMat > tdGradDivFluxdDOF_FD;
    tCMMasterDiffLinIso->eval_dGradDivFluxdDOF_FD( { MSI::Dof_Type::TEMP }, tdGradDivFluxdDOF_FD, 1E-6 );

    //check strain derivative
    bool tCheckGradDivFlux = fem::check( tdGradDivFluxdDOF, tdGradDivFluxdDOF_FD, tEpsilonRel );

    // set bool to check list
    tChecks( 6 ) = tCheckGradDivFlux;

//    // debug
//    moris::print(tdGradDivFluxdDOF, "tdGradDivFluxdDOF");
//    moris::print(tdGradDivFluxdDOF_FD, "tdGradDivFluxdDOF_FD");

    // clean up
    //------------------------------------------------------------------------------
    tFIs.clear();

    // return cell of checks
    return tChecks;

}/* TEST Function */

// ------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------- //
TEST_CASE( "CM_Diff_Lin_Iso_QUAD4", "[moris],[fem],[CM_Diff_Lin_Iso_QUAD4]" )
{
    //create a quad4 space element
    Matrix< DDRMat > tXHat = {
            { 0.0, 0.0},
            { 1.0, 0.0},
            { 1.0, 1.0},
            { 0.0, 1.0}};

    //create a line time element
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 1.0e-3;
    tTHat( 1 ) = 1.1e-3;

    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule(
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create an interpolation rule
    mtk::Interpolation_Rule tIPRule (
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // set coefficients for field interpolators
    Matrix< DDRMat > tUHat0 = {{3.9},{4.4},{4.9},{4.2},{4.9},{5.4},{5.9},{6.0}};
    Matrix< DDRMat > tParametricPoint = {{-0.4}, { 0.1}, {-0.6}};

    // run test
    moris::Cell<bool> tChecks = test_diffusion_constitutive_model(
            tXHat,
            tTHat,
            tGeomInterpRule,
            tIPRule,
            tUHat0,
            tParametricPoint);

    // checks
    bool tCheckdStress       = tChecks(0);
    bool tCheckdStrain       = tChecks(1);
    bool tCheckEnergy        = tChecks(2);
    bool tCheckEnergyDot     = tChecks(3);
    bool tCheckGradH         = tChecks(4);
    bool tCheckGradEnergyDot = tChecks(5);
    bool tCheckGradDivFlux   = tChecks(6);

    REQUIRE( tCheckdStress );
    REQUIRE( tCheckdStrain );
    REQUIRE( tCheckEnergy );
    REQUIRE( tCheckEnergyDot );
    REQUIRE( tCheckGradH );
    REQUIRE( tCheckGradEnergyDot );
    REQUIRE( tCheckGradDivFlux );

}

// ------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------- //
TEST_CASE( "CM_Diff_Lin_Iso_HEX27", "[moris],[fem],[CM_Diff_Lin_Iso_HEX27]" )
{
    // set number of spatial dimensions
    uint tSpatialDims = 3;

    //create a quad4 space element
    Matrix< DDRMat > tXHat = {
            { 0.0, 0.0, 0.0}, { 1.0, 0.0, 0.0}, { 1.0, 1.0, 0.0}, { 0.0, 1.0, 0.0},
            { 0.0, 0.0, 1.0}, { 1.0, 0.0, 1.0}, { 1.0, 1.0, 1.0}, { 0.0, 1.0, 1.0},
            { 0.5, 0.0, 0.0}, { 1.0, 0.5, 0.0}, { 0.5, 1.0, 0.0}, { 0.0, 0.5, 0.0},
            { 0.0, 0.0, 0.5}, { 1.0, 0.0, 0.5}, { 1.0, 1.0, 0.5}, { 0.0, 1.0, 0.5},
            { 0.5, 0.0, 1.0}, { 1.0, 0.5, 1.0}, { 0.5, 1.0, 1.0}, { 0.0, 0.5, 1.0},
            { 0.5, 0.5, 0.5}, { 0.5, 0.5, 0.0}, { 0.5, 0.5, 1.0},
            { 0.5, 0.0, 0.5}, { 1.0, 0.5, 0.5}, { 0.5, 1.0, 0.5}, { 0.0, 0.5, 0.5}};

    //create a line time element
    Matrix< DDRMat > tTHat( 3, 1 );
    tTHat( 0 ) = 1.00e-3;
    tTHat( 2 ) = 1.05e-3;
    tTHat( 1 ) = 1.10e-3;

    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule(
            mtk::Geometry_Type::HEX,
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
    Matrix< DDRMat > tUHat0 = {
            {4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.4},{4.9},{4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.8},{4.1},{4.3},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.8},{4.9},
            {5.1},{5.2},{5.3},{5.3},{5.5},{5.6},{5.7},{5.8},{5.9},{5.3},{5.2},{5.3},{5.4},{5.5},{5.2},{5.7},{5.8},{5.9},{5.1},{5.4},{5.3},{5.6},{5.5},{5.9},{5.7},{5.8},{5.9},
            {6.4},{6.2},{6.3},{6.4},{6.2},{6.6},{6.7},{6.1},{6.9},{6.1},{6.1},{6.3},{6.4},{6.9},{6.8},{6.7},{6.6},{6.5},{6.6},{6.2},{6.3},{6.4},{6.5},{6.6},{6.7},{6.8},{6.9}};
    Matrix< DDRMat > tParametricPoint = {{-0.4}, { 0.1}, {-0.6}, {0.3}};

    // run test
    moris::Cell<bool> tChecks = test_diffusion_constitutive_model(
            tXHat,
            tTHat,
            tGeomInterpRule,
            tIPRule,
            tUHat0,
            tParametricPoint,
            tSpatialDims);

    // checks
    bool tCheckdStress       = tChecks(0);
    bool tCheckdStrain       = tChecks(1);
    bool tCheckEnergy        = tChecks(2);
    bool tCheckEnergyDot     = tChecks(3);
    bool tCheckGradH         = tChecks(4);
    bool tCheckGradEnergyDot = tChecks(5);
    bool tCheckGradDivFlux   = tChecks(6);

    REQUIRE( tCheckdStress );
    REQUIRE( tCheckdStrain );
    REQUIRE( tCheckEnergy );
    REQUIRE( tCheckEnergyDot );
    REQUIRE( tCheckGradH );
    REQUIRE( tCheckGradEnergyDot );
    REQUIRE( tCheckGradDivFlux );
}

// ------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------- //
TEST_CASE( "CM_Diff_Lin_Iso_QUAD16", "[moris],[fem],[CM_Diff_Lin_Iso_QUAD16]" )
{

    //create a quad4 space element
    Matrix< DDRMat > tXHat = {
            { 0.0, 0.0},
            { 3.0, 0.0},
            { 3.0, 3.0},
            { 0.0, 3.0},
            { 1.0, 0.0},
            { 2.0, 0.0},
            { 3.0, 1.0},
            { 3.0, 2.0},
            { 2.0, 3.0},
            { 1.0, 3.0},
            { 0.0, 2.0},
            { 0.0, 1.0},
            { 1.0, 1.0},
            { 2.0, 1.0},
            { 2.0, 2.0},
            { 1.0, 2.0}};

    //create a line time element
    Matrix< DDRMat > tTHat( 4, 1 );
    tTHat( 0 ) = 1.00e-3;
    tTHat( 2 ) = 1.05e-3;
    tTHat( 3 ) = 1.10e-3;
    tTHat( 1 ) = 1.15e-3;

    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::CUBIC,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::CUBIC );

    // create an interpolation rule
    mtk::Interpolation_Rule tIPRule (
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::CUBIC,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::CUBIC );

    // set coefficients for field interpolators
    Matrix< DDRMat > tUHat0 = {
            {4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.8},{4.9},{4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},
            {5.1},{5.2},{5.3},{5.4},{5.5},{5.6},{5.7},{5.8},{5.9},{5.1},{5.2},{5.3},{5.4},{5.5},{5.6},{5.7},
            {6.1},{6.2},{6.3},{6.4},{6.5},{6.6},{6.7},{6.8},{6.9},{6.1},{6.2},{6.3},{6.4},{6.5},{6.6},{6.7},
            {7.1},{7.2},{7.3},{7.4},{7.5},{7.6},{7.7},{7.8},{7.9},{7.1},{7.2},{7.3},{7.4},{7.5},{7.6},{7.7},};

    Matrix< DDRMat > tParametricPoint = {{ 0.8}, {-0.9}, { 0.2}};

    // run test
    moris::Cell<bool> tChecks = test_diffusion_constitutive_model(
            tXHat,
            tTHat,
            tGeomInterpRule,
            tIPRule,
            tUHat0,
            tParametricPoint);

    // checks
    bool tCheckdStress       = tChecks(0);
    bool tCheckdStrain       = tChecks(1);
    bool tCheckEnergy        = tChecks(2);
    bool tCheckEnergyDot     = tChecks(3);
    bool tCheckGradH         = tChecks(4);
    bool tCheckGradEnergyDot = tChecks(5);
    bool tCheckGradDivFlux   = tChecks(6);

    REQUIRE( tCheckdStress );
    REQUIRE( tCheckdStrain );
    REQUIRE( tCheckEnergy );
    REQUIRE( tCheckEnergyDot );
    REQUIRE( tCheckGradH );
    REQUIRE( tCheckGradEnergyDot );
    REQUIRE( tCheckGradDivFlux );
}

