#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_CM_Factory.hpp" //FEM/INT/src

#define protected public
#define private   public
#include "cl_FEM_Constitutive_Model.hpp" //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#undef protected
#undef private

void tValFunctionCM_Diff_Lin_Iso
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 )
         + aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tConstValFunction_UT_CM_Diff_PC
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tDerFunctionCM_Diff_Lin_Iso
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

namespace moris
{
    namespace fem
    {
        moris::Cell<bool> test_phase_change_constitutive_model(
                Matrix< DDRMat > aXHat,
                Matrix< DDRMat > aTHat,
                Interpolation_Rule aGeomInterpRule,
                Interpolation_Rule aIPRule,
                Matrix< DDRMat > aUHat0,
                Matrix< DDRMat > aParametricPoint,
                uint aSpatialDim = 2)
        {
            // initialize cell of checks
            moris::Cell<bool> tChecks( 3, false );

            // real for check
            real tEpsilonRel = 5E-5;

            // create the properties --------------------------------------------------------------------- //

            // conductivity
            std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property >();
            tPropMasterConductivity->set_parameters( {{{ 1.1 }}, {{ 1.1 }}} );
//            tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
//            tPropMasterConductivity->set_val_function( tValFunctionCM_Diff_Lin_Iso );
//            tPropMasterConductivity->set_dof_derivative_functions( { tDerFunctionCM_Diff_Lin_Iso } );
            tPropMasterConductivity->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // density
            std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
            tPropMasterDensity->set_parameters( {{{ 1.2 }}} );
            tPropMasterDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterDensity->set_val_function( tValFunctionCM_Diff_Lin_Iso );
            tPropMasterDensity->set_dof_derivative_functions( { tDerFunctionCM_Diff_Lin_Iso } );
//            tPropMasterDensity->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // heat capacity
            std::shared_ptr< fem::Property > tPropMasterHeatCapacity = std::make_shared< fem::Property >();
            tPropMasterHeatCapacity->set_parameters( {{{ 1.3 }}} );
            tPropMasterHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterHeatCapacity->set_val_function( tValFunctionCM_Diff_Lin_Iso );
            tPropMasterHeatCapacity->set_dof_derivative_functions( { tDerFunctionCM_Diff_Lin_Iso } );
//            tPropMasterHeatCapacity->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // latent heat
            std::shared_ptr< fem::Property > tPropMasterLatentHeat = std::make_shared< fem::Property >();
            tPropMasterLatentHeat->set_parameters( {{{ 100.0}}} );
            //tPropMasterLatentHeat->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterLatentHeat->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // phase change temp
            std::shared_ptr< fem::Property > tPropMasterTmelt = std::make_shared< fem::Property >();
            tPropMasterTmelt->set_parameters( {{{ 5.0 }}} );
            //tPropMasterTupper->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterTmelt->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // phase change constant
            std::shared_ptr< fem::Property > tPropMasterPCconst = std::make_shared< fem::Property >();
            tPropMasterPCconst->set_parameters( {{{ 1.7 }}} );
            //tPropMasterPCconst->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterPCconst->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // phase state function type
            std::shared_ptr< fem::Property > tPropMasterPCfunction = std::make_shared< fem::Property >();
            tPropMasterPCfunction->set_parameters( {{{ 1 }}} );
            //tPropMasterPCfunction->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterPCfunction->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // temperature load
            std::shared_ptr< fem::Property > tPropMasterBodyLoad = nullptr;

            // define constitutive models ---------------------------------------------------------------- //
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIsoPC = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO_PC );
            tCMMasterDiffLinIsoPC->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterConductivity, "Conductivity" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterDensity     , "Density" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterHeatCapacity, "Heat_Capacity" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterLatentHeat  , "Latent_Heat" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterTmelt       , "PC_Temp" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterPCfunction  , "Phase_State_Function" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterPCconst     , "Phase_Change_Const" );
            tCMMasterDiffLinIsoPC->set_space_dim( aSpatialDim );

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
            tCMMasterDiffLinIsoPC->set_field_interpolator_manager( &tFIManager );

            // check Hdot ------------------------------------------------------------------
            //------------------------------------------------------------------------------

            // evaluate the constitutive model flux derivative
            Matrix< DDRMat > tdHdotdDOF = tCMMasterDiffLinIsoPC->dHdotdDOF( { MSI::Dof_Type::TEMP } );
            //print( tdFluxdDOF, "tdFluxdDOF");

            // evaluate the constitutive model stress derivative by FD
            Matrix< DDRMat > tdHdotdDOF_FD;
            tCMMasterDiffLinIsoPC->eval_dHdotdDOF_FD( { MSI::Dof_Type::TEMP }, tdHdotdDOF_FD, 1E-6 );

            //check stress derivative
            bool tCheckHdot = true;
            real tEpsilon = std::abs( tEpsilonRel * tdHdotdDOF( 0, 0 ) );
            for ( uint iStress = 0; iStress < tdHdotdDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdHdotdDOF.n_cols(); jStress++ )
                {
                    tCheckHdot = tCheckHdot && ( std::abs( tdHdotdDOF( iStress, jStress ) - tdHdotdDOF_FD( iStress, jStress ) ) < tEpsilon );
                }
            }
            //REQUIRE( tCheckHdot );
            tChecks(0) = tCheckHdot;

// debug
//moris::print(tdHdotdDOF, "tdHdotdDOF");
//moris::print(tdHdotdDOF_FD, "tdHdotdDOF_FD");


            // check gradHdot --------------------------------------------------------------
            //------------------------------------------------------------------------------

            // evaluate the constitutive model flux derivative
            Matrix< DDRMat > tdGradHdotdDOF = tCMMasterDiffLinIsoPC->dGradHdotdDOF( { MSI::Dof_Type::TEMP } );
            //print( tdFluxdDOF, "tdFluxdDOF");

            // evaluate the constitutive model stress derivative by FD
            Matrix< DDRMat > tdGradHdotdDOF_FD;
            tCMMasterDiffLinIsoPC->eval_dGradHdotdDOF_FD( { MSI::Dof_Type::TEMP }, tdGradHdotdDOF_FD, 1E-6 );

            //check stress derivative
            bool tCheckGradHdot = true;
            tEpsilon = std::abs( tEpsilonRel * tdGradHdotdDOF( 0, 0 ) );
            for ( uint iStress = 0; iStress < tdGradHdotdDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdGradHdotdDOF.n_cols(); jStress++ )
                {
                    tCheckGradHdot = tCheckGradHdot && ( std::abs( tdGradHdotdDOF( iStress, jStress ) - tdGradHdotdDOF_FD( iStress, jStress ) ) < tEpsilon );

// debug
//if ( std::abs( tdGradHdotdDOF( iStress, jStress ) - tdGradHdotdDOF_FD( iStress, jStress ) ) > tEpsilon )
//{
//    std::cout << "tdGradHdotdDOF: failed Jacobian check at: " << iStress << "x" << jStress
//            << " with difference: " << std::abs( tdGradHdotdDOF( iStress, jStress ) - tdGradHdotdDOF_FD( iStress, jStress ) )
//            << " , allowed tolerance: " << tEpsilon << "\n" << std::flush;
//}

                }
            }
            //REQUIRE( tCheckGradHdot );
            tChecks(1) = tCheckGradHdot;

// debug
//moris::print(tdGradHdotdDOF, "tdGradHdotdDOF");
//moris::print(tdGradHdotdDOF_FD, "tdGradHdotdDOF_FD");


            // check graddivflux -----------------------------------------------------------
            //------------------------------------------------------------------------------

            // evaluate the constitutive model strain derivative
            Matrix< DDRMat > tdGradDivFluxdDOF = tCMMasterDiffLinIsoPC->dGradDivFluxdDOF( { MSI::Dof_Type::TEMP } );
            //print( tdStraindDOF, "tdStraindDOF" );

            // evaluate the constitutive model strain derivative by FD
            Matrix< DDRMat > tdGradDivFluxdDOF_FD;
            tCMMasterDiffLinIsoPC->eval_dGradDivFluxdDOF_FD( { MSI::Dof_Type::TEMP }, tdGradDivFluxdDOF_FD, 1E-6 );

            //check strain derivative
            bool tCheckGradDivFlux = true;
            tEpsilon = std::abs( tEpsilonRel );
            for ( uint iStress = 0; iStress < tdGradDivFluxdDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdGradDivFluxdDOF.n_cols(); jStress++ )
                {
                    tCheckGradDivFlux = tCheckGradDivFlux && ( std::abs( tdGradDivFluxdDOF( iStress, jStress ) - tdGradDivFluxdDOF_FD( iStress, jStress ) ) < tEpsilon );
                }
            }
            //REQUIRE( tCheckGradDivFlux );
            tChecks(2) = tCheckGradDivFlux;
// debug
//moris::print(tdGradDivFluxdDOF, "tdGradDivFluxdDOF");
//moris::print(tdGradDivFluxdDOF_FD, "tdGradDivFluxdDOF_FD");


            //------------------------------------------------------------------------------

            // clean up
            tFIs.clear();

            // return cell of checks
            return tChecks;

        }/* TEST Function */


        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "CM_Diff_Lin_Iso_PC_QUAD4", "[moris],[fem],[CM_Diff_Lin_Iso_PC_QUAD4]" )
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
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create an interpolation rule
            Interpolation_Rule tIPRule (
                    mtk::Geometry_Type::QUAD,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0 = {{3.9},{4.4},{4.9},{4.2},{4.9},{5.4},{5.9},{6.0}};
            Matrix< DDRMat > tParametricPoint = {{-0.4}, { 0.1}, {-0.6}};

            // run test
            moris::Cell<bool> tChecks = test_phase_change_constitutive_model(
                            tXHat,
                            tTHat,
                            tGeomInterpRule,
                            tIPRule,
                            tUHat0,
                            tParametricPoint);

            // checks
            bool tCheckHdot = tChecks(0);
            bool tCheckGradHdot = tChecks(1);
            bool tCheckGradDivFlux = tChecks(2);
            REQUIRE( tCheckHdot );
            REQUIRE( tCheckGradHdot );
            REQUIRE( tCheckGradDivFlux );

        }

        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "CM_Diff_Lin_Iso_PC_HEX27", "[moris],[fem],[CM_Diff_Lin_Iso_PC_HEX27]" )
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
            Interpolation_Rule tGeomInterpRule(
                    mtk::Geometry_Type::HEX,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC );

            // create an interpolation rule
            Interpolation_Rule tIPRule (
                    mtk::Geometry_Type::HEX,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC );


            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0 = {
                    {4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.4},{4.9},{4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.8},{4.1},{4.3},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.8},{4.9},
                    {5.1},{5.2},{5.3},{5.3},{5.5},{5.6},{5.7},{5.8},{5.9},{5.3},{5.2},{5.3},{5.4},{5.5},{5.2},{5.7},{5.8},{5.9},{5.1},{5.4},{5.3},{5.6},{5.5},{5.9},{5.7},{5.8},{5.9},
                    {6.4},{6.2},{6.3},{6.4},{6.2},{6.6},{6.7},{6.1},{6.9},{6.1},{6.1},{6.3},{6.4},{6.9},{6.8},{6.7},{6.6},{6.5},{6.6},{6.2},{6.3},{6.4},{6.5},{6.6},{6.7},{6.8},{6.9}};
            Matrix< DDRMat > tParametricPoint = {{-0.4}, { 0.1}, {-0.6}, {0.3}};

            // run test
            moris::Cell<bool> tChecks = test_phase_change_constitutive_model(
                            tXHat,
                            tTHat,
                            tGeomInterpRule,
                            tIPRule,
                            tUHat0,
                            tParametricPoint,
                            tSpatialDims);

            // checks
            bool tCheckHdot = tChecks(0);
            bool tCheckGradHdot = tChecks(1);
            bool tCheckGradDivFlux = tChecks(2);
            REQUIRE( tCheckHdot );
            REQUIRE( tCheckGradHdot );
            REQUIRE( tCheckGradDivFlux );

        }

        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "CM_Diff_Lin_Iso_PC_QUAD16", "[moris],[fem],[CM_Diff_Lin_Iso_PC_QUAD16]" )
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
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::CUBIC,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::CUBIC );

            // create an interpolation rule
            Interpolation_Rule tIPRule (
                    mtk::Geometry_Type::QUAD,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::CUBIC,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::CUBIC );


            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0 = {
                    {4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.8},{4.9},{4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},
                    {5.1},{5.2},{5.3},{5.4},{5.5},{5.6},{5.7},{5.8},{5.9},{5.1},{5.2},{5.3},{5.4},{5.5},{5.6},{5.7},
                    {6.1},{6.2},{6.3},{6.4},{6.5},{6.6},{6.7},{6.8},{6.9},{6.1},{6.2},{6.3},{6.4},{6.5},{6.6},{6.7},
                    {7.1},{7.2},{7.3},{7.4},{7.5},{7.6},{7.7},{7.8},{7.9},{7.1},{7.2},{7.3},{7.4},{7.5},{7.6},{7.7},};
            Matrix< DDRMat > tParametricPoint = {{-0.4}, { 0.1}, {-0.6}};

            // run test
            moris::Cell<bool> tChecks = test_phase_change_constitutive_model(
                            tXHat,
                            tTHat,
                            tGeomInterpRule,
                            tIPRule,
                            tUHat0,
                            tParametricPoint);

            // checks
            bool tCheckHdot = tChecks(0);
            bool tCheckGradHdot = tChecks(1);
            bool tCheckGradDivFlux = tChecks(2);
            REQUIRE( tCheckHdot );
            REQUIRE( tCheckGradHdot );
            REQUIRE( tCheckGradDivFlux );

        }


    }/* namespace fem */
}/* namespace moris */
