
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
         + aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tDerFunctionCM_Diff_Lin_Iso
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

namespace moris
{
    namespace fem
    {
        moris::Cell<bool> test_diffusion_constitutive_model(
                Matrix< DDRMat > aXHat,
                Matrix< DDRMat > aTHat,
                Interpolation_Rule aGeomInterpRule,
                Interpolation_Rule aIPRule,
                Matrix< DDRMat > aUHat0,
                Matrix< DDRMat > aParametricPoint)
        {
            // real for check
            real tEpsilon = 1E-6;

            // initialize cell of checks
            moris::Cell<bool> tChecks( 2, false );

            // real for check
            real tEpsilonRel = 1E-6;

            // create the properties --------------------------------------------------------------------- //
            std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property >();
            tPropMasterConductivity->set_parameters( {{{ 1.0}}, {{1.0 }}} );
            tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterConductivity->set_val_function( tValFunctionCM_Diff_Lin_Iso );
            tPropMasterConductivity->set_dof_derivative_functions( { tDerFunctionCM_Diff_Lin_Iso } );

            // define constitutive models ---------------------------------------------------------------- //
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMMasterDiffLinIso->set_property( tPropMasterConductivity, "Conductivity" );
            tCMMasterDiffLinIso->set_space_dim( 2 );

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
            bool tCheckdStress = true;
            for ( uint iStress = 0; iStress < tdFluxdDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdFluxdDOF.n_cols(); jStress++ )
                {
                    tCheckdStress = tCheckdStress && ( tdFluxdDOF( iStress, jStress ) - tdFluxdDOF_FD( iStress, jStress ) < tEpsilon );
                }
            }
            //REQUIRE( tCheckdStress );
            tChecks(0) = tCheckdStress;

// debug
//moris::print(tdFluxdDOF, "tdFluxdDOF");
//moris::print(tdFluxdDOF_FD, "tdFluxdDOF_FD");


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
            bool tCheckdStrain = true;
            for ( uint iStress = 0; iStress < tdStraindDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdStraindDOF.n_cols(); jStress++ )
                {
                    tCheckdStrain = tCheckdStrain && ( tdStraindDOF( iStress, jStress ) - tdStraindDOF_FD( iStress, jStress ) < tEpsilon );
                }
            }
            //REQUIRE( tCheckdStrain );
            tChecks(1) = tCheckdStrain;

// debug
//moris::print(tdStraindDOF, "tdStraindDOF");
//moris::print(tdStraindDOF_FD, "tdStraindDOF_FD");


            // check constitutive matrix----------------------------------------------------
            //------------------------------------------------------------------------------
            // evaluate the constitutive model constitutive matrix
            Matrix< DDRMat > tConst = tCMMasterDiffLinIso->constitutive();
            //print( tConst, "tConst");

            // check traction---------------------------------------------------------------
            //------------------------------------------------------------------------------
            // define a normal
            Matrix< DDRMat > tNormal = {{1.0},{0.0}};

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

            // clean up
            //------------------------------------------------------------------------------
            tFIs.clear();

            // return cell of checks
            return tChecks;

        }/* TEST Function */



        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "CM_Diff_Lin_Iso_Linear", "[moris],[fem],[CM_Diff_Lin_Iso_Linear]" )
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
            Interpolation_Rule tGeomInterpRule(
                    mtk::Geometry_Type::QUAD,
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
            moris::Cell<bool> tChecks = test_diffusion_constitutive_model(
                            tXHat,
                            tTHat,
                            tGeomInterpRule,
                            tIPRule,
                            tUHat0,
                            tParametricPoint);

            // checks
            bool tCheckdStress = tChecks(0);
            bool tCheckdStrain = tChecks(1);
            REQUIRE( tCheckdStress );
            REQUIRE( tCheckdStrain );

        }

        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "CM_Diff_Lin_Iso_Quadratic", "[moris],[fem],[CM_Diff_Lin_Iso_Quadratic]" )
        {
            //create a quad4 space element
            Matrix< DDRMat > tXHat = {
                    { 0.0, 0.0},
                    { 1.0, 0.0},
                    { 1.0, 1.0},
                    { 0.0, 1.0},
                    { 0.5, 0.0},
                    { 1.0, 0.5},
                    { 0.5, 1.0},
                    { 0.0, 0.5},
                    { 0.5, 0.5}};

            //create a line time element
            Matrix< DDRMat > tTHat( 3, 1 );
            tTHat( 0 ) = 1.00e-3;
            tTHat( 2 ) = 1.05e-3;
            tTHat( 1 ) = 1.10e-3;

            //create a space geometry interpolation rule
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC );

            // create an interpolation rule
            Interpolation_Rule tIPRule (
                    mtk::Geometry_Type::QUAD,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::QUADRATIC );


            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0 = {
                    {4.1},{4.2},{4.3},{4.4},{4.5},{4.6},{4.7},{4.8},{4.9},
                    {5.1},{5.2},{5.3},{5.4},{5.5},{5.6},{5.7},{5.8},{5.9},
                    {6.1},{6.2},{6.3},{6.4},{6.5},{6.6},{6.7},{6.8},{6.9},};
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
            bool tCheckdStress = tChecks(0);
            bool tCheckdStrain = tChecks(1);
            REQUIRE( tCheckdStress );
            REQUIRE( tCheckdStrain );


        }

        // ------------------------------------------------------------------------------------- //
        // ------------------------------------------------------------------------------------- //
        TEST_CASE( "CM_Diff_Lin_Iso_Cubic", "[moris],[fem],[CM_Diff_Lin_Iso_Cubic]" )
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
            moris::Cell<bool> tChecks = test_diffusion_constitutive_model(
                            tXHat,
                            tTHat,
                            tGeomInterpRule,
                            tIPRule,
                            tUHat0,
                            tParametricPoint);

            // checks
            bool tCheckdStress = tChecks(0);
            bool tCheckdStrain = tChecks(1);
            REQUIRE( tCheckdStress );
            REQUIRE( tCheckdStrain );
        }


    }/* namespace fem */
}/* namespace moris */
