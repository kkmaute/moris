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
    aPropMatrix = aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

namespace moris
{
    namespace fem
    {
        TEST_CASE( "CM_Diff_Lin_Iso_PC", "[moris],[fem],[CM_Diff_Lin_Iso_PC]" )
        {
            // real for check
            real tEpsilon = 1E-6;

            // create the properties --------------------------------------------------------------------- //

            // conductivity
            std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property >();
            tPropMasterConductivity->set_parameters( {{{ 1.0}}, {{1.0 }}} );
            tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterConductivity->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // density
            std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
            tPropMasterDensity->set_parameters( {{{ 1.0}}} );
            tPropMasterDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterDensity->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // heat capacity
            std::shared_ptr< fem::Property > tPropMasterHeatCapacity = std::make_shared< fem::Property >();
            tPropMasterHeatCapacity->set_parameters( {{{ 1.0}}} );
            tPropMasterHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterHeatCapacity->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // latent heat
            std::shared_ptr< fem::Property > tPropMasterLatentHeat = std::make_shared< fem::Property >();
            tPropMasterLatentHeat->set_parameters( {{{ 100.0}}} );
            tPropMasterLatentHeat->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterLatentHeat->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // lower phase change temp
            std::shared_ptr< fem::Property > tPropMasterTlower = std::make_shared< fem::Property >();
            tPropMasterTlower->set_parameters( {{{ 1.0 }}} );
            tPropMasterTlower->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterTlower->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // upper phase change temp
            std::shared_ptr< fem::Property > tPropMasterTupper = std::make_shared< fem::Property >();
            tPropMasterTupper->set_parameters( {{{ 2.0 }}} );
            tPropMasterTupper->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterTupper->set_val_function( tConstValFunction_UT_CM_Diff_PC );

            // phase change constant
            std::shared_ptr< fem::Property > tPropMasterPCconst = std::make_shared< fem::Property >();
            tPropMasterPCconst->set_parameters( {{{ 0.0 }}} );
            tPropMasterPCconst->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterPCconst->set_val_function( tConstValFunction_UT_CM_Diff_PC );

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
            tCMMasterDiffLinIsoPC->set_property( tPropMasterTlower      , "Lower_PC_Temp" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterTupper      , "Upper_PC_Temp" );
            tCMMasterDiffLinIsoPC->set_property( tPropMasterPCconst     , "Phase_Change_Const" );
            tCMMasterDiffLinIsoPC->set_space_dim( 2 );

            //create a quad4 space element
            Matrix< DDRMat > tXHat( 4, 2 );
            tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
            tXHat( 1, 0 ) = 1.0; tXHat( 1, 1 ) = 0.0;
            tXHat( 2, 0 ) = 1.0; tXHat( 2, 1 ) = 1.0;
            tXHat( 3, 0 ) = 0.0; tXHat( 3, 1 ) = 1.0;

            //create a line time element
            Matrix< DDRMat > tTHat( 2, 1 );
            tTHat( 0 ) = 0.0;
            tTHat( 1 ) = 5.0;

            //create a space geometry interpolation rule
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

            //create a space and a time geometry interpolator
            Geometry_Interpolator tGI( tGeomInterpRule );

            //set the coefficients xHat, tHat
            tGI.set_coeff( tXHat, tTHat );
            tGI.set_space_time({{0.0}, {0.0}, {-1.0}});

            // create an interpolation rule
            Interpolation_Rule tIPRule ( mtk::Geometry_Type::QUAD,
                                         Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR,
                                         Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // create a TEMP field interpolator
            Cell< Field_Interpolator* > tFIs( 1, nullptr );
            tFIs( 0 ) = new Field_Interpolator ( 1, tIPRule, & tGI, { MSI::Dof_Type::TEMP } );

            // set coefficients for field interpolators
            Matrix< DDRMat > tUHat0 = {{1.0},{1.0},{2.0},{2.0},{1.0},{1.0},{2.0},{2.0}};
            tFIs( 0 )->set_coeff( tUHat0 );
            tFIs( 0 )->set_space_time({{0.0}, {0.0}, {-1.0}});

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
            for ( uint iStress = 0; iStress < tdHdotdDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdHdotdDOF.n_cols(); jStress++ )
                {
                    tCheckHdot = tCheckHdot && ( tdHdotdDOF( iStress, jStress ) - tdHdotdDOF_FD( iStress, jStress ) < tEpsilon );
                }
            }
            REQUIRE( tCheckHdot );



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
            for ( uint iStress = 0; iStress < tdGradHdotdDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdGradHdotdDOF.n_cols(); jStress++ )
                {
                    tCheckGradHdot = tCheckGradHdot && ( tdGradHdotdDOF( iStress, jStress ) - tdGradHdotdDOF_FD( iStress, jStress ) < tEpsilon );
                }
            }
            REQUIRE( tCheckGradHdot );


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
            for ( uint iStress = 0; iStress < tdGradDivFluxdDOF.n_rows(); iStress++ )
            {
                for( uint jStress = 0; jStress < tdGradDivFluxdDOF.n_cols(); jStress++ )
                {
                    tCheckGradDivFlux = tCheckGradDivFlux && ( tdGradDivFluxdDOF( iStress, jStress ) - tdGradDivFluxdDOF_FD( iStress, jStress ) < tEpsilon );
                }
            }
            REQUIRE( tCheckGradDivFlux );

            // clean up
            //------------------------------------------------------------------------------
            tFIs.clear();

        }/* TEST_CASE */


    }/* namespace fem */
}/* namespace moris */
