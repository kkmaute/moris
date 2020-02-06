
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

//moris::Matrix< moris::DDRMat > tValFunctionCM_Diff_Lin_Iso( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
//                                                            moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
//                                                            moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
//                                                            moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
//{
//    return aParameters( 0 ) + aParameters( 1 ) * aDofFI( 0 )->val();
//}
//
//moris::Matrix< moris::DDRMat > tDerFunctionCM_Diff_Lin_Iso( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
//                                                            moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
//                                                            moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
//                                                            moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
//{
//    return aParameters( 1 ) * aDofFI( 0 )->N();
//}

moris::Matrix< moris::DDRMat > tValFunctionCM_Diff_Lin_Iso
( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager *         aFIManager )
{
    return aParameters( 0 )
         + aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

moris::Matrix< moris::DDRMat > tDerFunctionCM_Diff_Lin_Iso
( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager *         aFIManager )
{
    return aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

namespace moris
{
    namespace fem
    {
        TEST_CASE( "CM_Diff_Lin_Iso", "[moris],[fem],[CM_Diff_Lin_Iso]" )
        {
            // real for check
            real tEpsilon = 1E-6;

            // create the properties
            std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property >();
            tPropMasterConductivity->set_parameters( {{{ 1.0}}, {{1.0 }}} );
            tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tPropMasterConductivity->set_val_function( tValFunctionCM_Diff_Lin_Iso );
            tPropMasterConductivity->set_dof_derivative_functions( { tDerFunctionCM_Diff_Lin_Iso } );

            // define constitutive models
            fem::CM_Factory tCMFactory;

            std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
            tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
            tCMMasterDiffLinIso->set_property( tPropMasterConductivity, "Conductivity" );
            tCMMasterDiffLinIso->set_space_dim( 2 );

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

            // create a set
            //MSI::Equation_Set * tSet = new fem::Set();
            fem::Set tSet;

            // set size for the set EqnObjDofTypeList
            tSet.mEqnObjDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tSet.mDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tSet.mDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            // set size and populate the set master dof type map
            tSet.mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tSet.mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

            // create a field interpolator manager
            Cell< Cell < MSI::Dof_Type > > tDofTypes = {{ MSI::Dof_Type::TEMP }};
            Field_Interpolator_Manager tFIManager( tDofTypes, &tSet );
            tFIManager.mFI = tFIs;

            // set field interpolators
            tCMMasterDiffLinIso->set_field_interpolator_manager( &tFIManager );
            tCMMasterDiffLinIso->set_geometry_interpolator( &tGI );

            //tCMMasterDiffLinIso->get_global_dof_type_list();

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
            REQUIRE( tCheckdStress );

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
            REQUIRE( tCheckdStrain );

            // check constitutive matrix----------------------------------------------------
            //------------------------------------------------------------------------------
            // evaluate the constitutive model constitutive matrix
            Matrix< DDRMat > tConst = tCMMasterDiffLinIso->constitutive();
            //print( tConst, "tConst");

            // evaluate the constitutive model constitutive matrix derivative
            Matrix< DDRMat > tdConstdDOF = tCMMasterDiffLinIso->dConstdDOF( { MSI::Dof_Type::TEMP } );
            //print( tdConstdDOF, "tdConstdDOF" );

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
            Matrix< DDRMat > tTestTraction = tCMMasterDiffLinIso->testTraction( tNormal );
            //print( tTestTraction, "tTestTraction");

            // evaluate the constitutive model test traction derivative
            Matrix< DDRMat > tdTestTractiondDOF = tCMMasterDiffLinIso->dTestTractiondDOF( { MSI::Dof_Type::TEMP }, tNormal,{{0.0}} );
            //print( tdTestTractiondDOF, "tdTestTractiondDOF" );

            // check test strain------------------------------------------------------------
            //------------------------------------------------------------------------------
            // evaluate the constitutive model test strain
            Matrix< DDRMat > tTestStrain = tCMMasterDiffLinIso->testStrain();
            //print( tTestStrain, "tTestStrain");

            // clean up
            //------------------------------------------------------------------------------
            // delete the field interpolator pointers
//            for( Field_Interpolator* tFI : tFIs )
//            {
//                delete tFI;
//            }
            tFIs.clear();

            //delete tSet;

        }/* TEST_CASE */


    }/* namespace fem */
}/* namespace moris */
