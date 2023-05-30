/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Input.cpp
 *
 */

#include "catch.hpp"
#include "paths.hpp"
#include "cl_FEM_Model.hpp" //FEM/INT/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
//#include "fn_Parsing_Tools.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
        // This test case tests the FEM inputs
//        TEST_CASE( "FEM Input", "[moris],[fem],[FEM_Input]" )
//        {
//            // create a cell of cell of parameter list for fem
//            moris::Cell< moris::Cell< ParameterList > > tParameterList( 5 );
//
//            //------------------------------------------------------------------------------
//            // fill the property part of the parameter list
//            uint tNumProperties = 2;
//            tParameterList( 0 ).resize( tNumProperties );
//
//            // create parameter list for property 1
//            tParameterList( 0 )( 0 ) = prm::create_property_parameter_list();
//            tParameterList( 0 )( 0 ).set( "property_name",            "Property1" );
//            tParameterList( 0 )( 0 ).set( "dof_dependencies",         "UX,UY;TEMP" );
//            tParameterList( 0 )( 0 ).set( "function_parameters",      "1.0, 2.0, 3.0; 4.0, 5.0, 6.0; 7.0, 8.0, 9.0 / 10.0, 11.0, 12.0; 13.0, 14.0, 15.0; 16.0, 17.0, 18.0" );
//            tParameterList( 0 )( 0 ).set( "value_function",           "Func1" );
//            tParameterList( 0 )( 0 ).set( "dof_derivative_functions", "DofDerFunc1,DofDerFunc2" );
//            tParameterList( 0 )( 0 ).set( "dv_derivative_functions",  "DvDerFunc1" );
//
//            // create parameter list for property 2
//            tParameterList( 0 )( 1 ) = prm::create_property_parameter_list();
//            tParameterList( 0 )( 1 ).set( "property_name",    "Property2" );
//            tParameterList( 0 )( 1 ).set( "dof_dependencies", "UX,UY;TEMP" );
//
//            //------------------------------------------------------------------------------
//            // fill the constitutive model part of the parameter list
//            uint tNumCMs = 2;
//            tParameterList( 1 ).resize( tNumCMs );
//
//            // create parameter list for constitutive model 1
//            tParameterList( 1 )( 0 ) = prm::create_constitutive_model_parameter_list();
//            tParameterList( 1 )( 0 ).set( "constitutive_name", "CM1" );
//            tParameterList( 1 )( 0 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
//            tParameterList( 1 )( 0 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
//            tParameterList( 1 )( 0 ).set( "properties",        "Property1,Conductivity" );
//
//            // create parameter list for constitutive model 2
//            tParameterList( 1 )( 1 ) = prm::create_constitutive_model_parameter_list();
//            tParameterList( 1 )( 1 ).set( "constitutive_name", "CM2" );
//            tParameterList( 1 )( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::DIFF_LIN_ISO ) );
//            tParameterList( 1 )( 1 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
//            tParameterList( 1 )( 1 ).set( "properties",        "Property2,Conductivity" );
//
//            //------------------------------------------------------------------------------
//            // fill the stabilization parameter part of the parameter list
//            uint tNumSPs = 2;
//            tParameterList( 2 ).resize( tNumSPs );
//
//            // create parameter list for stabilization parameter 1
//            tParameterList( 2 )( 0 ) = prm::create_stabilization_parameter_parameter_list();
//            tParameterList( 2 )( 0 ).set( "stabilization_name",  "SP1" );
//            tParameterList( 2 )( 0 ).set( "stabilization_type",  static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
//            tParameterList( 2 )( 0 ).set( "function_parameters", "1.0, 2.0, 3.0; 4.0, 5.0, 6.0; 7.0, 8.0, 9.0 / 10.0, 11.0, 12.0; 13.0, 14.0, 15.0; 16.0, 17.0, 18.0" );
//            tParameterList( 2 )( 0 ).set( "leader_properties",   "Property1,Material" );
//
//            // create parameter list for stabilization parameter 2
//            tParameterList( 2 )( 1 ) = prm::create_stabilization_parameter_parameter_list();
//            tParameterList( 2 )( 1 ).set( "stabilization_name", "SP2" );
//            tParameterList( 2 )( 1 ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::NITSCHE_INTERFACE ) );
//            tParameterList( 2 )( 1 ).set( "leader_properties",  "Property1,Material" );
//            tParameterList( 2 )( 1 ).set( "follower_properties",   "Property2,Material" );
//
//            //------------------------------------------------------------------------------
//            // fill the IWG part of the parameter list
//            uint tNumIWGs = 3;
//            tParameterList( 3 ).resize( tNumIWGs );
//
//            // create parameter list for IWG 1
//            tParameterList( 3 )( 0 ) = prm::create_IWG_parameter_list();
//            tParameterList( 3 )( 0 ).set( "IWG_name",                   "IWG1" );
//            tParameterList( 3 )( 0 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_BULK ) );
//            tParameterList( 3 )( 0 ).set( "dof_residual",               "TEMP" );
//            tParameterList( 3 )( 0 ).set( "leader_dof_dependencies",    "TEMP" );
//            tParameterList( 3 )( 0 ).set( "leader_properties",          "Property1,Load" );
//            tParameterList( 3 )( 0 ).set( "leader_constitutive_models", "CM1,DiffLinIso" );
//            tParameterList( 3 )( 0 ).set( "mesh_set_names",             "MeshSet1,MeshSet2" );
//
//            // create parameter list for IWG 2
//            tParameterList( 3 )( 1 ) = prm::create_IWG_parameter_list();
//            tParameterList( 3 )( 1 ).set( "IWG_name",                   "IWG2" );
//            tParameterList( 3 )( 1 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::SPATIALDIFF_INTERFACE ) );
//            tParameterList( 3 )( 1 ).set( "dof_residual",               "TEMP" );
//            tParameterList( 3 )( 1 ).set( "leader_dof_dependencies",    "TEMP" );
//            tParameterList( 3 )( 1 ).set( "follower_dof_dependencies",     "TEMP" );
//            tParameterList( 3 )( 1 ).set( "leader_constitutive_models", "CM1,DiffLinIso" );
//            tParameterList( 3 )( 1 ).set( "follower_constitutive_models",  "CM2,DiffLinIso" );
//            tParameterList( 3 )( 1 ).set( "stabilization_parameters",   "SP2,NitscheInterface" );
//            tParameterList( 3 )( 1 ).set( "mesh_set_names",             "MeshSet2" );
//
//
//            // create parameter list for IWG 3
//            tParameterList( 3 )( 2 ) = prm::create_IWG_parameter_list();
//            tParameterList( 3 )( 2 ).set( "IWG_name",                   "IWG3" );
//            tParameterList( 3 )( 2 ).set( "IWG_type",                   static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
//            tParameterList( 3 )( 2 ).set( "dof_residual",               "UX,UY,UZ" );
//            tParameterList( 3 )( 2 ).set( "leader_dof_dependencies",    "UX,UY,UZ" );
//            tParameterList( 3 )( 2 ).set( "leader_constitutive_models", "CM1,ElastLinIso" );
//            tParameterList( 3 )( 2 ).set( "leader_properties",          "Property2,Dirichlet;Property1,Select" );
//            tParameterList( 3 )( 2 ).set( "mesh_set_names",             "MeshSet2,MeshSet3,MeshSet4" );
//
//            //------------------------------------------------------------------------------
//            // fill the IQI part of the parameter list
//            uint tNumIQIs = 2;
//            tParameterList( 4 ).resize( tNumIQIs );
//
//            // create parameter list for IQI 1
//            tParameterList( 4 )( 0 ) = prm::create_IQI_parameter_list();
//            tParameterList( 4 )( 0 ).set( "IQI_name",                   "IQI1" );
//            tParameterList( 4 )( 0 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
//            tParameterList( 4 )( 0 ).set( "leader_dof_dependencies",    "TEMP" );
//            tParameterList( 4 )( 0 ).set( "leader_constitutive_models", "CM1,Elast" );
//            tParameterList( 4 )( 0 ).set( "mesh_set_names",             "MeshSet5,MeshSet2" );
//
//            // create parameter list for IQI 2
//            tParameterList( 4 )( 1 ) = prm::create_IQI_parameter_list();
//            tParameterList( 4 )( 1 ).set( "IQI_name",                   "IQI2" );
//            tParameterList( 4 )( 1 ).set( "IQI_type",                   static_cast< uint >( fem::IQI_Type::STRAIN_ENERGY ) );
//            tParameterList( 4 )( 1 ).set( "leader_dof_dependencies",    "TEMP" );
//            tParameterList( 4 )( 1 ).set( "leader_constitutive_models", "CM2,Elast" );
//            tParameterList( 4 )( 1 ).set( "mesh_set_names",             "MeshSet5" );
//
//            //------------------------------------------------------------------------------
//            // open input file
//            std::string tMeshFilePath = moris::get_base_moris_dir();
//            tMeshFilePath = tMeshFilePath + "projects/FEM/INT/test/data/FEM_input_test.so";
//            std::shared_ptr< Library_IO > tLibrary = std::make_shared< Library_IO >( tMeshFilePath );
//
//            // create a FEM model
//            FEM_Model tFEMModel;
//            tFEMModel.set_parameter_list( tParameterList );
//            tFEMModel.set_space_dim( 2 );
//            tFEMModel.initialize( tLibrary );
//
////            // parsing tool debug
////            std::string tString = " 1.0, 2.0, 3.0; 4.0, 5.0, 6.0; 7.0, 8.0, 9.0";
////            moris::Cell< Matrix< DDRMat > > tTest;
////            string_to_cell_mat_2( tString, tTest );
////            print( tTest, "tTest" );
////
////            std::string tString2 = " 1.0, 2.0, 3.0; 4.0, 5.0, 6.0; 7.0, 8.0, 9.0/ 10.0, 11.0, 12.0; 13.0, 14.0, 15.0; 16.0, 17.0, 18.0";
////            moris::Cell< Matrix< DDRMat > > tTest2;
////            string_to_cell_mat_2( tString2, tTest2 );
////            print( tTest2, "tTest2" );
//
////            std::string tString = "MeshSet5,MeshSet2";
////            moris::Cell< std::string > tStrings;
////            string_to_cell( tString, tStrings );
////            std::cout<<"1"<<tStrings( 0 )<<std::endl;
////            std::cout<<"2"<<tStrings( 1 )<<std::endl;
//
//        }/* END_TEST_CASE */
    }/* end_fem_namespace */
}/* end_moris_namespace */

