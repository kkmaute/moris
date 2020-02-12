/*
 * cl_FEM_Parameters.cpp
 *
 *  Created on: Feb 6, 2020
 *      Author: noel
 */

#include "cl_FEM_Parameters.hpp" //FEM/INT/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src

#include "assert.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_unique.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
//    // creates a FEM parameter list with default inputs
//    moris::Cell< moris::Cell< ParameterList > > create_fem_parameter_list()
//    {
//        uint tNumProperties              = 0;
//        uint tNumConstitutiveModels      = 0;
//        uint tNumStabilizationParameters = 0;
//        uint tNumIWGs                    = 0;
//        uint tNumIQIs                    = 0;
//
//        moris::Cell< moris::Cell< ParameterList > > tFEMParameterList;
//        tFEMParameterList( 0 ).resize( tNumProperties );              // Properties
//        tFEMParameterList( 1 ).resize( tNumConstitutiveModels );      // Constitutive models
//        tFEMParameterList( 2 ).resize( tNumStabilizationParameters ); // Stabilization parameters
//        tFEMParameterList( 3 ).resize( tNumIWGs );                    // IWGs
//        tFEMParameterList( 4 ).resize( tNumIQIs );                    // IQIs
//
//        // create property parameter list
//        for( uint iProp = 0; iProp < tNumProperties; iProp++ )
//        {
//            tFEMParameterList( 0 )( iProp ) = create_property_parameter_list();
//        }
//
//        // create constitutive model parameter list
//        for( uint iCM = 0; iCM < tNumConstitutiveModels; iCM++ )
//        {
//            tFEMParameterList( 1 )( iCM ) = create_constitutive_model_parameter_list();
//        }
//
//        // create stabilization parameter parameter list
//        for( uint iSP = 0; iSP < tNumStabilizationParameters; iSP++ )
//        {
//            tFEMParameterList( 2 )( iSP ) = create_stabilization_parameter_parameter_list();
//        }
//
//        // create IWG parameter list
//        for( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
//        {
//            tFEMParameterList( 3 )( iIWG ) = create_IWG_parameter_list();
//        }
//
//        // create IQI parameter list
//        for( uint iIQI = 0; iIQI < tNumIWGs; iIQI++ )
//        {
//            tFEMParameterList( 4 )( iIQI ) = create_IQI_parameter_list();
//        }
//
//        return tFEMParameterList;
//    }

//------------------------------------------------------------------------------
        ParameterList create_property_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "property_name",            std::string( "undefined" ) );
            tParameterList.insert( "function_parameters",      std::string( "" ) );
            tParameterList.insert( "value_function",           std::string( "" ) );
            tParameterList.insert( "dof_derivative_functions", std::string( "" ) );
            tParameterList.insert( "dv_derivative_functions",  std::string( "" ) );
            tParameterList.insert( "dof_dependencies",         std::string( "" ) );
            tParameterList.insert( "dv_dependencies",          std::string( "" ) );

            return tParameterList;
        }

//------------------------------------------------------------------------------
        // creates a constitutive model parameter list with default inputs
        ParameterList create_constitutive_model_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "constitutive_name", std::string( "undefined" ) );
            tParameterList.insert( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::UNDEFINED ) );
            tParameterList.insert( "dof_dependencies",  std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "dv_dependencies",   std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "properties",        std::string( "" ) );

            return tParameterList;
        }

//------------------------------------------------------------------------------
        // creates a stabilization parameter parameter list with default inputs
        ParameterList create_stabilization_parameter_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "stabilization_name",         std::string( "undefined" ) );
            tParameterList.insert( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::UNDEFINED ) );
            tParameterList.insert( "function_parameters",        std::string( "" ) );
            tParameterList.insert( "master_dof_dependencies",    std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "slave_dof_dependencies",     std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "master_dv_dependencies",     std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "slave_dv_dependencies",      std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "master_properties",          std::string( "" ) );
            tParameterList.insert( "slave_properties",           std::string( "" ) );
            tParameterList.insert( "master_constitutive_models", std::string( "" ) );
            tParameterList.insert( "slave_constitutive_models",  std::string( "" ) );

            return tParameterList;
        }

//------------------------------------------------------------------------------
        // creates an IWG parameter list with default inputs
        ParameterList create_IWG_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "IWG_name",                   std::string( "undefined" ) );
            tParameterList.insert( "IWG_type",                   static_cast< uint >( fem::IWG_Type::UNDEFINED ) );
            tParameterList.insert( "dof_residual",               std::string( "" ) );
            tParameterList.insert( "master_dof_dependencies",    std::string( "" ) );
            tParameterList.insert( "slave_dof_dependencies",     std::string( "" ) );
            tParameterList.insert( "master_dv_dependencies",     std::string( "" ) );
            tParameterList.insert( "slave_dv_dependencies",      std::string( "" ) );
            tParameterList.insert( "master_properties",          std::string( "" ) );
            tParameterList.insert( "slave_properties",           std::string( "" ) );
            tParameterList.insert( "master_constitutive_models", std::string( "" ) );
            tParameterList.insert( "slave_constitutive_models",  std::string( "" ) );
            tParameterList.insert( "stabilization_parameters",   std::string( "" ) );

            return tParameterList;
        }

//------------------------------------------------------------------------------
        // creates an IQI parameter list with default inputs
        ParameterList create_IQI_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "IQI_name",                   std::string( "undefined" ) );
            tParameterList.insert( "IQI_type",                   static_cast< uint >( fem::IQI_Type::UNDEFINED ) );
            tParameterList.insert( "master_dof_dependencies",    std::string( "" ) );
            tParameterList.insert( "slave_dof_dependencies",     std::string( "" ) );
            tParameterList.insert( "master_dv_dependencies",     std::string( "" ) );
            tParameterList.insert( "slave_dv_dependencies",      std::string( "" ) );
            tParameterList.insert( "master_properties",          std::string( "" ) );
            tParameterList.insert( "slave_properties",           std::string( "" ) );
            tParameterList.insert( "master_constitutive_models", std::string( "" ) );
            tParameterList.insert( "slave_constitutive_models",  std::string( "" ) );
            tParameterList.insert( "stabilization_parameters",   std::string( "" ) );

            return tParameterList;
        }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */


