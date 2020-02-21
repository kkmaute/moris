/*
 * cl_PRM_FEM_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: noel
 */

#ifndef PROJECTS_PRM_SRC_CL_PRM_FEM_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_CL_PRM_FEM_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_XML_Parser.hpp"
#include "cl_Param_List.hpp"
#include "cl_FEM_Enums.hpp" //FEM/INT/src

#include "cl_FEM_Enums.hpp" //FEM/INT/src

#include "cl_VIS_Output_Enums.hpp"


namespace moris
{
    namespace prm
    {

//------------------------------------------------------------------------------
        /*
         * creates a property parameter list with default inputs
         * @param [ out ] ParameterList a property parameter list
         */
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
        /*
         * creates a constitutive model parameter list with default inputs
         * @param [ out ] ParameterList a CM parameter list
         */
        ParameterList create_constitutive_model_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "constitutive_name", std::string( "undefined" ) );
            tParameterList.insert( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::UNDEFINED ) );
            tParameterList.insert( "dof_dependencies",  std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "dv_dependencies",   std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "properties",        std::string( "" ) );
            tParameterList.insert( "model_type",        static_cast< uint >( fem::Model_Type::UNDEFINED ) );

            return tParameterList;
        }

//------------------------------------------------------------------------------
        /*
         * creates a stabilization parameter parameter list with default inputs
         * @param [ out ] ParameterList a SP parameter list
         */
        ParameterList create_stabilization_parameter_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "stabilization_name",         std::string( "undefined" ) );
            tParameterList.insert( "stabilization_type",         static_cast< uint >( fem::Stabilization_Type::UNDEFINED ) );
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
        /*
         * creates an IWG parameter list with default inputs
         * @param [ out ] ParameterList a IWG parameter list
         */
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

            tParameterList.insert( "mesh_set_names",             std::string( "" ) );

            return tParameterList;
        }

//------------------------------------------------------------------------------
        /*
         * creates an IQI parameter list with default inputs
         * @param [ out ] ParameterList a IQI parameter list
         */
        ParameterList create_IQI_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "IQI_name",                   std::string( "undefined" ) );
            tParameterList.insert( "IQI_type",                   static_cast< uint >( fem::IQI_Type::UNDEFINED ) );
            tParameterList.insert( "IQI_output_type",            static_cast< uint >( vis::Output_Type::UNDEFINED ) );
            tParameterList.insert( "master_dof_dependencies",    std::string( "" ) );
            tParameterList.insert( "slave_dof_dependencies",     std::string( "" ) );
            tParameterList.insert( "master_dv_dependencies",     std::string( "" ) );
            tParameterList.insert( "slave_dv_dependencies",      std::string( "" ) );
            tParameterList.insert( "master_properties",          std::string( "" ) );
            tParameterList.insert( "slave_properties",           std::string( "" ) );
            tParameterList.insert( "master_constitutive_models", std::string( "" ) );
            tParameterList.insert( "slave_constitutive_models",  std::string( "" ) );
            tParameterList.insert( "stabilization_parameters",   std::string( "" ) );
            tParameterList.insert( "vectorial_field_index",      -1 );

            tParameterList.insert( "mesh_set_names",             std::string( "" ) );

            return tParameterList;
        }

//------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_CL_PRM_FEM_PARAMETERS_HPP_ */
