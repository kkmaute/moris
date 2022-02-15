/*
 * fn_PRM_FEM_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: noel
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_FEM_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_FEM_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_XML_Parser.hpp"
#include "cl_Param_List.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
//FEM/VIS/src
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
        inline
        ParameterList create_property_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "property_name",              std::string( "undefined" ) );
            tParameterList.insert( "function_parameters",        std::string( "" ) );
            tParameterList.insert( "value_function",             std::string( "" ) );
            tParameterList.insert( "dof_derivative_functions",   std::string( "" ) );
            tParameterList.insert( "dv_derivative_functions",    std::string( "" ) );
            tParameterList.insert( "space_derivative_functions", std::string( "" ) );
            tParameterList.insert( "dof_dependencies",           std::string( "" ) );
            tParameterList.insert( "dv_dependencies",            std::string( "" ) );
            tParameterList.insert( "field_dependencies",         std::string( "" ) );

            return tParameterList;
        }

        //------------------------------------------------------------------------------
        /*
         * creates a field parameter list with default inputs
         * @param [ out ] ParameterList a property parameter list
         */
        inline
        ParameterList create_fem_field_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "field_name",             std::string( "undefined" ) );
            tParameterList.insert( "field_entity_type",      std::string( "" ) );
            tParameterList.insert( "field_type",             std::string( "" ) );
            tParameterList.insert( "field_create_from_file", std::string( "" ) );
            tParameterList.insert( "IQI_Name",               std::string( "" ) );
            tParameterList.insert( "field_output_to_file",   std::string( "" ) );

            return tParameterList;
        }

        //------------------------------------------------------------------------------
        /*
        * creates a material model parameter list with default inputs
        * @param [ out ] ParameterList a MM parameter list
        */
        inline
        ParameterList create_material_model_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "material_name", std::string( "undefined" ) );
            tParameterList.insert( "material_type", static_cast< uint >( fem::Material_Type::UNDEFINED ) );
            tParameterList.insert( "dof_dependencies",  std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "properties",        std::string( "" ) );
            tParameterList.insert( "phase_name",        std::string( "" ) );

            return tParameterList;
        }

        //------------------------------------------------------------------------------
        /*
         * creates a constitutive model parameter list with default inputs
         * @param [ out ] ParameterList a CM parameter list
         */
        inline
        ParameterList create_constitutive_model_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "constitutive_name", std::string( "undefined" ) );
            tParameterList.insert( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::UNDEFINED ) );
            tParameterList.insert( "dof_dependencies",  std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "dv_dependencies",   std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "properties",        std::string( "" ) );
            tParameterList.insert( "material_model",    std::string( "" ) );
            tParameterList.insert( "model_type",        static_cast< uint >( fem::Model_Type::UNDEFINED ) );
            tParameterList.insert( "phase_name",        std::string( "" ) );

            return tParameterList;
        }

        //------------------------------------------------------------------------------
        /*
         * creates a stabilization parameter parameter list with default inputs
         * @param [ out ] ParameterList a SP parameter list
         */
        inline
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
            tParameterList.insert( "cluster_measures",           std::pair< std::string, std::string >( "", "" ) );

            tParameterList.insert( "master_phase_name",          std::string( "" ) );
            tParameterList.insert( "slave_phase_name",           std::string( "" ) );

            return tParameterList;
        }

        //------------------------------------------------------------------------------
        /*
         * creates an IWG parameter list with default inputs
         * @param [ out ] ParameterList a IWG parameter list
         */
        inline
        ParameterList create_IWG_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "IWG_name",                   std::string( "undefined" ) );
            tParameterList.insert( "IWG_bulk_type",              static_cast< uint >( fem::Element_Type::BULK ) );
            tParameterList.insert( "IWG_type",                   static_cast< uint >( fem::IWG_Type::UNDEFINED ) );
            tParameterList.insert( "dof_residual",               std::string( "" ) );

            tParameterList.insert( "master_phase_name",          "" );
            tParameterList.insert( "slave_phase_name",           "" );
            tParameterList.insert( "master_dof_dependencies",    std::string( "" ) );
            tParameterList.insert( "slave_dof_dependencies",     std::string( "" ) );
            tParameterList.insert( "master_dv_dependencies",     std::string( "" ) );
            tParameterList.insert( "slave_dv_dependencies",      std::string( "" ) );
            tParameterList.insert( "master_field_types",         std::string( "" ) );
            tParameterList.insert( "slave_field_types",          std::string( "" ) );

            tParameterList.insert( "master_properties",          std::string( "" ) );
            tParameterList.insert( "slave_properties",           std::string( "" ) );
            tParameterList.insert( "master_material_model",      std::string( "" ) );
            tParameterList.insert( "slave_material_model",       std::string( "" ) );
            tParameterList.insert( "master_constitutive_models", std::string( "" ) );
            tParameterList.insert( "slave_constitutive_models",  std::string( "" ) );
            tParameterList.insert( "stabilization_parameters",   std::string( "" ) );

            tParameterList.insert( "master_physics",             std::string( "" ) );
            tParameterList.insert( "slave_physics",              std::string( "" ) );
            tParameterList.insert( "stabilizations",             std::string( "" ) );

            tParameterList.insert( "ghost_order",                MORIS_UINT_MAX );

            tParameterList.insert( "mesh_set_names",             std::string( "" ) );
            tParameterList.insert( "master_phase_name",          std::string( "" ) );
            tParameterList.insert( "slave_phase_name",           std::string( "" ) );
            tParameterList.insert( "side_ordinals",              std::string( "" ) );
            tParameterList.insert( "neighbor_phases",            std::string( "" ) );

            tParameterList.insert( "time_continuity",            false );
            tParameterList.insert( "time_boundary",              false );

            return tParameterList;
        }

        //------------------------------------------------------------------------------
        /*
         * creates an IQI parameter list with default inputs
         * @param [ out ] ParameterList a IQI parameter list
         */
        inline
        ParameterList create_IQI_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "IQI_name",                   std::string( "undefined" ) );
            tParameterList.insert( "IQI_type",                   static_cast< uint >( fem::IQI_Type::UNDEFINED ) );
            tParameterList.insert( "IQI_bulk_type",              static_cast< uint >( fem::Element_Type::BULK ) );
            tParameterList.insert( "dof_quantity",               std::string( "" ) );
            tParameterList.insert( "vectorial_field_index",      -1 );

            tParameterList.insert( "master_phase_name",          "" );
            tParameterList.insert( "slave_phase_name",           "" );
            tParameterList.insert( "master_dof_dependencies",    std::string( "" ) );
            tParameterList.insert( "slave_dof_dependencies",     std::string( "" ) );
            tParameterList.insert( "master_dv_dependencies",     std::string( "" ) );
            tParameterList.insert( "slave_dv_dependencies",      std::string( "" ) );
            tParameterList.insert( "master_field_types",         std::string( "" ) );
            tParameterList.insert( "slave_field_types",          std::string( "" ) );

            tParameterList.insert( "function_parameters",        std::string( "" ) );

            tParameterList.insert( "master_properties",          std::string( "" ) );
            tParameterList.insert( "slave_properties",           std::string( "" ) );
            tParameterList.insert( "master_constitutive_models", std::string( "" ) );
            tParameterList.insert( "slave_constitutive_models",  std::string( "" ) );
            tParameterList.insert( "stabilization_parameters",   std::string( "" ) );

            tParameterList.insert( "master_physics",             std::string( "" ) );
            tParameterList.insert( "slave_physics",              std::string( "" ) );
            tParameterList.insert( "stabilizations",             std::string( "" ) );

            tParameterList.insert( "master_phase_name",          std::string( "" ) );
            tParameterList.insert( "slave_phase_name",           std::string( "" ) );
            tParameterList.insert( "mesh_set_names",             std::string( "" ) );
            tParameterList.insert( "side_ordinals",              std::string( "" ) );
            tParameterList.insert( "neighbor_phases",            std::string( "" ) );

            tParameterList.insert( "time_continuity",            false );
            tParameterList.insert( "time_boundary",              false );

            tParameterList.insert( "normalization",              "none" ); // options: time, design, vector of reference values

            return tParameterList;
        }

        //------------------------------------------------------------------------------
        /*
         * creates phase parameter list with default inputs
         * @param [ out ] ParameterList a phase parameter list
         */
        inline
        ParameterList create_phase_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "phase_name",               "undefined" );
            tParameterList.insert( "phase_indices",            "" );
            tParameterList.insert( "dof_dependencies",         "" );
            tParameterList.insert( "dv_dependencies",          "" );
            tParameterList.insert( "constitutive_models",      "" );
            tParameterList.insert( "stabilization_parameters", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------
        /*
         * creates computation fem parameter list with default inputs
         * @param [ out ] ParameterList a computation parameter list
         */
        inline
        ParameterList create_computation_parameter_list()
        {
            ParameterList tParameterList;

            // bool true for printing physics
            tParameterList.insert( "print_physics_model", false );

            // bool true for analytical forward analysis, false for finite difference
            // decide if dRdu (jacobian) and dQIdu are computed by A/FD
            tParameterList.insert( "is_analytical_forward", true );

            // enum for finite difference scheme for forward analysis
            tParameterList.insert( "finite_difference_scheme_forward",
                    static_cast< uint >( fem::FDScheme_Type::POINT_1_FORWARD ) );

            // real for relative perturbation size for finite difference for forward analysis
            tParameterList.insert( "finite_difference_perturbation_size_forward", 1e-6 );

            // bool true for analytical sensitivity analysis, false for finite difference
            // decide if dRdp and dQIdp are computed by A/FD
            tParameterList.insert( "is_analytical_sensitivity", false );

            // enum for finite difference scheme for sensitivity analysis
            tParameterList.insert( "finite_difference_scheme",
                    static_cast< uint >( fem::FDScheme_Type::POINT_1_FORWARD ) );

            // real for relative perturbation size for finite difference
            tParameterList.insert( "finite_difference_perturbation_size", 1e-6 );

            // enum for finite difference perturbation strategy (relative, absolute)
            tParameterList.insert( "finite_difference_perturbation_strategy",
                    static_cast< uint >( fem::Perturbation_Type::RELATIVE ) );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_FN_PRM_FEM_PARAMETERS_HPP_ */
