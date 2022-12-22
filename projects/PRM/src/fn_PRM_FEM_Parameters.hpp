/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_FEM_Parameters.hpp
 *
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_FEM_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_FEM_PARAMETERS_HPP_


#include "cl_Param_List.hpp"

#include "cl_FEM_Enums.hpp"
#include "cl_VIS_Output_Enums.hpp"

namespace moris
{
    namespace prm
    {
        //------------------------------------------------------------------------------

        /**
         * creates a property parameter list with default inputs
         * @param [ out ] ParameterList a property parameter list
         */
        inline ParameterList
        create_property_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "property_name", "undefined" );
            tParameterList.insert( "function_parameters", "" );
            tParameterList.insert( "value_function", "" );
            tParameterList.insert( "dof_derivative_functions", "" );
            tParameterList.insert( "dv_derivative_functions", "" );
            tParameterList.insert( "space_derivative_functions", "" );
            tParameterList.insert( "dof_dependencies", "" );
            tParameterList.insert( "dv_dependencies", "" );
            tParameterList.insert( "field_dependencies", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        /**
         * creates a field parameter list with default inputs
         * @param [ out ] ParameterList a property parameter list
         */
        inline ParameterList
        create_fem_field_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "field_name", "undefined" );
            tParameterList.insert( "field_entity_type", "" );
            tParameterList.insert( "field_type", "" );
            tParameterList.insert( "field_create_from_file", "" );
            tParameterList.insert( "IQI_Name", "" );
            tParameterList.insert( "field_output_to_file", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        /**
         * creates a material model parameter list with default inputs
         * @param [ out ] ParameterList a MM parameter list
         */
        inline ParameterList
        create_material_model_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "material_name", "undefined" );
            tParameterList.insert( "material_type", (uint)( fem::Material_Type::UNDEFINED ) );
            tParameterList.insert( "dof_dependencies", std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "properties", "" );
            tParameterList.insert( "phase_name", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        /**
         * creates a constitutive model parameter list with default inputs
         * @param [ out ] ParameterList a CM parameter list
         */
        inline ParameterList
        create_constitutive_model_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "constitutive_name", "undefined" );
            tParameterList.insert( "constitutive_type", (uint)( fem::Constitutive_Type::UNDEFINED ) );
            tParameterList.insert( "function_parameters", "" );
            tParameterList.insert( "dof_dependencies", std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "dv_dependencies", std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "properties", "" );
            tParameterList.insert( "material_model", "" );
            tParameterList.insert( "model_type", (uint)( fem::Model_Type::UNDEFINED ) );
            tParameterList.insert( "phase_name", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        /**
         * creates a stabilization parameter parameter list with default inputs
         * @param [ out ] ParameterList a SP parameter list
         */
        inline ParameterList
        create_stabilization_parameter_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "stabilization_name", "undefined" );
            tParameterList.insert( "stabilization_type", (uint)( fem::Stabilization_Type::UNDEFINED ) );
            tParameterList.insert( "function_parameters", "" );
            tParameterList.insert( "master_dof_dependencies", std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "slave_dof_dependencies", std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "master_dv_dependencies", std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "slave_dv_dependencies", std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "master_properties", "" );
            tParameterList.insert( "slave_properties", "" );
            tParameterList.insert( "master_constitutive_models", "" );
            tParameterList.insert( "slave_constitutive_models", "" );
            tParameterList.insert( "cluster_measures", std::pair< std::string, std::string >( "", "" ) );
            tParameterList.insert( "stabilization_dof_type", "" );

            tParameterList.insert( "master_phase_name", "" );
            tParameterList.insert( "slave_phase_name", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        /**
         * creates an IWG parameter list with default inputs
         * @param [ out ] ParameterList a IWG parameter list
         */
        inline ParameterList
        create_IWG_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "IWG_name", "undefined" );
            tParameterList.insert( "IWG_bulk_type", (uint)( fem::Element_Type::BULK ) );
            tParameterList.insert( "IWG_type", (uint)( fem::IWG_Type::UNDEFINED ) );
            tParameterList.insert( "dof_residual", "" );

            tParameterList.insert( "master_phase_name", "" );
            tParameterList.insert( "slave_phase_name", "" );
            tParameterList.insert( "master_dof_dependencies", "" );
            tParameterList.insert( "slave_dof_dependencies", "" );
            tParameterList.insert( "master_dv_dependencies", "" );
            tParameterList.insert( "slave_dv_dependencies", "" );
            tParameterList.insert( "master_field_types", "" );
            tParameterList.insert( "slave_field_types", "" );

            tParameterList.insert( "master_properties", "" );
            tParameterList.insert( "slave_properties", "" );
            tParameterList.insert( "master_material_model", "" );
            tParameterList.insert( "slave_material_model", "" );
            tParameterList.insert( "master_constitutive_models", "" );
            tParameterList.insert( "slave_constitutive_models", "" );
            tParameterList.insert( "stabilization_parameters", "" );

            tParameterList.insert( "master_physics", "" );
            tParameterList.insert( "slave_physics", "" );
            tParameterList.insert( "stabilizations", "" );

            tParameterList.insert( "ghost_order", MORIS_UINT_MAX );

            tParameterList.insert( "mesh_set_names", "" );
            tParameterList.insert( "master_phase_name", "" );
            tParameterList.insert( "slave_phase_name", "" );
            tParameterList.insert( "side_ordinals", "" );
            tParameterList.insert( "neighbor_phases", "" );

            tParameterList.insert( "time_continuity", false );
            tParameterList.insert( "time_boundary", false );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        /**
         * creates an IQI parameter list with default inputs
         * @param [ out ] ParameterList a IQI parameter list
         */
        inline ParameterList
        create_IQI_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "IQI_name", "undefined" );
            tParameterList.insert( "IQI_type", (uint)( fem::IQI_Type::UNDEFINED ) );
            tParameterList.insert( "IQI_bulk_type", (uint)( fem::Element_Type::BULK ) );
            tParameterList.insert( "dof_quantity", "" );
            tParameterList.insert( "vectorial_field_index", -1 );

            tParameterList.insert( "master_phase_name", "" );
            tParameterList.insert( "slave_phase_name", "" );
            tParameterList.insert( "master_dof_dependencies", "" );
            tParameterList.insert( "slave_dof_dependencies", "" );
            tParameterList.insert( "master_dv_dependencies", "" );
            tParameterList.insert( "slave_dv_dependencies", "" );
            tParameterList.insert( "master_field_types", "" );
            tParameterList.insert( "slave_field_types", "" );

            tParameterList.insert( "function_parameters", "" );

            tParameterList.insert( "master_properties", "" );
            tParameterList.insert( "slave_properties", "" );
            tParameterList.insert( "master_constitutive_models", "" );
            tParameterList.insert( "slave_constitutive_models", "" );
            tParameterList.insert( "stabilization_parameters", "" );

            tParameterList.insert( "master_physics", "" );
            tParameterList.insert( "slave_physics", "" );
            tParameterList.insert( "stabilizations", "" );

            tParameterList.insert( "master_phase_name", "" );
            tParameterList.insert( "slave_phase_name", "" );
            tParameterList.insert( "mesh_set_names", "" );
            tParameterList.insert( "side_ordinals", "" );
            tParameterList.insert( "neighbor_phases", "" );

            tParameterList.insert( "time_continuity", false );
            tParameterList.insert( "time_boundary", false );

            tParameterList.insert( "normalization", "none" );    // options: time, design, vector of reference values

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        /**
         * creates phase parameter list with default inputs
         * @param [ out ] ParameterList a phase parameter list
         */
        inline ParameterList
        create_phase_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "phase_name", "undefined" );
            tParameterList.insert( "phase_indices", "" );
            tParameterList.insert( "dof_dependencies", "" );
            tParameterList.insert( "dv_dependencies", "" );
            tParameterList.insert( "constitutive_models", "" );
            tParameterList.insert( "stabilization_parameters", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

        /**
         * creates computation fem parameter list with default inputs
         * @param [ out ] ParameterList a computation parameter list
         */
        inline ParameterList
        create_computation_parameter_list()
        {
            ParameterList tParameterList;

            // bool true for printing physics
            tParameterList.insert( "print_physics_model", false );

            // bool true for analytical forward analysis, false for finite difference
            // decide if dRdu (jacobian) and dQIdu are computed by A/FD
            tParameterList.insert( "is_analytical_forward", true );

            // enum for finite difference scheme for forward analysis
            tParameterList.insert( "finite_difference_scheme_forward", (uint)( fem::FDScheme_Type::POINT_1_FORWARD ) );

            // real for relative perturbation size for finite difference for forward analysis
            tParameterList.insert( "finite_difference_perturbation_size_forward", 1e-6 );

            // bool true for analytical sensitivity analysis, false for finite difference
            // decide if dRdp and dQIdp are computed by A/FD
            tParameterList.insert( "is_analytical_sensitivity", false );

            // enum for finite difference scheme for sensitivity analysis
            tParameterList.insert( "finite_difference_scheme", (uint)( fem::FDScheme_Type::POINT_1_FORWARD ) );

            // real for relative perturbation size for finite difference
            tParameterList.insert( "finite_difference_perturbation_size", 1e-6 );

            // enum for finite difference perturbation strategy (relative, absolute)
            tParameterList.insert( "finite_difference_perturbation_strategy", (uint)( fem::Perturbation_Type::RELATIVE ) );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif /* PROJECTS_PRM_SRC_FN_PRM_FEM_PARAMETERS_HPP_ */
