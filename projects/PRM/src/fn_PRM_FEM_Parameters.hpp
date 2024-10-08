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

#include "cl_Parameter_List.hpp"

#include "cl_FEM_Enums.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_VIS_Output_Enums.hpp"

namespace moris::prm
{
    //------------------------------------------------------------------------------

    /**
     * creates a property parameter list with default inputs
     * @param [ out ] ParameterList a property parameter list
     */
    inline Parameter_List
    create_property_parameter_list()
    {
        Parameter_List tParameterList;

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
    inline Parameter_List
    create_fem_field_parameter_list()
    {
        Parameter_List tParameterList;

        tParameterList.insert( "field_name", "undefined" );
        tParameterList.insert( "field_entity_type", "" );
        tParameterList.insert( "field_type", "" );
        tParameterList.insert( "field_create_from_file", "" );
        tParameterList.insert( "field_file_time_index", 0 );
        tParameterList.insert( "field_file_field_index", 0 );
        tParameterList.insert( "field_file_update", false );
        tParameterList.insert( "IQI_Name", "" );
        tParameterList.insert( "field_output_to_file", "" );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * creates a material model parameter list with default inputs
     * @param [ out ] ParameterList a MM parameter list
     */
    inline Parameter_List
    create_material_model_parameter_list()
    {
        Parameter_List tParameterList;

        tParameterList.insert( "material_name", "undefined" );
        tParameterList.insert( "material_type", fem::Material_Type::UNDEFINED, fem::Material_Type::PERFECT_GAS, fem::Material_Type::VAN_DER_WAALS_FLUID );
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
    inline Parameter_List
    create_constitutive_model_parameter_list()
    {
        Parameter_List tParameterList;

        tParameterList.insert( "constitutive_name", "undefined" );
        tParameterList.insert( "constitutive_type", fem::Constitutive_Type::UNDEFINED, fem::Constitutive_Type::DIFF_LIN_ISO, fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
        tParameterList.insert( "function_parameters", "" );
        tParameterList.insert( "dof_dependencies", std::pair< std::string, std::string >( "", "" ) );
        tParameterList.insert( "dv_dependencies", std::pair< std::string, std::string >( "", "" ) );
        tParameterList.insert( "properties", "" );
        tParameterList.insert( "material_model", "" );
        tParameterList.insert( "model_type", fem::Model_Type::UNDEFINED, fem::Model_Type::PLANE_STRESS, fem::Model_Type::DEVIATORIC );
        tParameterList.insert( "phase_name", "" );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * creates a stabilization parameter parameter list with default inputs
     * @param [ out ] ParameterList a SP parameter list
     */
    inline Parameter_List
    create_stabilization_parameter_parameter_list()
    {
        Parameter_List tParameterList;

        tParameterList.insert( "stabilization_name", "undefined" );
        tParameterList.insert( "stabilization_type", fem::Stabilization_Type::UNDEFINED, fem::Stabilization_Type::DIRICHLET_NITSCHE, fem::Stabilization_Type::LAGRANGE_MULTIPLIER_L2 );
        tParameterList.insert( "function_parameters", "" );
        tParameterList.insert( "leader_dof_dependencies", std::pair< std::string, std::string >( "", "" ) );
        tParameterList.insert( "follower_dof_dependencies", std::pair< std::string, std::string >( "", "" ) );
        tParameterList.insert( "leader_dv_dependencies", std::pair< std::string, std::string >( "", "" ) );
        tParameterList.insert( "follower_dv_dependencies", std::pair< std::string, std::string >( "", "" ) );
        tParameterList.insert( "leader_properties", "" );
        tParameterList.insert( "follower_properties", "" );
        tParameterList.insert( "leader_constitutive_models", "" );
        tParameterList.insert( "follower_constitutive_models", "" );
        tParameterList.insert( "cluster_measures", std::pair< std::string, std::string >( "", "" ) );
        tParameterList.insert( "stabilization_dof_type", "" );

        tParameterList.insert( "leader_phase_name", "" );
        tParameterList.insert( "follower_phase_name", "" );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * creates an IWG parameter list with default inputs
     * @param [ out ] ParameterList a IWG parameter list
     */
    inline Parameter_List
    create_IWG_parameter_list()
    {
        Parameter_List tParameterList;

        tParameterList.insert( "IWG_name", "undefined" );
        tParameterList.insert( "IWG_bulk_type", fem::Element_Type::BULK, fem::Element_Type::BULK, fem::Element_Type::TIME_BOUNDARY );
        tParameterList.insert( "IWG_type", fem::IWG_Type::UNDEFINED, fem::IWG_Type::L2, fem::IWG_Type::USER_DEFINED );
        tParameterList.insert( "dof_residual", "" );

        tParameterList.insert( "leader_phase_name", "" );
        tParameterList.insert( "follower_phase_name", "" );
        tParameterList.insert( "leader_dof_dependencies", "" );
        tParameterList.insert( "follower_dof_dependencies", "" );
        tParameterList.insert( "leader_dv_dependencies", "" );
        tParameterList.insert( "follower_dv_dependencies", "" );
        tParameterList.insert( "leader_field_types", "" );
        tParameterList.insert( "follower_field_types", "" );

        tParameterList.insert( "function_parameters", "" );

        tParameterList.insert( "leader_properties", "" );
        tParameterList.insert( "follower_properties", "" );
        tParameterList.insert( "leader_material_model", "" );
        tParameterList.insert( "follower_material_model", "" );
        tParameterList.insert( "leader_constitutive_models", "" );
        tParameterList.insert( "follower_constitutive_models", "" );
        tParameterList.insert( "stabilization_parameters", "" );

        tParameterList.insert( "leader_physics", "" );
        tParameterList.insert( "follower_physics", "" );
        tParameterList.insert( "stabilizations", "" );

        tParameterList.insert( "ghost_order", MORIS_UINT_MAX );

        tParameterList.insert( "mesh_set_names", "" );
        tParameterList.insert( "leader_phase_name", "" );
        tParameterList.insert( "follower_phase_name", "" );
        tParameterList.insert( "side_ordinals", "" );
        tParameterList.insert( "neighbor_phases", "" );

        tParameterList.insert( "analytical_jacobian", true );

        tParameterList.insert( "time_continuity", false );
        tParameterList.insert( "time_boundary", false );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

    /**
     * creates an IQI parameter list with default inputs
     * @param [ out ] ParameterList a IQI parameter list
     */
    inline Parameter_List
    create_IQI_parameter_list()
    {
        Parameter_List tParameterList;

        tParameterList.insert( "IQI_name", "undefined" );
        tParameterList.insert( "IQI_type", fem::IQI_Type::UNDEFINED, fem::IQI_Type::VOLUME, fem::IQI_Type::LINEAR_ELASTICITY_DAMAGE );
        tParameterList.insert( "IQI_bulk_type", fem::Element_Type::BULK, fem::Element_Type::BULK, fem::Element_Type::TIME_BOUNDARY );
        tParameterList.insert( "dof_quantity", "" );
        tParameterList.insert( "vectorial_field_index", -1 );

        tParameterList.insert( "leader_phase_name", "" );
        tParameterList.insert( "follower_phase_name", "" );
        tParameterList.insert( "leader_dof_dependencies", "" );
        tParameterList.insert( "follower_dof_dependencies", "" );
        tParameterList.insert( "leader_dv_dependencies", "" );
        tParameterList.insert( "follower_dv_dependencies", "" );
        tParameterList.insert( "leader_field_types", "" );
        tParameterList.insert( "follower_field_types", "" );

        tParameterList.insert( "function_parameters", "" );

        tParameterList.insert( "leader_properties", "" );
        tParameterList.insert( "follower_properties", "" );
        tParameterList.insert( "leader_constitutive_models", "" );
        tParameterList.insert( "follower_constitutive_models", "" );
        tParameterList.insert( "stabilization_parameters", "" );

        tParameterList.insert( "leader_physics", "" );
        tParameterList.insert( "follower_physics", "" );
        tParameterList.insert( "stabilizations", "" );

        tParameterList.insert( "leader_phase_name", "" );
        tParameterList.insert( "follower_phase_name", "" );
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
    inline Parameter_List
    create_phase_parameter_list()
    {
        Parameter_List tParameterList;

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
    inline Parameter_List
    create_computation_parameter_list()
    {
        Parameter_List tParameterList;

        // bool true for printing physics
        tParameterList.insert( "print_physics_model", false );

        // bool true for analytical forward analysis, false for finite difference
        // decide if dRdu (jacobian) and dQIdu are computed by A/FD
        tParameterList.insert( "is_analytical_forward", true );

        // enum for finite difference scheme for forward analysis
        tParameterList.insert( "finite_difference_scheme_forward", fem::FDScheme_Type::POINT_1_FORWARD, fem::FDScheme_Type::POINT_1_FORWARD, fem::FDScheme_Type::POINT_5 );

        // real for relative perturbation size for finite difference for forward analysis
        tParameterList.insert( "finite_difference_perturbation_size_forward", 1e-6 );

        // integration order that should be used for the raytracing of nonconformal side sets
        tParameterList.insert( "nonconformal_integration_order", static_cast< uint >( mtk::Integration_Order::UNDEFINED ) );

        // Maximum length of a ray in the negative direction (w.r.t the outward pointing normal) that is used for the raytracing of nonconformal side sets.
        // This prevents that rays of thin-walled bodies get mapped "through" the body itself on the other side.
        // E.g. if you have two thin bars that are positioned on top of each other, the rays of the top-side of the upper bar should not be mapped to the lower bar.
        // Setting this value to a very small value might lead to bodies penetrating each other without being detected (especially for large load stepping increments).
        tParameterList.insert( "nonconformal_max_negative_ray_length", -0.05 );

        // determine the maximum ray length in the positive direction. This prevents mapping of elements that are too far apart to be even considered as neighbors.
        tParameterList.insert( "nonconformal_max_positive_ray_length", 0.1 );

        // bool true for analytical sensitivity analysis, false for finite difference
        // decide if dRdp and dQIdp are computed by A/FD
        tParameterList.insert( "is_analytical_sensitivity", false );

        // enum for finite difference scheme for sensitivity analysis
        tParameterList.insert( "finite_difference_scheme", fem::FDScheme_Type::POINT_1_FORWARD, fem::FDScheme_Type::POINT_1_FORWARD, fem::FDScheme_Type::POINT_5 );

        // real for relative perturbation size for finite difference
        tParameterList.insert( "finite_difference_perturbation_size", 1e-6 );

        // enum for finite difference perturbation strategy (relative, absolute)
        tParameterList.insert( "finite_difference_perturbation_strategy", fem::Perturbation_Type::RELATIVE, fem::Perturbation_Type::RELATIVE, fem::Perturbation_Type::ABSOLUTE );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm

#endif /* PROJECTS_PRM_SRC_FN_PRM_FEM_PARAMETERS_HPP_ */
