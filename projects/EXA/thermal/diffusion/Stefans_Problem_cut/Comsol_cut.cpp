/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Comsol_cut.cpp
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {
#endif

namespace moris
{
    // global variables
    extern uint gInterpolationOrder;

    // Geometry Parameters
    real tXlength = 0.028;
    real tYlength = 0.0015;
    real tXcenter = 0.5 * tXlength;
    real tYcenter = 0.5 * tYlength;
    real tEps     = 1.0e-4;

    // mesh
    // std::string tNumElemsPerDim = "520,6";
    // std::string tHMRDomainDimensions = "0.2912,0.0030";
    std::string tNumElemsPerDim      = "26,3";
    std::string tHMRDomainDimensions = "0.030,0.0030";

    // time solver parameters
    sint tStep = 15;
    real tTmax = 480.0;

    // ramp up of Dirichlet BC (number of time slabs to ramp up the value on the BC)
    real tRampUp = 9.0;

    // Constant function for properties
    void
    Func_Const(
            Matrix< DDRMat >&                aPropMatrix,
            Vector< Matrix< DDRMat > >&      aParameters,
            fem::Field_Interpolator_Manager* aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void
    Func_Initial_Condition(
            Matrix< DDRMat >&                aPropMatrix,
            Vector< Matrix< DDRMat > >&      aParameters,
            fem::Field_Interpolator_Manager* aFIManager )
    {
        aPropMatrix = { { 313.0 } };
    }

    void
    Func_Wall_Condition(
            Matrix< DDRMat >&                aPropMatrix,
            Vector< Matrix< DDRMat > >&      aParameters,
            fem::Field_Interpolator_Manager* aFIManager )
    {
        real tT = aFIManager->get_IP_geometry_interpolator()->valt()( 0 );

        real value = std::min( 350.0, 313.0 + ( 350 - 313 ) * tT / ( tRampUp * tTmax / tStep ) );

        aPropMatrix         = { { 0.0 } };
        aPropMatrix( 0, 0 ) = value;
    }

    bool
    Output_Criterion( tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    // Level set functions for Rectangle
    real Top_Boundary(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        real tLSval = 0.5 * tYlength - aCoordinates( 1 ) + tYcenter;

        // clean return value to return non-zero value
        return tLSval;
    }

    real Bottom_Boundary(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        real tLSval = 0.5 * tYlength + aCoordinates( 1 ) - tYcenter;

        // clean return value to return non-zero value
        return tLSval;
    }

    real Left_Boundary(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        real tLSval = 0.5 * tXlength + aCoordinates( 0 ) - tXcenter;

        // clean return value to return non-zero value
        return tLSval;
    }

    real Right_Boundary(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        real tLSval = 0.5 * tXlength - aCoordinates( 0 ) + tXcenter;

        // clean return value to return non-zero value
        return tLSval;
    }

    void
    HMRParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 ).set( "domain_dimensions", tHMRDomainDimensions );
        tParameterlist( 0 ).set( "domain_offset", "-0.0015,-0.0008" );
        tParameterlist( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 ).set( "lagrange_output_meshes", "0" );

        switch ( gInterpolationOrder )
        {
            case 1:
            {
                tParameterlist( 0 ).set( "lagrange_orders", "1" );
                tParameterlist( 0 ).set( "bspline_orders", "1" );
                break;
            }
            case 2:
            {
                tParameterlist( 0 ).set( "lagrange_orders", "2" );
                tParameterlist( 0 ).set( "bspline_orders", "2" );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "EXA::Comsol_cut: This 2D Example can only be run with Linear or Quadratic" );
            }
        }

        tParameterlist( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 ).set( "refinement_buffer", 3 );
        tParameterlist( 0 ).set( "staircase_buffer", 3 );
        tParameterlist( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 ).set( "adaptive_refinement_level", 1 );
    }

    void
    OPTParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        tParameterlist( 0 ).set( "is_optimization_problem", false );
    }

    void
    XTKParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        tParameterlist( 0 ).set( "decompose", true );
        tParameterlist( 0 ).set( "decomposition_type", "conformal" );
        tParameterlist( 0 ).set( "enrich", true );
        tParameterlist( 0 ).set( "basis_rank", "bspline" );
        tParameterlist( 0 ).set( "enrich_mesh_indices", "0" );
        tParameterlist( 0 ).set( "ghost_stab", true );
        tParameterlist( 0 ).set( "multigrid", false );
        tParameterlist( 0 ).set( "verbose", true );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( prm::create_gen_parameter_list() );


        // Geometry parameter lists
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Bottom_Boundary" );

        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Top_Boundary" );

        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Right_Boundary" );

        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name", "Left_Boundary" );
    }

    void
    FEMParameterList( Vector< Submodule_Parameter_Lists >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropDensity" );
        tParameterList( 0 ).set( "function_parameters", "0.75" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropConductivity" );
        tParameterList( 0 ).set( "function_parameters", "2.1e-7" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 3
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropHeatCapacity" );
        tParameterList( 0 ).set( "function_parameters", "2.4" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 3
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropLatentHeat" );
        tParameterList( 0 ).set( "function_parameters", "175" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 9
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropPCTemp" );
        tParameterList( 0 ).set( "function_parameters", "314.5" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropPhaseState" );
        tParameterList( 0 ).set( "function_parameters", "2.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 11
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropPCconst" );
        tParameterList( 0 ).set( "function_parameters", "3.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 4
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropImposedTemp" );
        tParameterList( 0 ).set( "value_function", "Func_Wall_Condition" );

        // create parameter list for property 6
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropWeightCurrent" );
        tParameterList( 0 ).set( "function_parameters", "100.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 7
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropWeightPrevious" );
        tParameterList( 0 ).set( "function_parameters", "100.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 8
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropInitialCondition" );
        tParameterList( 0 ).set( "value_function", "Func_Initial_Condition" );

        //------------------------------------------------------------------------------


        // create parameter list for constitutive model 1
        tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMDiffusionPhaseChange" );
        tParameterList( 1 ).set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO_PC );
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity,Density;"
                "PropHeatCapacity,HeatCapacity;"
                "PropLatentHeat,LatentHeat;"
                "PropPCTemp,PCTemp;"
                "PropPhaseState,PhaseStateFunction;"
                "PropPCconst,PhaseChangeConst" );

        //------------------------------------------------------------------------------


        // Ghost
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPGPTemp" );
        tParameterList( 2 ).set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        tParameterList( 2 ).set( "function_parameters", "0.01" );
        tParameterList( 2 ).set( "leader_properties", "PropConductivity,Material" );

        // GGLS parameter
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPGGLSDiffusion" );
        tParameterList( 2 ).set( "stabilization_type", fem::Stabilization_Type::GGLS_DIFFUSION );
        tParameterList( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 2 ).set( "leader_properties",
                "PropConductivity , Conductivity;"
                "PropDensity      , Density;"
                "PropHeatCapacity , HeatCapacity;"
                "PropLatentHeat   , LatentHeat;"
                "PropPCTemp       , PCTemp;"
                "PropPhaseState   , PhaseStateFunction;"
                "PropPCconst      , PhaseChangeConst" );

        // Dirichlet SP
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( 2 ).set( "function_parameters", "1000.0" );
        tParameterList( 2 ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // Bulk
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGDiffusionBulk" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPGGLSDiffusion,GGLSParam" );
        tParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p15,HMR_dummy_c_p15" );

        // Dirichlet BC
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGOutletTemp" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "leader_properties", "PropImposedTemp,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", "iside_b0_15_b1_14" );

        // Time Continuity
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGTimeContinuityTemp" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::TIME_CONTINUITY_DOF );
        tParameterList( 3 ).set( "dof_residual", "TEMP" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 3 ).set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        tParameterList( 3 ).set( "mesh_set_names", "HMR_dummy_n_p15,HMR_dummy_c_p15" );
        tParameterList( 3 ).set( "time_continuity", true );

        //------------------------------------------------------------------------------
        // create parameter list for IQI 4
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 ).set( "dof_quantity", "TEMP" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "TEMP" );
        tParameterList( 4 ).set( "vectorial_field_index", 0 );
        tParameterList( 4 ).set( "mesh_set_names", "HMR_dummy_n_p15,HMR_dummy_c_p15" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 8 );


        tParameterlist( 0 ).add_parameter_list( prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        tParameterlist( 1 ).add_parameter_list( prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).add_parameter_list( prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 ).set( "NLA_rel_res_norm_drop", 2.0e-05 );
        tParameterlist( 2 ).set( "NLA_relaxation_parameter", 0.96 );
        tParameterlist( 2 ).set( "NLA_max_iter", 20 );

        tParameterlist( 3 ).add_parameter_list( prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_DofTypes", "TEMP" );

        tParameterlist( 4 ).add_parameter_list( prm::create_time_solver_algorithm_parameter_list() );
        tParameterlist( 4 ).set( "TSA_Num_Time_Steps", tStep );
        tParameterlist( 4 ).set( "TSA_Time_Frame", tTmax );

        tParameterlist( 5 ).add_parameter_list( prm::create_time_solver_parameter_list() );
        tParameterlist( 5 ).set( "TSA_DofTypes", "TEMP" );
        tParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        tParameterlist( 5 ).set( "TSA_time_level_per_type", "TEMP,2" );

        tParameterlist( 6 ).add_parameter_list( prm::create_solver_warehouse_parameterlist() );

        tParameterlist( 7 ).add_parameter_list( prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Comsol_cut.exo" ) );
        tParameterlist( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 ).set( "Set_Names", "HMR_dummy_n_p15,HMR_dummy_c_p15" );
        tParameterlist( 0 ).set( "Field_Names", "TEMP" );
        tParameterlist( 0 ).set( "Field_Type", "NODAL" );
        tParameterlist( 0 ).set( "IQI_Names", "IQIBulkTEMP" );
    }

    void
    MORISGENERALParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
