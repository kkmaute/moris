/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Comsol_conform.cpp
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
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

#ifdef __cplusplus
extern "C" {

// global variables
extern uint gTestCaseIndex;
extern uint gInterpolationOrder;
extern bool gPrintReferenceValues;

#endif
//------------------------------------------------------------------------------
namespace moris
{

    // Geometry Parameters
    // moris::real tXlength = 0.28;
    moris::real tXlength = 0.028;
    moris::real tYlength = 0.0015;
    moris::real tXcenter = 0.5 * tXlength;
    moris::real tYcenter = 0.5 * tYlength;
    moris::real tEps     = 1.0e-4;

    // mesh
    // std::string tNumElemsPerDim = "500,1";
    // std::string tHMRDomainDimensions = "0.28,0.0005";
    std::string tNumElemsPerDim      = "50,1";
    std::string tHMRDomainDimensions = "0.028,0.0005";

    // time solver parameters
    moris::sint tStep = 10;
    moris::real tTmax = 480.0;

    // ramp up of Dirichlet BC (number of time slabs to ramp up the value on the BC)
    moris::real tRampUp = 6.0;

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void
    Func_Initial_Condition(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = { { 313.0 } };
    }

    void
    Func_Wall_Condition(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        real tT = aFIManager->get_IP_geometry_interpolator()->valt()( 0 );

        real value = std::min( 350.0, 313.0 + ( 350 - 313 ) * tT / ( tRampUp * tTmax / tStep ) );

        aPropMatrix         = { { 0.0 } };
        aPropMatrix( 0, 0 ) = value;
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    // Dummy level set function for conform mesh
    moris::real Dummy_LS(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const Vector< real > & aGeometryParameters )
    {
        moris::real tLSval = aCoordinates( 0 ) + 1.0;

        // clean return value to return non-zero value
        return tLSval;
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tHMRDomainDimensions );
        aParameterLists( 0 ).set( "domain_offset", "0.0,0.0" );
        aParameterLists( 0 ).set( "domain_sidesets", "1,2,3,4" );
        aParameterLists( 0 ).set( "lagrange_output_meshes", "0" );

        switch ( gInterpolationOrder )
        {
            case 1:
            {
                aParameterLists( 0 ).set( "lagrange_orders", "1" );
                aParameterLists( 0 ).set( "bspline_orders", "1" );
                break;
            }
            case 2:
            {
                aParameterLists( 0 ).set( "lagrange_orders", "2" );
                aParameterLists( 0 ).set( "bspline_orders", "2" );
                break;
            }
            default:
            {
                MORIS_ERROR( false, "EXA::Comsol_conform: This 2D Example can only be run with Linear or Quadratic" );
            }
        }

        aParameterLists( 0 ).set( "lagrange_pattern", "0" );
        aParameterLists( 0 ).set( "bspline_pattern", "0" );

        aParameterLists( 0 ).set( "lagrange_to_bspline", "0" );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );
        aParameterLists( 0 ).set( "refinement_buffer", 3 );
        aParameterLists( 0 ).set( "staircase_buffer", 3 );
        aParameterLists( 0 ).set( "initial_refinement", "0" );
        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        aParameterLists( 0 ).set( "adaptive_refinement_level", 1 );
    }

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        aParameterLists( 0 ).set( "is_optimization_problem", false );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", "conformal" );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", "bspline" );
        aParameterLists( 0 ).set( "enrich_mesh_indices", "0" );
        aParameterLists( 0 ).set( "ghost_stab", true );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );

        // Geometry parameter lists
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Dummy_LS" );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", "0.75" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropConductivity" );
        aParameterLists( 0 ).set( "function_parameters", "2.1e-7" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropHeatCapacity" );
        aParameterLists( 0 ).set( "function_parameters", "2.4" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 3
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropLatentHeat" );
        aParameterLists( 0 ).set( "function_parameters", "175" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 9
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPCTemp" );
        aParameterLists( 0 ).set( "function_parameters", "314.5" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPhaseState" );
        aParameterLists( 0 ).set( "function_parameters", "2.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 11
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropPCconst" );
        aParameterLists( 0 ).set( "function_parameters", "3.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropImposedTemp" );
        aParameterLists( 0 ).set( "value_function", "Func_Wall_Condition" );

        // create parameter list for property 6
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightCurrent" );
        aParameterLists( 0 ).set( "function_parameters", "100.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 7
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropWeightPrevious" );
        aParameterLists( 0 ).set( "function_parameters", "100.0" );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 8
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropInitialCondition" );
        aParameterLists( 0 ).set( "value_function", "Func_Initial_Condition" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        aParameterLists( 1 ).set( "constitutive_name", "CMDiffusionPhaseChange" );
        aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO_PC ) ;
        aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 1 ).set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity,Density;"
                "PropHeatCapacity,HeatCapacity;"
                "PropLatentHeat,LatentHeat;"
                "PropPCTemp,PCTemp;"
                "PropPhaseState,PhaseStateFunction;"
                "PropPCconst,PhaseChangeConst" );

        //------------------------------------------------------------------------------

        // Ghost
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGPTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GHOST_DISPL ) ;
        aParameterLists( 2 ).set( "function_parameters", "0.01" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );

        // GGLS parameter
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPGGLSDiffusion" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::GGLS_DIFFUSION ) ;
        aParameterLists( 2 ).set( "leader_dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists( 2 ).set( "leader_properties",
                "PropConductivity , Conductivity;"
                "PropDensity      , Density;"
                "PropHeatCapacity , HeatCapacity;"
                "PropLatentHeat   , LatentHeat;"
                "PropPCTemp       , PCTemp;"
                "PropPhaseState   , PhaseStateFunction;"
                "PropPCconst      , PhaseChangeConst" );

        // Dirichlet SP
        aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        aParameterLists( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists( 2 ).set( "function_parameters", "1000.0" );
        aParameterLists( 2 ).set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // Bulk
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPGGLSDiffusion,GGLSParam" );
        aParameterLists( 3 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        // Dirichlet BC
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGOutletTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties", "PropImposedTemp,Dirichlet" );
        aParameterLists( 3 ).set( "leader_constitutive_models", "CMDiffusionPhaseChange,Diffusion" );
        aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists( 3 ).set( "mesh_set_names", "SideSet_4_n_p1" );

        // Time Continuity
        aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        aParameterLists( 3 ).set( "IWG_name", "IWGTimeContinuityTemp" );
        aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
        aParameterLists( 3 ).set( "dof_residual", "TEMP" );
        aParameterLists( 3 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 3 ).set( "leader_properties",
                "PropWeightCurrent,WeightCurrent;"
                "PropWeightPrevious,WeightPrevious;"
                "PropInitialCondition,InitialCondition" );
        aParameterLists( 3 ).set( "mesh_set_names", "HMR_dummy_n_p1" );
        aParameterLists( 3 ).set( "time_continuity", true );

        //------------------------------------------------------------------------------
        // Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "vectorial_field_index", 0 );
        aParameterLists( 4 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        // Max Temperature IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQIMaxTEMP" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::MAX_DOF ) ;
        aParameterLists( 4 ).set( "dof_quantity", "TEMP" );
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "function_parameters", "313.0/30.0" );
        aParameterLists( 4 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        // Latent Heat IQI
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQILatentHeatAbsorption" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::LATENT_HEAT_ABSORPTION ) ;
        aParameterLists( 4 ).set( "leader_dof_dependencies", "TEMP" );
        aParameterLists( 4 ).set( "leader_properties",
                "PropDensity,Density;"
                "PropLatentHeat,LatentHeat;"
                "PropPCTemp,PCTemp;"
                "PropPhaseState,PhaseStateFunction;"
                "PropPCconst,PhaseChangeConst" );
        aParameterLists( 4 ).set( "mesh_set_names", "HMR_dummy_n_p1" );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", 1e-7 );
        aParameterLists( 2 ).set( "NLA_relaxation_strategy",  sol::SolverRelaxationType::InvResNormAdaptive ) ;
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", 0.5 );
        aParameterLists( 2 ).set( "NLA_relaxation_damping", 0.5 );
        aParameterLists( 2 ).set( "NLA_max_iter", 40 );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "TEMP" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tStep );
        aParameterLists( 4 ).set( "TSA_Time_Frame", tTmax );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "TEMP" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists( 5 ).set( "TSA_time_level_per_type", "TEMP,2" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "Comsol_conform_" + std::to_string( gTestCaseIndex ) + ".exo" ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", "HMR_dummy_n_p1" );
        aParameterLists( 0 ).set( "Field_Names", "TEMP,MAX_DOF,LATENT_HEAT_ABSORPTION" );
        aParameterLists( 0 ).set( "Field_Type", "NODAL,GLOBAL,GLOBAL" );
        aParameterLists( 0 ).set( "IQI_Names", "IQIBulkTEMP,IQIMaxTEMP,IQILatentHeatAbsorption" );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }
    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
