/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Shape_Sensitivity_Sweep.cpp
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
//------------------------------------------------------------------------------
namespace moris
{
    std::string tMeshSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tDBCSets  = "iside_b0_1_b1_3";
    std::string tNBCSets  = "iside_b0_1_b1_0";

    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-7;

    // FD step size in sweep
    std::string tFDsweep = "1.0e-7";

    //--------------------------------------------------------------------------------------------------------------

    // Constant function for properties
    void
    Func_Const( moris::Matrix< moris::DDRMat >       &aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > > &aParameters,
            moris::fem::Field_Interpolator_Manager   *aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    Output_Criterion( moris::tsa::Time_Solver *aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );

        return tConstraintTypes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_objectives( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tObjectives = { { aCriteria( 0 ) } };

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tConstraints = { { aCriteria( 1 ) } };

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, 2 );
        tDObjectiveDCriteria( 0 ) = 1;
        tDObjectiveDCriteria( 1 ) = 0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( const Vector< real > &aADVs, const Vector< real > &aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, 2 );
        tDConstraintDCriteria( 0 ) = 0;
        tDConstraintDCriteria( 1 ) = 1.0;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", "2, 2" );
        tParameterlist( 0 ).set( "domain_dimensions", "2.0, 2.0" );
        tParameterlist( 0 ).set( "domain_offset", "-1.0, -1.0" );
        tParameterlist( 0 ).set( "domain_sidesets", "1,2,3,4" );
        tParameterlist( 0 ).set( "lagrange_output_meshes", "0" );

        tParameterlist( 0 ).set( "lagrange_orders", "1" );
        tParameterlist( 0 ).set( "lagrange_pattern", "0" );
        tParameterlist( 0 ).set( "bspline_orders", "1" );
        tParameterlist( 0 ).set( "bspline_pattern", "0" );

        tParameterlist( 0 ).set( "lagrange_to_bspline", "0" );

        tParameterlist( 0 ).set( "truncate_bsplines", 1 );
        tParameterlist( 0 ).set( "refinement_buffer", 3 );
        tParameterlist( 0 ).set( "staircase_buffer", 3 );
        tParameterlist( 0 ).set( "initial_refinement", "0" );
        tParameterlist( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 ).set( "use_multigrid", 0 );

        tParameterlist( 0 ).set( "adaptive_refinement_level", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
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
        tParameterlist( 0 ).set( "verbose", false );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( moris::prm::create_gen_parameter_list() );
        tParameterlist( 0 ).set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );

        // Geometry parameter lists

        // vertical line
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        tParameterlist( 1 ).set( "center_x", 0.8, 0.8, 0.8 );
        tParameterlist( 1 ).set( "center_y", 0.3 );
        tParameterlist( 1 ).set( "normal_x", 1.0 );
        tParameterlist( 1 ).set( "normal_y", 0.0, 0.0, 0.0 );

        // oblique line
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
        tParameterlist( 1 ).set( "center_x", -0.6, -0.6, -0.6 );
        tParameterlist( 1 ).set( "center_y", -0.3, -0.3, -0.3 );
        tParameterlist( 1 ).set( "normal_x", 1.0, 1.0, 1.0 );
        tParameterlist( 1 ).set( "normal_y", 1.0, 1.0, 1.0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Vector< Submodule_Parameter_Lists > &tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropDensity" );
        tParameterList( 0 ).set( "function_parameters", "1.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 2
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropYoungs" );
        tParameterList( 0 ).set( "function_parameters", "1.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 5
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropFlux" );
        tParameterList( 0 ).set( "function_parameters", "10.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 4
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropDirichletU" );
        tParameterList( 0 ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 10
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropTraction" );
        tParameterList( 0 ).set( "function_parameters", "0.0;0.1" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        // create parameter list for property 7
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropPoisson" );
        tParameterList( 0 ).set( "function_parameters", "0.0" );
        tParameterList( 0 ).set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMStrucLinIso1" );
        tParameterList( 1 ).set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        tParameterList( 1 ).set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPNitscheTemp" );
        tParameterList( 2 ).set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        tParameterList( 2 ).set( "function_parameters", "100.0" );
        tParameterList( 2 ).set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------
        // create parameter list for IWG 1
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGBulkU_1" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        tParameterList( 3 ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 ).set( "mesh_set_names", tMeshSets );

        // create parameter list for IWG 2
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGDirichletU" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        tParameterList( 3 ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 ).set( "leader_properties", "PropDirichletU,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", tDBCSets );

        // create parameter list for IWG 3
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGTraction" );
        tParameterList( 3 ).set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        tParameterList( 3 ).set( "dof_residual", "UX,UY" );
        tParameterList( 3 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 3 ).set( "leader_properties", "PropTraction,Traction" );
        tParameterList( 3 ).set( "mesh_set_names", tNBCSets );

        //------------------------------------------------------------------------------
        // create parameter list for IQI 4
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIDispX" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 ).set( "vectorial_field_index", 0 );
        tParameterList( 4 ).set( "mesh_set_names", tMeshSets );

        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIDispY" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::DOF );
        tParameterList( 4 ).set( "dof_quantity", "UX,UY" );
        tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 ).set( "vectorial_field_index", 1 );
        tParameterList( 4 ).set( "mesh_set_names", tMeshSets );

        // create parameter list for IQI 4
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkStrainEnergy" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        tParameterList( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        tParameterList( 4 ).set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        tParameterList( 4 ).set( "mesh_set_names", tMeshSets );

        // create parameter list for IQI 4
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQIBulkVolume" );
        tParameterList( 4 ).set( "IQI_type", fem::IQI_Type::VOLUME );
        tParameterList( 4 ).set( "leader_properties", "PropDensity,Density" );
        tParameterList( 4 ).set( "mesh_set_names", tMeshSets );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );

        tParameterList( 5 ).set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        tParameterList( 5 ).set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 8 );


        tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
        tParameterlist( 0 ).set( "preconditioners", "0" );

        tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set( "NLA_DofTypes", "UX,UY" );

        tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        tParameterlist( 4 ).set( "TSA_Num_Time_Steps", 1 );
        tParameterlist( 4 ).set( "TSA_Time_Frame", 1.0 );

        tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 ).set( "TSA_DofTypes", "UX,UY" );
        tParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        tParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        tParameterlist( 6 ).set( "Sensitivity_Analysis_Type", sol::SensitivityAnalysisType::ADJOINT );

        tParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        tParameterlist( 7 ).set( "ifpack_prec_type", "ILU" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", "shape_sensitivities.exo" ) );
        tParameterlist( 0 ).set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        tParameterlist( 0 ).set( "Set_Names", tMeshSets );
        tParameterlist( 0 ).set( "Field_Names", "UX,UY" );
        tParameterlist( 0 ).set( "Field_Type", "NODAL,NODAL" );
        tParameterlist( 0 ).set( "IQI_Names", "IQIDispX,IQIDispY" );
        tParameterlist( 0 ).set( "Save_Frequency", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 ).set( "is_optimization_problem", true );
        tParameterlist( 0 ).set( "problem", "user_defined" );
        tParameterlist( 0 ).set( "library", "Shape_Sensitivity_Sweep.so" );

        tParameterlist( 2 ).add_parameter_list( moris::prm::create_sweep_parameter_list() );
        tParameterlist( 2 ).set( "hdf5_path", "shape_opt_test.hdf5" );
        tParameterlist( 2 ).set( "evaluate_objective_gradients", true );
        tParameterlist( 2 ).set( "evaluate_constraint_gradients", true );
        tParameterlist( 2 ).set( "num_evaluations_per_adv", "1" );
        tParameterlist( 2 ).set( "include_bounds", false );
        tParameterlist( 2 ).set( "finite_difference_type", "all" );
        tParameterlist( 2 ).set( "finite_difference_epsilons", tFDsweep );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MORISGENERALParameterList( Vector< Submodule_Parameter_Lists > &tParameterlist )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
