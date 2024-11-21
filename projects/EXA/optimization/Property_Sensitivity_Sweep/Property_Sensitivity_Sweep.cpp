/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * SIMP.cpp
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
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // FD in adjoint
    real tFEMFdEpsilon = 1.0e-7;

    // FD step size in sweep
    std::string tFDsweep = "1.0e-7";

    std::string tMeshSets = "HMR_dummy_n_p0";

    //--------------------------------------------------------------------------------------------------------------

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tYoungsFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        const moris::Matrix< moris::DDRMat >& tHCT  = aParameters( 0 );
        moris::real                           tBeta = aParameters( 0 )( 1 );
        aPropMatrix                                 = tHCT * std::pow( aFIManager->get_field_interpolators_for_type( gen::PDV_Type::DENSITY )->val()( 0 ), tBeta );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tDerYoungsFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        const moris::Matrix< moris::DDRMat >& tHCT  = aParameters( 0 );
        moris::real                           tBeta = aParameters( 0 )( 1 );
        aPropMatrix                                 = tBeta * tHCT * std::pow( aFIManager->get_field_interpolators_for_type( gen::PDV_Type::DENSITY )->val()( 0 ), tBeta - 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tDensityFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::DENSITY )->val();
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    tDerDensityFunc(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aFIManager->get_field_interpolators_for_type( gen::PDV_Type::DENSITY )->N();
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris::real
    Func_Right_Plane(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        return aCoordinates( 0 ) - 1000;
    }

    //--------------------------------------------------------------------------------------------------------------

    moris::real
    sensitivity_func(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        return 0;
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
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );
        tObjectives( 0 ) = aCriteria( 0 ) / 2.7;

        return tObjectives;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = aCriteria( 1 ) - 0.8;

        return tConstraints;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );

        return tDObjectiveDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria( 1, 2 );
        tDObjectiveDCriteria( 0 ) = 1.0 / 2.7;
        tDObjectiveDCriteria( 1 ) = 0;

        return tDObjectiveDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );

        return tDConstraintDADV;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria( 1, 2 );
        tDConstraintDCriteria( 0 ) = 0;
        tDConstraintDCriteria( 1 ) = 1.0;

        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", "2,2" );
        aParameterLists.set( "domain_dimensions", "1,1" );
        aParameterLists.set( "domain_offset", "0.0,0.0" );
        aParameterLists.set( "domain_sidesets", "1,2,3,4" );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "truncate_bsplines", 1 );
        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );
        aParameterLists.set( "initial_refinement", "0" );

        aParameterLists.set( "use_multigrid", 0 );
        aParameterLists.set( "severity_level", 1 );

        aParameterLists.set( "adaptive_refinement_level", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );
        aParameterLists.set( "PDV_types", "DENSITY" );

        // Geometry parameter lists (dummy geometry)
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Right_Plane" );

        // Density property
        aParameterLists( GEN::PROPERTIES ).add_parameter_list( gen::Field_Type::CONSTANT );
        aParameterLists.set( "name", "density" );
        aParameterLists.set( "constant", 0.4 );
        aParameterLists.set( "discretization_mesh_index", 0 );
        aParameterLists.set( "discretization_lower_bound", 0.001 );
        aParameterLists.set( "discretization_upper_bound", 1.0 );
        aParameterLists.set( "pdv_type", "DENSITY" );
        aParameterLists.set( "pdv_mesh_set_names", "HMR_dummy_n_p0,SideSet_4_n_p0,SideSet_2_n_p0" );
        aParameterLists.set( "pdv_mesh_set_indices", "" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();

        //------------------------------------------------------------------------------

        // create parameter list for property 1
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "tDensityFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerDensityFunc" );
        aParameterLists.set( "dv_dependencies", "DENSITY" );

        // create parameter list for property 2
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "1.0,3.0" );
        aParameterLists.set( "value_function", "tYoungsFunc" );
        aParameterLists.set( "dv_derivative_functions", "tDerYoungsFunc" );
        aParameterLists.set( "dv_dependencies", "DENSITY" );

        // create parameter list for property 5
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropFlux" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 4
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 10
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropNeumann" );
        aParameterLists.set( "function_parameters", "1.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        // create parameter list for property 7
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso1" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------

        // create parameter list for stabilization parameter 1
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheTemp" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        //------------------------------------------------------------------------------

        // create parameter list for IWG 1
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU_1" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", "HMR_dummy_n_p0" );

        // create parameter list for IWG 2
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheTemp,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "SideSet_4_n_p0" );

        // create parameter list for IWG 3
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGNeumannFlux" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropNeumann,Traction" );
        aParameterLists.set( "mesh_set_names", "SideSet_2_n_p0" );

        //------------------------------------------------------------------------------

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDispX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tMeshSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDispY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tMeshSets );

        /*
                aParameterLists( FEM::IQI ).add_parameter_list();
                aParameterLists.set( "IQI_name", "IQINodeX" );
                aParameterLists.set( "IQI_type", fem::IQI_Type::NODE_SENSITIVITY );
                aParameterLists.set( "vectorial_field_index", 0 );
                aParameterLists.set( "function_parameters", "0" );
                aParameterLists.set( "mesh_set_names", tMeshSets );

                aParameterLists( FEM::IQI ).add_parameter_list();
                aParameterLists.set( "IQI_name", "IQINodeY" );
                aParameterLists.set( "IQI_type", fem::IQI_Type::NODE_SENSITIVITY );
                aParameterLists.set( "vectorial_field_index", 1 );
                aParameterLists.set( "function_parameters", "0" );
                aParameterLists.set( "mesh_set_names", tMeshSets );

                aParameterLists( FEM::IQI ).add_parameter_list();
                aParameterLists.set( "IQI_name", "IQIDispXDV" );
                aParameterLists.set( "IQI_type", fem::IQI_Type::DOF_SENSITIVITY );
                aParameterLists.set( "dof_quantity", "UX,UY" );
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
                aParameterLists.set( "vectorial_field_index", 0 );
                aParameterLists.set( "function_parameters", "0" );
                aParameterLists.set( "mesh_set_names", tMeshSets );

                aParameterLists( FEM::IQI ).add_parameter_list();
                aParameterLists.set( "IQI_name", "IQIDispYDV" );
                aParameterLists.set( "IQI_type", fem::IQI_Type::DOF_SENSITIVITY );
                aParameterLists.set( "dof_quantity", "UX,UY" );
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
                aParameterLists.set( "vectorial_field_index", 1 );
                aParameterLists.set( "function_parameters", "0" );
                aParameterLists.set( "mesh_set_names", tMeshSets );
        */
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso1,Elast" );
        aParameterLists.set( "mesh_set_names", tMeshSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tMeshSets );

        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
        aParameterLists.set( "finite_difference_perturbation_size", tFEMFdEpsilon );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", 1 );
        aParameterLists.set( "TSA_Time_Frame", 1.0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::SOLVER_WAREHOUSE );
        aParameterLists.set( "Sensitivity_Analysis_Type", sol::SensitivityAnalysisType::ADJOINT );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Property_Sensitivity_Sweep.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tMeshSets );
        aParameterLists.set( "Field_Names", "UX,UY" );              //,DPX,DPY,DUX,DUY" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL" );         //,NODAL,NODAL,NODAL,NODAL" );
        aParameterLists.set( "IQI_Names", "IQIDispX,IQIDispY" );    //,IQINodeX,IQINodeY,IQIDispXDV,IQIDispYDV" );
        // aParameterLists.set( "Analysis_Type", "FORWARD,FORWARD,SENSITIVITY,SENSITIVITY,SENSITIVITY,SENSITIVITY" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", "./Property_Sensitivity_Sweep.so" );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "hdf5_path", "property_opt_test.hdf5" );
        aParameterLists.set( "evaluate_objective_gradients", true );
        aParameterLists.set( "evaluate_constraint_gradients", true );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "include_bounds", false );
        aParameterLists.set( "finite_difference_type", "all" );
        aParameterLists.set( "finite_difference_epsilons", tFDsweep );
    }

    //--------------------------------------------------------------------------------------------------------------

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
