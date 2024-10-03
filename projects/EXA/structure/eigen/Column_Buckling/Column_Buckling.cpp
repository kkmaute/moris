/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * single_element.cpp
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
// #include "cl_DLA_Eigen_Solver.hpp"
#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_MSI_Parameters.hpp"
#include "fn_PRM_SOL_Parameters.hpp"
#include "fn_PRM_VIS_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "fn_PRM_OPT_Parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

//------------------------------------------------------------------------------

// interpolation order
extern std::string tOrder;

// spatial dimensions
extern uint tDim;

// stress used in geometric stiffness (PK2 or Linear)
extern std::string tStressType;

// output file
extern std::string tOutputFileName;

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------------

namespace moris
{
    void create_trilinos_parameter_list( Vector< Submodule_Parameter_Lists >& );
    void create_petsc_parameter_list( Vector< Submodule_Parameter_Lists >& );

    /* ------------------------------------------------------------------------ */
    // General

    // use Petsc or Trilinos eigensolver (Trilinos does not work for this problem yet)
    bool tUsePetsc = true;

    // Traction load
    std::string tTraction = tDim == 3 ? "0.0;-1.0;0.0" : "0.0;-1.0";

    // Prescribed displacements
    std::string tDirichlet = tDim == 3 ? "0.0;0.0;0.0" : "0.0;0.0";

    // Number of eigenvectors
    sint tNumEigenVectors = 5;

    std::string tDofString = tDim == 3 ? "UX,UY,UZ" : "UX,UY";

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    std::string tNumElemsPerDim = tDim == 3 ? "5,50,5" : "5,50";
    std::string tDomainDims     = tDim == 3 ? "1,10,1" : "1,10";
    std::string tDomainOffset   = tDim == 3 ? "0.0,0.0,0.0" : "0.0,0.0";
    std::string tDomainSidesets = tDim == 3 ? "1,2,3,4,5,6" : "1,2,3,4";

    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    // Bulk Phases
    std::string tBulk = "HMR_dummy_n_p0";

    // boundaries
    std::string tBottom = "SideSet_1_n_p0";
    std::string tTop    = "SideSet_3_n_p0";

    /* ------------------------------------------------------------------------ */
    // material parameters, kg is scaled with a factor 1e-6

    // Copper material properties
    std::string tYoungsModulus = "100";
    std::string tPoissonRatio  = "0.3";

    /* ------------------------------------------------------------------------ */
    // Output Config
    bool tOutputCriterion = true;

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop    = 1.0e-10;
    moris::real tNLA_relaxation_parameter = 1.0;

    /* ------------------------------------------------------------------------ */
    // DUMMY FUNCTIONS
    /* ------------------------------------------------------------------------ */

    // Output criterion for VIS mesh
    bool Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return tOutputCriterion;
    }

    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void OPTParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        tParameterlist( 0 ).set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void HMRParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 ).set( "lagrange_output_meshes", std::string( "0" ) );

        tParameterlist( 0 ).set( "lagrange_orders", tOrder );
        tParameterlist( 0 ).set( "lagrange_pattern", std::string( "0" ) );
        tParameterlist( 0 ).set( "bspline_orders", tOrder );
        tParameterlist( 0 ).set( "bspline_pattern", std::string( "0" ) );

        tParameterlist( 0 ).set( "truncate_bsplines", 1 );

        tParameterlist( 0 ).set( "use_number_aura", 1 );

        tParameterlist( 0 ).set( "initial_refinement", "0" );

        tParameterlist( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 ).set( "write_lagrange_output_mesh_to_exodus", "LBAProblem.exo" );
    }

    /* ------------------------------------------------------------------------ */

    void XTKParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        tParameterlist( 0 ).set( "decompose", true );
        tParameterlist( 0 ).set( "decomposition_type", std::string( "conformal" ) );
        tParameterlist( 0 ).set( "enrich", true );
        tParameterlist( 0 ).set( "basis_rank", std::string( "bspline" ) );
        tParameterlist( 0 ).set( "enrich_mesh_indices", std::string( "0" ) );
        tParameterlist( 0 ).set( "ghost_stab", false );
        tParameterlist( 0 ).set( "multigrid", false );
        tParameterlist( 0 ).set( "verbose", true );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh", false );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void GENParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
    }

    /* ------------------------------------------------------------------------ */

    void FEMParameterList( Vector< Submodule_Parameter_Lists >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Youngs Modulus Shell
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", "PropYoungsModulus" );
        tParameterList( 0 ).set( "function_parameters", tYoungsModulus );

        // Poisson Ratio Shell
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropPoissonRatio" ) );
        tParameterList( 0 ).set( "function_parameters", tPoissonRatio );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // Compressive load
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropTraction" ) );
        tParameterList( 0 ).set( "function_parameters", tTraction );

        // Dirichlet structure
        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name", std::string( "PropDirichletStruct" ) );
        tParameterList( 0 ).set( "function_parameters", tDirichlet );

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        // linear elasticity - shell - 1
        tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMStrucLinIso" );
        tParameterList( 1 ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
        tParameterList( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofString, "Displacement" ) );
        tParameterList( 1 ).set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );

        tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMStrucNonLinIso" );
        tParameterList( 1 ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
        tParameterList( 1 ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF ) );
        tParameterList( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( tDofString, "Displacement" ) );
        tParameterList( 1 ).set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        // Displacements - Shell - back wall
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name", "SPNitscheStruc" );
        tParameterList( 2 ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        tParameterList( 2 ).set( "function_parameters", "100.0" );
        tParameterList( 2 ).set( "leader_properties", "PropYoungsModulus,Material" );

        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // linear elasticity
        ///*
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGStructBulk" );
        tParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        tParameterList( 3 ).set( "dof_residual", tDofString );
        tParameterList( 3 ).set( "leader_dof_dependencies", tDofString );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        tParameterList( 3 ).set( "mesh_set_names", tBulk );
        //*/

        //------------------------------------------------------------------------------

        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGStructBulkLBA" );
        tParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_NON_LINEAR_GEOMETRIC_STIFFNESS ) );
        tParameterList( 3 ).set( "dof_residual", tDofString );
        tParameterList( 3 ).set( "leader_dof_dependencies", tDofString );
        if ( tStressType == "PK2" )
        {
            tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucNonLinIso,ElastLinIso" );
        }
        else
        {
            tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        }
        tParameterList( 3 ).set( "mesh_set_names", tBulk );

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // Traction BC ( Compressive force at top )
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGNeumannTraction" );
        tParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        tParameterList( 3 ).set( "dof_residual", tDofString );
        tParameterList( 3 ).set( "leader_dof_dependencies", tDofString );
        tParameterList( 3 ).set( "leader_properties", "PropTraction,Traction" );
        tParameterList( 3 ).set( "mesh_set_names", tTop );

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------
        // Fixed bottom edge
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name", "IWGDirichletStruct" );
        tParameterList( 3 ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
        tParameterList( 3 ).set( "dof_residual", tDofString );
        tParameterList( 3 ).set( "leader_dof_dependencies", tDofString );
        tParameterList( 3 ).set( "leader_properties", "PropDirichletStruct,Dirichlet" );
        tParameterList( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        tParameterList( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        tParameterList( 3 ).set( "mesh_set_names", tBottom );

        ////////////////////////////////////////////////////////////////////////////////
        // Volume IQI
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name", "IQITotalVolume" );
        tParameterList( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        tParameterList( 4 ).set( "leader_dof_dependencies", tDofString );
        tParameterList( 4 ).set( "mesh_set_names", tBulk );

        // Eigen values and Eigenvectors

        for ( sint i = 0; i < tNumEigenVectors; ++i )
        {
            // eigen value
            tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            tParameterList( 4 ).set( "IQI_name", "IQIEIGENVALUE" + std::to_string( i + 1 ) );
            tParameterList( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VALUE ) );
            tParameterList( 4 ).set( "leader_dof_dependencies", tDofString );
            tParameterList( 4 ).set( "vectorial_field_index", i );
            tParameterList( 4 ).set( "function_parameters", "0" );
            tParameterList( 4 ).set( "mesh_set_names", tBulk );

            // X-displacement
            tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            tParameterList( 4 ).set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "X" );
            tParameterList( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
            tParameterList( 4 ).set( "dof_quantity", tDofString );
            tParameterList( 4 ).set( "leader_dof_dependencies", tDofString );
            tParameterList( 4 ).set( "vectorial_field_index", 0 );
            tParameterList( 4 ).set( "function_parameters", std::to_string( i ) );
            tParameterList( 4 ).set( "mesh_set_names", tBulk );

            // Y-displacement
            tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            tParameterList( 4 ).set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "Y" );
            tParameterList( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
            tParameterList( 4 ).set( "dof_quantity", tDofString );
            tParameterList( 4 ).set( "leader_dof_dependencies", tDofString );
            tParameterList( 4 ).set( "vectorial_field_index", 1 );
            tParameterList( 4 ).set( "function_parameters", std::to_string( i ) );
            tParameterList( 4 ).set( "mesh_set_names", tBulk );

            if ( tDim == 3 )
            {
                // Z-displacement
                tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
                tParameterList( 4 ).set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "Z" );
                tParameterList( 4 ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
                tParameterList( 4 ).set( "dof_quantity", tDofString );
                tParameterList( 4 ).set( "leader_dof_dependencies", tDofString );
                tParameterList( 4 ).set( "vectorial_field_index", 2 );
                tParameterList( 4 ).set( "function_parameters", std::to_string( i ) );
                tParameterList( 4 ).set( "mesh_set_names", tBulk );
                    }
        }

        // create computation parameter list
        tParameterList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 8 );


        if ( tUsePetsc )
        {
            create_petsc_parameter_list( tParameterlist );
        }
        else
        {
            create_trilinos_parameter_list( tParameterlist );
        }
    }

    void
    create_trilinos_parameter_list( Vector< Submodule_Parameter_Lists >& aParameterlist )
    {
        //---------------------------------------------------------------------------------------------------------------------------------------
        aParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        aParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::EIGEN_SOLVER ) );
        aParameterlist( 0 ).set( "Eigen_Algorithm", "EIGALG_BLOCK_DAVIDSON" );
        aParameterlist( 0 ).set( "Verbosity", false );
        aParameterlist( 0 ).set( "Which", "SM" );
        aParameterlist( 0 ).set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterlist( 0 ).set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterlist( 0 ).set( "Num_Blocks", 2 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterlist( 0 ).set( "MaxSubSpaceDims", 75 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterlist( 0 ).set( "Initial_Guess", 0 );
        aParameterlist( 0 ).set( "MaxRestarts", 20 );
        aParameterlist( 0 ).set( "Convergence_Tolerance", 1e-01 );
        aParameterlist( 0 ).set( "Relative_Convergence_Tolerance", true );
        aParameterlist( 0 ).set( "preconditioners", "0" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterlist( 1 ).set( "DLA_Linear_solver_algorithms", "0" );

        aParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterlist( 1 ).set( "DLA_Linear_solver_algorithms", "1" );
        aParameterlist( 1 ).set( "RHS_Matrix_Type", "GeomStiffMat" );    // MassMat or IdentityMat

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterlist( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterlist( 2 ).set( "NLA_max_iter", 2 );
        aParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", true );

        aParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterlist( 2 ).set( "NLA_Linear_solver", 1 );
        aParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterlist( 2 ).set( "NLA_max_iter", 1 );
        aParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", false );
        aParameterlist( 2 ).set( "NLA_is_eigen_problem", true );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0,1" );
        aParameterlist( 3 ).set( "NLA_DofTypes", tDofString );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterlist( 5 ).set( "TSA_DofTypes", tDofString );
        aParameterlist( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterlist( 6 ).set( "SOL_save_operator_to_matlab", "LBAMat" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        // Ifpack Preconditioner parameters
        aParameterlist( 7 ).set( "ifpack_prec_type", "Amesos" );
        aParameterlist( 7 ).set( "amesos: solver type", "Amesos_Pardiso" );    // Amesos_Umfpack or Amesos_Pardiso

        // Preconditioner parameters
        aParameterlist( 7 ).set( "overlap-level", 0 );
        aParameterlist( 7 ).set( "schwarz: combine mode", "add" );    // for Amesos_Umfpack and Amesos_Pardiso provide this parameter with "add" mode
    }

    void
    create_petsc_parameter_list( Vector< Submodule_Parameter_Lists >& aParameterlist )
    {
        //---------------------------------------------------------------------------------------------------------------------------------------


        aParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        aParameterlist( 0 ).set( "KSPType", "gmres" );
        aParameterlist( 0 ).set( "preconditioners", "0" );

        // find max eigen value
        aParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::SLEPC_SOLVER ) );
        aParameterlist( 0 ).set( "Eigen_Algorithm", "krylovschur" );
        aParameterlist( 0 ).set( "Which", std::string( "LM" ) );
        aParameterlist( 0 ).set( "Num_Eig_Vals", 5 );
        aParameterlist( 0 ).set( "STType", "shift_invert" );
        aParameterlist( 0 ).set( "sub_linear_solver", "0" );    // 10 shift_invert
        aParameterlist( 0 ).set( "is_symmetric", false );       // 10 shift_invert
        aParameterlist( 0 ).set( "Update_Flag", true );         // 10 shift_invert
        aParameterlist( 0 ).set( "Verbosity", false );

        // precondioerr
        aParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC ) );
        aParameterlist( 7 ).set( "PCType", "mumps" );

        //---------------------------------------------------------------------------------------------------------------------------------------


        aParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterlist( 1 ).set( "DLA_Linear_solver_algorithms", "0" );

        aParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterlist( 1 ).set( "DLA_Linear_solver_algorithms", "1" );
        aParameterlist( 1 ).set( "RHS_Matrix_Type", "GeomStiffMat" );    // MassMat or IdentityMat

        //---------------------------------------------------------------------------------------------------------------------------------------


        aParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterlist( 2 ).set( "NLA_Linear_solver", 0 );
        aParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterlist( 2 ).set( "NLA_max_iter", 2 );
        aParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", true );

        aParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterlist( 2 ).set( "NLA_Linear_solver", 1 );
        aParameterlist( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterlist( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterlist( 2 ).set( "NLA_max_iter", 1 );
        aParameterlist( 2 ).set( "NLA_combined_res_jac_assembly", false );
        aParameterlist( 2 ).set( "NLA_is_eigen_problem", true );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterlist( 3 ).set( "NLA_Nonlinear_solver_algorithms", "0,1" );
        aParameterlist( 3 ).set( "NLA_DofTypes", tDofString );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterlist( 5 ).set( "TSA_DofTypes", tDofString );
        aParameterlist( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterlist( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterlist( 6 ).set( "SOL_TPL_Type", sol::MapType::Petsc );
    }

    void MSIParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        tParameterlist( 0 ).set( "number_eigen_vectors", tNumEigenVectors );
    }

    void VISParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 ).set( "Set_Names", tBulk );

        std::string Field_Names = "VOLUME";
        std::string Field_Type  = "GLOBAL";
        std::string IQI_Names   = "IQITotalVolume";

        for ( sint i = 0; i < tNumEigenVectors; ++i )
        {
            Field_Names += ",";
            Field_Type += ",";
            IQI_Names += ",";

            Field_Names += ",EVAL" + std::to_string( i + 1 ) + ",EVEC" + std::to_string( i + 1 ) + "X,EVEC" + std::to_string( i + 1 ) + "Y";
            Field_Type += ",GLOBAL,NODAL,NODAL";
            IQI_Names += ",IQIEIGENVALUE" + std::to_string( i + 1 ) + ",IQIEIGENVEC" + std::to_string( i + 1 ) + "X,IQIEIGENVEC" + std::to_string( i + 1 ) + "Y";
            if ( tDim == 3 )
            {
                Field_Names += ",EVEC" + std::to_string( i + 1 ) + "Z";
                Field_Type += ",NODAL";
                IQI_Names += ",IQIEIGENVEC" + std::to_string( i + 1 ) + "Z";
            }
        }

        tParameterlist( 0 ).set( "Field_Names", Field_Names );
        tParameterlist( 0 ).set( "Field_Type", Field_Type );
        tParameterlist( 0 ).set( "IQI_Names", IQI_Names );

        tParameterlist( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Submodule_Parameter_Lists >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
