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
#include "parameters.hpp"
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
    void create_trilinos_parameter_list( Module_Parameter_Lists& );
    void create_petsc_parameter_list( Module_Parameter_Lists& );

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

    Vector< uint > tNumElemsPerDim = tDim == 3 ? Vector< uint >{ 5, 50, 5 } : Vector< uint >{ 5, 50 };
    Vector< real > tDomainDims     = tDim == 3 ? Vector< real >{ 1.0, 10.0, 1.0 } : Vector< real >{ 1.0, 10.0 };

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

    void OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "lagrange_output_meshes", std::string( "0" ) );

        aParameterLists.set( "lagrange_orders", tOrder );
        aParameterLists.set( "lagrange_pattern", std::string( "0" ) );
        aParameterLists.set( "bspline_orders", tOrder );
        aParameterLists.set( "bspline_pattern", std::string( "0" ) );


        aParameterLists.set( "initial_refinement", "0" );


        aParameterLists.set( "lagrange_mesh_output_file_name", "LBAProblem.exo" );
    }

    /* ------------------------------------------------------------------------ */

    void XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", std::string( "conformal" ) );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", std::string( "bspline" ) );
        aParameterLists.set( "enrich_mesh_indices", std::string( "0" ) );
        aParameterLists.set( "ghost_stab", false );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */

    void FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Youngs Modulus Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungsModulus" );
        aParameterLists.set( "function_parameters", tYoungsModulus );

        // Poisson Ratio Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPoissonRatio" ) );
        aParameterLists.set( "function_parameters", tPoissonRatio );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // Compressive load
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropTraction" ) );
        aParameterLists.set( "function_parameters", tTraction );

        // Dirichlet structure
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropDirichletStruct" ) );
        aParameterLists.set( "function_parameters", tDirichlet );

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        // linear elasticity - shell - 1
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso" );
        aParameterLists.set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
        aParameterLists.set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofString, "Displacement" ) );
        aParameterLists.set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucNonLinIso" );
        aParameterLists.set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
        aParameterLists.set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF ) );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( tDofString, "Displacement" ) );
        aParameterLists.set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );

        ////////////////////////////////////////////////////////////////////////////////

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        // Displacements - Shell - back wall
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheStruc" );
        aParameterLists.set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungsModulus,Material" );

        ////////////////////////////////////////////////////////////////////////////////
        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // linear elasticity
        ///*
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGStructBulk" );
        aParameterLists.set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
        aParameterLists.set( "dof_residual", tDofString );
        aParameterLists.set( "leader_dof_dependencies", tDofString );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tBulk );
        //*/

        //------------------------------------------------------------------------------

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGStructBulkLBA" );
        aParameterLists.set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_NON_LINEAR_GEOMETRIC_STIFFNESS ) );
        aParameterLists.set( "dof_residual", tDofString );
        aParameterLists.set( "leader_dof_dependencies", tDofString );
        if ( tStressType == "PK2" )
        {
            aParameterLists.set( "leader_constitutive_models", "CMStrucNonLinIso,ElastLinIso" );
        }
        else
        {
            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        }
        aParameterLists.set( "mesh_set_names", tBulk );

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // Traction BC ( Compressive force at top )
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGNeumannTraction" );
        aParameterLists.set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
        aParameterLists.set( "dof_residual", tDofString );
        aParameterLists.set( "leader_dof_dependencies", tDofString );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tTop );

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------
        // Fixed bottom edge
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletStruct" );
        aParameterLists.set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
        aParameterLists.set( "dof_residual", tDofString );
        aParameterLists.set( "leader_dof_dependencies", tDofString );
        aParameterLists.set( "leader_properties", "PropDirichletStruct,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tBottom );

        ////////////////////////////////////////////////////////////////////////////////
        // Volume IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQITotalVolume" );
        aParameterLists.set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );
        aParameterLists.set( "leader_dof_dependencies", tDofString );
        aParameterLists.set( "mesh_set_names", tBulk );

        // Eigen values and Eigenvectors

        for ( sint i = 0; i < tNumEigenVectors; ++i )
        {
            // eigen value
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIEIGENVALUE" + std::to_string( i + 1 ) );
            aParameterLists.set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VALUE ) );
            aParameterLists.set( "leader_dof_dependencies", tDofString );
            aParameterLists.set( "vectorial_field_index", i );
            aParameterLists.set( "function_parameters", "0" );
            aParameterLists.set( "mesh_set_names", tBulk );

            // X-displacement
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "X" );
            aParameterLists.set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
            aParameterLists.set( "dof_quantity", tDofString );
            aParameterLists.set( "leader_dof_dependencies", tDofString );
            aParameterLists.set( "vectorial_field_index", 0 );
            aParameterLists.set( "function_parameters", std::to_string( i ) );
            aParameterLists.set( "mesh_set_names", tBulk );

            // Y-displacement
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "Y" );
            aParameterLists.set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
            aParameterLists.set( "dof_quantity", tDofString );
            aParameterLists.set( "leader_dof_dependencies", tDofString );
            aParameterLists.set( "vectorial_field_index", 1 );
            aParameterLists.set( "function_parameters", std::to_string( i ) );
            aParameterLists.set( "mesh_set_names", tBulk );

            if ( tDim == 3 )
            {
                // Z-displacement
                aParameterLists( FEM::IQI ).add_parameter_list();
                aParameterLists.set( "IQI_name", "IQIEIGENVEC" + std::to_string( i + 1 ) + "Z" );
                aParameterLists.set( "IQI_type", static_cast< uint >( fem::IQI_Type::EIGEN_VECTOR ) );
                aParameterLists.set( "dof_quantity", tDofString );
                aParameterLists.set( "leader_dof_dependencies", tDofString );
                aParameterLists.set( "vectorial_field_index", 2 );
                aParameterLists.set( "function_parameters", std::to_string( i ) );
                aParameterLists.set( "mesh_set_names", tBulk );
                    }
        }

        // create computation parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        if ( tUsePetsc )
        {
            create_petsc_parameter_list( aParameterLists );
        }
        else
        {
            create_trilinos_parameter_list( aParameterLists );
        }
    }

    void
    create_trilinos_parameter_list( Module_Parameter_Lists& aParameterLists )
    {
        //---------------------------------------------------------------------------------------------------------------------------------------
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::EIGEN_SOLVER );
        aParameterLists.set( "Eigen_Algorithm", "EIGALG_BLOCK_DAVIDSON" );
        aParameterLists.set( "Verbosity", false );
        aParameterLists.set( "Which", "SM" );
        aParameterLists.set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterLists.set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterLists.set( "Num_Blocks", 2 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterLists.set( "MaxSubSpaceDims", 75 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterLists.set( "Initial_Guess", 0 );
        aParameterLists.set( "MaxRestarts", 20 );
        aParameterLists.set( "Convergence_Tolerance", 1e-01 );
        aParameterLists.set( "Relative_Convergence_Tolerance", true );
        aParameterLists.set( "preconditioners", "0" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "1" );
        aParameterLists.set( "RHS_Matrix_Type", "GeomStiffMat" );    // MassMat or IdentityMat

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", 2 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 1 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", 1 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );
        aParameterLists.set( "NLA_is_eigen_problem", true );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0,1" );
        aParameterLists.set( "NLA_DofTypes", tDofString );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", tDofString );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_save_operator_to_matlab", "LBAMat" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::IFPACK );
        // Ifpack Preconditioner parameters
        aParameterLists.set( "ifpack_prec_type", "Amesos" );
        aParameterLists.set( "amesos: solver type", "Amesos_Pardiso" );    // Amesos_Umfpack or Amesos_Pardiso

        // Preconditioner parameters
        aParameterLists.set( "overlap-level", 0 );
        aParameterLists.set( "schwarz: combine mode", "add" );    // for Amesos_Umfpack and Amesos_Pardiso provide this parameter with "add" mode
    }

    void
    create_petsc_parameter_list( Module_Parameter_Lists& aParameterLists )
    {
        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::PETSC );
        aParameterLists.set( "KSPType", "gmres" );
        aParameterLists.set( "preconditioners", "0" );

        // find max eigen value
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::SLEPC_SOLVER );
        aParameterLists.set( "Eigen_Algorithm", "krylovschur" );
        aParameterLists.set( "Which", std::string( "LM" ) );
        aParameterLists.set( "Num_Eig_Vals", 5 );
        aParameterLists.set( "STType", "shift_invert" );
        aParameterLists.set( "sub_linear_solver", "0" );    // 10 shift_invert
        aParameterLists.set( "is_symmetric", false );       // 10 shift_invert
        aParameterLists.set( "Update_Flag", true );         // 10 shift_invert
        aParameterLists.set( "Verbosity", false );

        // precondioerr
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::PETSC );
        aParameterLists.set( "PCType", "mumps" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "1" );
        aParameterLists.set( "RHS_Matrix_Type", "GeomStiffMat" );    // MassMat or IdentityMat

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", 2 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", true );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 1 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", 1 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );
        aParameterLists.set( "NLA_is_eigen_problem", true );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0,1" );
        aParameterLists.set( "NLA_DofTypes", tDofString );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", tDofString );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        //---------------------------------------------------------------------------------------------------------------------------------------

        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_TPL_Type", sol::MapType::Petsc );
    }

    void MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_eigen_vectors", tNumEigenVectors );
    }

    void VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        aParameterLists.set( "Set_Names", tBulk );

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

        aParameterLists.set( "Field_Names", Field_Names );
        aParameterLists.set( "Field_Type", Field_Type );
        aParameterLists.set( "IQI_Names", IQI_Names );

        aParameterLists.set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
