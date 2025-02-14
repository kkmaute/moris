/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Cantilever_Eigen.cpp
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
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"
#include "fn_stringify_matrix.hpp"

#include "AztecOO.h"

//--------------------------------------------------------------------------------

// global variable - interpolation order
extern uint gOrder;

// global variable - Preconditioner solver
extern std::string gPrecSolver;

// global variable - TestCaseIndex
extern uint gTestCaseIndex;

//--------------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    void create_trilinos_parameter_list( Module_Parameter_Lists& ), create_petsc_parameter_list( Module_Parameter_Lists& );
    /* ------------------------------------------------------------------------ */
    // General

    // file name
    std::string tName = "Cantilever_Eigen";

    // activate structural portion
    bool tHaveStruct     = true;
    bool tTimeContinuity = true;

    // Traction load
    std::string tTraction = "100.0";

    /* ------------------------------------------------------------------------ */
    // Solver config

    int         tNLA_max_iter             = 1;
    moris::real tNLA_rel_res_norm_drop    = 1.0e-10;
    moris::real tNLA_relaxation_parameter = 1.0;

    int         tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame     = 1;

    /* ------------------------------------------------------------------------ */
    // HMR parameters

    // Interpolation order
    std::string tOrder = "1";

    Vector< uint > tNumElemsPerDim = { 40, 1200 };
    Vector< real > tDomainDims     = { 0.01, 0.30 };

    /* ------------------------------------------------------------------------ */
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    // Bulk Phases
    std::string tBulk = "HMR_dummy_n_p0";

    // boundaries
    std::string tBottom = "SideSet_1_n_p0";
    std::string tRight  = "SideSet_2_n_p0";
    std::string tTop    = "SideSet_3_n_p0";
    std::string tLeft   = "SideSet_4_n_p0";

    /* ------------------------------------------------------------------------ */
    // material parameters, kg is scaled with a factor 1e-6

    // material properties
    std::string tDensity       = "1.0";
    std::string tYoungsModulus = "1";
    std::string tPoissonRatio  = "0.0";

    /* ------------------------------------------------------------------------ */
    // boundary conditions

    // bedding to supress RBM
    std::string tBedding = std::to_string( 1.0 * 1.0e-5 );

    moris::real tInitialStruct = 0.0;

    /* ------------------------------------------------------------------------ */
    // Output Config
    std::string tOutputFileName = tName + "_"+ gPrecSolver + ".exo" ;

    std::string tLibraryName     = tName + ".so";
    std::string tHDF5Path        = tName + ".hdf5";
    std::string tGENOutputFile   = "GEN_" + tName + ".exo";
    bool        tOutputCriterion = true;

    //------------------------------------------------------------------------------
    //-------------------------------- FUNCTIONS -----------------------------------
    //------------------------------------------------------------------------------

    /* ------------------------------------------------------------------------ */
    // GEOMETRY (LEVEL-SET) FUNCTIONS
    /* ------------------------------------------------------------------------ */

    //-----------------------------------------------------------------------------

    // Back Wall
    moris::real
    Back_Wall(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const Vector< real >& aGeometryParameters )
    {
        real tVal = aCoordinates( 0 ) - 100.0;

        // clean return value to return non-zero value
        return std::abs( tVal ) < 1.0e-8 ? 1.0e-8 : tVal;
    }

    /* ------------------------------------------------------------------------ */
    // PROPERTY FUNCTIONS (incl. INITIAL & BOUNDARY CONDITIONS)
    /* ------------------------------------------------------------------------ */

    // Constant function for properties
    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // initial structure
    void
    Func_Initial_Condition_Struct(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        aPropMatrix = { { tInitialStruct }, { tInitialStruct } };
    }

    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&                aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*        aFIManager )
    {
        // aPropMatrix = { { aParameters( 0 )( 0 ) },{ 0.0 } };
        aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) } };
    }

    /* ------------------------------------------------------------------------ */
    // DUMMY FUNCTIONS
    /* ------------------------------------------------------------------------ */

    // Output criterion for VIS mesh
    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return tOutputCriterion;
    }

    // Dummy function for unused sensitivities if needed
    moris::Matrix< DDRMat >
    Func_Dummy_Sensitivity(
            const moris::Matrix< DDRMat >&     aCoordinates,
            const Vector< real >& aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue = { { 0.0 } };
        return aReturnValue;
    }

    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists.set( "domain_dimensions", tDomainDims );
        aParameterLists.set( "lagrange_output_meshes", std::string( "0" ) );

        aParameterLists.set( "lagrange_orders", tOrder );
        aParameterLists.set( "lagrange_pattern", std::string( "0" ) );
        aParameterLists.set( "bspline_orders", tOrder );
        aParameterLists.set( "bspline_pattern", std::string( "0" ) );


        aParameterLists.set( "initial_refinement", "0" );


        aParameterLists.set( "lagrange_mesh_output_file_name", "Cantilever_Eigen_HMR.exo" );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", std::string( "conformal" ) );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", std::string( "bspline" ) );
        aParameterLists.set( "enrich_mesh_indices", std::string( "0" ) );
        aParameterLists.set( "ghost_stab", false );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "output_mesh_file", tGENOutputFile );

        // Dummy Geometry
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Back_Wall" );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Density Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", tDensity );
        aParameterLists.set( "value_function", "Func_Const" );

        // Youngs Modulus Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungsModulus" );
        aParameterLists.set( "function_parameters", tYoungsModulus );
        aParameterLists.set( "value_function", "Func_Const" );

        // Poisson Ratio Shell
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPoissonRatio" ) );
        aParameterLists.set( "function_parameters", tPoissonRatio );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // pressure load
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropTraction" ) );
        aParameterLists.set( "function_parameters", tTraction );
        aParameterLists.set( "value_function", std::string( "Func_Neumann_U" ) );

        // pressure load
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropPressure" ) );
        aParameterLists.set( "function_parameters", tTraction );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // Dirichlet structure
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropDirichletStruct" ) );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // time continuity weights
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropWeightCurrent" ) );
        aParameterLists.set( "function_parameters", std::string( tDensity ) );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropWeightPrevious" ) );
        aParameterLists.set( "function_parameters", std::string( tDensity ) );
        aParameterLists.set( "value_function", std::string( "Func_Const" ) );

        // Initial Structure
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", std::string( "PropInitialConditionStruct" ) );
        aParameterLists.set( "value_function", std::string( "Func_Initial_Condition_Struct" ) );

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // linear elasticity
            aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
            aParameterLists.set( "constitutive_name", "CMStrucLinIso" );
            aParameterLists.set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
            aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
            aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
            aParameterLists.set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );
            }

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // Displacements - Shell - back wall
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPNitscheStruc" );
            aParameterLists.set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
            aParameterLists.set( "function_parameters", "100.0" );
            aParameterLists.set( "leader_properties", "PropYoungsModulus,Material" );
            }

        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // linear elasticity
        if ( tHaveStruct )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGStructShell" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
            aParameterLists.set( "dof_residual", "UX,UY" );

            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            // aParameterLists.set( "leader_properties",          "PropBedding,Bedding" );
            aParameterLists.set( "mesh_set_names", tBulk );
            }

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // pressure pulling on outside of Shell
        if ( tHaveStruct )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGNeumannPressure" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
            aParameterLists.set( "dof_residual", "UX,UY" );

            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists.set( "leader_properties", "PropTraction,Traction" );
            aParameterLists.set( "mesh_set_names", tLeft );
            }

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - shell - back wall
        if ( tHaveStruct )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGDirichletStruct" );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) ;
            aParameterLists.set( "dof_residual", "UX,UY" );

            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists.set( "leader_properties", "PropDirichletStruct,Dirichlet" );
            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            aParameterLists.set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
            aParameterLists.set( "mesh_set_names", tBottom );
            }

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity Structure
        if ( tHaveStruct )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", std::string( "IWGTimeContinuityStruct" ) );
            aParameterLists.set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
            aParameterLists.set( "dof_residual", std::string( "UX,UY" ) );

            {
                aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists.set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialConditionStruct,InitialCondition" ) );
            aParameterLists.set( "mesh_set_names", tBulk );
            aParameterLists.set( "time_continuity", tTimeContinuity );
            }

        //Volume IQI - TotalDomain - use once to find total volume to compute max dof
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQITotalVolume" );
        aParameterLists.set( "IQI_type",  fem::IQI_Type::VOLUME ) ;

        if ( tHaveStruct )
        {
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        }
        aParameterLists.set( "mesh_set_names", tBulk );

        if ( tHaveStruct )
        {
            // X-displacement
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkDISPX" );
            aParameterLists.set( "IQI_type", ( fem::IQI_Type::EIGEN_VECTOR ) );
            aParameterLists.set( "function_parameters", "1" );
            aParameterLists.set( "dof_quantity", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "vectorial_field_index", 0 );
            aParameterLists.set( "mesh_set_names", tBulk );

            // Y-displacement
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "IQIBulkDISPY" );
            aParameterLists.set( "IQI_type", ( fem::IQI_Type::EIGEN_VECTOR ) );
            aParameterLists.set( "function_parameters", "1" );
            aParameterLists.set( "dof_quantity", "UX,UY" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists.set( "vectorial_field_index", 1 );
            aParameterLists.set( "mesh_set_names", tBulk );
            }

        // create computation parameter list
        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        gPrecSolver == "Slepc" ? create_petsc_parameter_list( aParameterLists ) : create_trilinos_parameter_list( aParameterLists );

    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_eigen_vectors", 5 );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names", tBulk );

        if ( tHaveStruct )
        {
            aParameterLists.set( "Field_Names", "UX,UY" );
            aParameterLists.set( "Field_Type", "NODAL,NODAL" );
            aParameterLists.set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY" );
        }

        aParameterLists.set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
    void
    create_trilinos_parameter_list( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::EIGEN_SOLVER );
        aParameterLists.set( "Eigen_Algorithm", "EIGALG_BLOCK_DAVIDSON" );
        aParameterLists.set( "Verbosity", false );
        aParameterLists.set( "Which", "SM" );
        aParameterLists.set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterLists.set( "NumFreeDofs", 1000 );      // For 2D problem of rectangular elements number of free dofs = 2*node_x*node_y
        aParameterLists.set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterLists.set( "Num_Blocks", 2 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterLists.set( "MaxSubSpaceDims", 75 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterLists.set( "Initial_Guess", 0 );
        aParameterLists.set( "MaxRestarts", 20 );
        aParameterLists.set( "Convergence_Tolerance", 1e-05 );
        aParameterLists.set( "Relative_Convergence_Tolerance", true );
        aParameterLists.set( "preconditioners", "0" );

        // Ml Preconditioner parameters
        // aParameterLists.set( "ml_prec_type", "NSSA"); // options: SA, NSSA, DD, DD-ML

        // Print eigenvector parameter
        //aParameterLists.set( "Print_vector", "LINSOL_EXPORT_MATLAB" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );
        aParameterLists.set( "RHS_Matrix_Type", "MassMat" );    // MassMat or IdentityMat

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "UX,1;UY,1" );

        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_save_operator_to_matlab", "MassMat" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::IFPACK );
        // Ifpack Preconditioner parameters
        aParameterLists.set( "ifpack_prec_type", "Amesos" );
        aParameterLists.set( "amesos: solver type", gPrecSolver );    // Amesos_Umfpack or Amesos_Pardiso
        // Preconditioner parameters
        aParameterLists.set( "overlap-level", 0 );
        aParameterLists.set( "schwarz: combine mode", "add" );    // for Amesos_Umfpack and Amesos_Pardiso provide this parameter with "add" mode
    }

    //---------------------------------------------------------------------------

    void
    create_petsc_parameter_list( Module_Parameter_Lists& aParameterLists )
    {

        // 1 slpec solver and the associated linear solver object
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::PETSC );
        aParameterLists.set( "KSPType", "gmres" );
        aParameterLists.set( "preconditioners", "0" );    // 10 shift_invert

        // find max eigen value
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::SLEPC_SOLVER );
        aParameterLists.set( "Eigen_Algorithm", "krylovschur" );
        aParameterLists.set( "Which", std::string( "LM" ) );
        aParameterLists.set( "Num_Eig_Vals", 5 );
        aParameterLists.set( "STType", "shift_invert" );
        aParameterLists.set( "sub_linear_solver", "0" );    // 10 shift_invert
        aParameterLists.set( "is_symmetric", true );       // 10 shift_invert
        aParameterLists.set( "Update_Flag", true );         // 10 shift_invert
        aParameterLists.set( "Verbosity", false );

        // precondioerr
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::PETSC );
        aParameterLists.set( "PCType", "mumps" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "1" );
        aParameterLists.set( "RHS_Matrix_Type", "MassMat" );    // MassMat or IdentityMat

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists.set( "TSA_Time_Frame", tTSA_Time_Frame );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "UX,1;UY,1" );

        // aParameterLists.set( "SOL_save_operator_to_matlab", "MassMat" );
        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_TPL_Type",  sol::MapType::Petsc ) ;
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
