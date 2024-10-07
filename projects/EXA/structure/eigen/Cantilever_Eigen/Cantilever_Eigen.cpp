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

    std::string tNumElemsPerDim = "40,1200";
    std::string tDomainDims     = "0.01, 0.30";
    std::string tDomainOffset   = "0.0,0.0";
    std::string tDomainSidesets = "1,2,3,4";

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
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_opt_problem_parameter_list() );
        aParameterLists( 0 ).set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        aParameterLists( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        aParameterLists( 0 ).set( "domain_dimensions", tDomainDims );
        aParameterLists( 0 ).set( "domain_offset", tDomainOffset );
        aParameterLists( 0 ).set( "domain_sidesets", tDomainSidesets );
        aParameterLists( 0 ).set( "lagrange_output_meshes", std::string( "0" ) );

        aParameterLists( 0 ).set( "lagrange_orders", tOrder );
        aParameterLists( 0 ).set( "lagrange_pattern", std::string( "0" ) );
        aParameterLists( 0 ).set( "bspline_orders", tOrder );
        aParameterLists( 0 ).set( "bspline_pattern", std::string( "0" ) );

        aParameterLists( 0 ).set( "truncate_bsplines", 1 );

        aParameterLists( 0 ).set( "use_number_aura", 1 );

        aParameterLists( 0 ).set( "initial_refinement", "0" );

        aParameterLists( 0 ).set( "initial_refinement_pattern", "0" );

        aParameterLists( 0 ).set( "use_multigrid", 0 );
        aParameterLists( 0 ).set( "severity_level", 0 );

        aParameterLists( 0 ).set( "write_lagrange_output_mesh_to_exodus", "Cantilever_Eigen_HMR.exo" );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        aParameterLists( 0 ).set( "decompose", true );
        aParameterLists( 0 ).set( "decomposition_type", std::string( "conformal" ) );
        aParameterLists( 0 ).set( "enrich", true );
        aParameterLists( 0 ).set( "basis_rank", std::string( "bspline" ) );
        aParameterLists( 0 ).set( "enrich_mesh_indices", std::string( "0" ) );
        aParameterLists( 0 ).set( "ghost_stab", false );
        aParameterLists( 0 ).set( "multigrid", false );
        aParameterLists( 0 ).set( "verbose", true );
        aParameterLists( 0 ).set( "print_enriched_ig_mesh", true );
        aParameterLists( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
        aParameterLists( 0 ).set( "output_mesh_file", tGENOutputFile );

        // Dummy Geometry
        aParameterLists( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists( 1 ).set( "field_function_name", "Back_Wall" );
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
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropDensity" );
        aParameterLists( 0 ).set( "function_parameters", tDensity );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Youngs Modulus Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", "PropYoungsModulus" );
        aParameterLists( 0 ).set( "function_parameters", tYoungsModulus );
        aParameterLists( 0 ).set( "value_function", "Func_Const" );

        // Poisson Ratio Shell
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropPoissonRatio" ) );
        aParameterLists( 0 ).set( "function_parameters", tPoissonRatio );
        aParameterLists( 0 ).set( "value_function", std::string( "Func_Const" ) );

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // pressure load
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropTraction" ) );
        aParameterLists( 0 ).set( "function_parameters", tTraction );
        aParameterLists( 0 ).set( "value_function", std::string( "Func_Neumann_U" ) );

        // pressure load
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropPressure" ) );
        aParameterLists( 0 ).set( "function_parameters", tTraction );
        aParameterLists( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // Dirichlet structure
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropDirichletStruct" ) );
        aParameterLists( 0 ).set( "function_parameters", "0.0;0.0" );
        aParameterLists( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // time continuity weights
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropWeightCurrent" ) );
        aParameterLists( 0 ).set( "function_parameters", std::string( tDensity ) );
        aParameterLists( 0 ).set( "value_function", std::string( "Func_Const" ) );

        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropWeightPrevious" ) );
        aParameterLists( 0 ).set( "function_parameters", std::string( tDensity ) );
        aParameterLists( 0 ).set( "value_function", std::string( "Func_Const" ) );

        // Initial Structure
        aParameterLists( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        aParameterLists( 0 ).set( "property_name", std::string( "PropInitialConditionStruct" ) );
        aParameterLists( 0 ).set( "value_function", std::string( "Func_Initial_Condition_Struct" ) );

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // linear elasticity
            aParameterLists( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
            aParameterLists( 1 ).set( "constitutive_name", "CMStrucLinIso" );
            aParameterLists( 1 ).set( "model_type",  fem::Model_Type::PLANE_STRESS ) ;
            aParameterLists( 1 ).set( "constitutive_type",  fem::Constitutive_Type::STRUC_LIN_ISO ) ;
            aParameterLists( 1 ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
            aParameterLists( 1 ).set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );
            }

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // Displacements - Shell - back wall
            aParameterLists( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            aParameterLists( 2 ).set( "stabilization_name", "SPNitscheStruc" );
            aParameterLists( 2 ).set( "stabilization_type",  fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
            aParameterLists( 2 ).set( "function_parameters", "100.0" );
            aParameterLists( 2 ).set( "leader_properties", "PropYoungsModulus,Material" );
            }

        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // linear elasticity
        if ( tHaveStruct )
        {
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGStructShell" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_BULK ) ;
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );

            {
                aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            // aParameterLists( 3 ).set( "leader_properties",          "PropBedding,Bedding" );
            aParameterLists( 3 ).set( "mesh_set_names", tBulk );
            }

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // pressure pulling on outside of Shell
        if ( tHaveStruct )
        {
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGNeumannPressure" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_NEUMANN ) ;
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );

            {
                aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists( 3 ).set( "leader_properties", "PropTraction,Traction" );
            aParameterLists( 3 ).set( "mesh_set_names", tLeft );
            }

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - shell - back wall
        if ( tHaveStruct )
        {
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", "IWGDirichletStruct" );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) ;
            aParameterLists( 3 ).set( "dof_residual", "UX,UY" );

            {
                aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists( 3 ).set( "leader_properties", "PropDirichletStruct,Dirichlet" );
            aParameterLists( 3 ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            aParameterLists( 3 ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
            aParameterLists( 3 ).set( "mesh_set_names", tBottom );
            }

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity Structure
        if ( tHaveStruct )
        {
            aParameterLists( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            aParameterLists( 3 ).set( "IWG_name", std::string( "IWGTimeContinuityStruct" ) );
            aParameterLists( 3 ).set( "IWG_type",  fem::IWG_Type::TIME_CONTINUITY_DOF ) ;
            aParameterLists( 3 ).set( "dof_residual", std::string( "UX,UY" ) );

            {
                aParameterLists( 3 ).set( "leader_dof_dependencies", "UX,UY" );
            }

            aParameterLists( 3 ).set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialConditionStruct,InitialCondition" ) );
            aParameterLists( 3 ).set( "mesh_set_names", tBulk );
            aParameterLists( 3 ).set( "time_continuity", tTimeContinuity );
            }

        //Volume IQI - TotalDomain - use once to find total volume to compute max dof
        aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        aParameterLists( 4 ).set( "IQI_name", "IQITotalVolume" );
        aParameterLists( 4 ).set( "IQI_type",  fem::IQI_Type::VOLUME ) ;

        if ( tHaveStruct )
        {
            aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
        }
        aParameterLists( 4 ).set( "mesh_set_names", tBulk );

        if ( tHaveStruct )
        {
            // X-displacement
            aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPX" );
            aParameterLists( 4 ).set( "IQI_type", ( fem::IQI_Type::EIGEN_VECTOR ) );
            aParameterLists( 4 ).set( "function_parameters", "1" );
            aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
            aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 4 ).set( "vectorial_field_index", 0 );
            aParameterLists( 4 ).set( "mesh_set_names", tBulk );

            // Y-displacement
            aParameterLists( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
            aParameterLists( 4 ).set( "IQI_name", "IQIBulkDISPY" );
            aParameterLists( 4 ).set( "IQI_type", ( fem::IQI_Type::EIGEN_VECTOR ) );
            aParameterLists( 4 ).set( "function_parameters", "1" );
            aParameterLists( 4 ).set( "dof_quantity", "UX,UY" );
            aParameterLists( 4 ).set( "leader_dof_dependencies", "UX,UY" );
            aParameterLists( 4 ).set( "vectorial_field_index", 1 );
            aParameterLists( 4 ).set( "mesh_set_names", tBulk );
            }

        // create computation parameter list
        aParameterLists( 5 ).add_parameter_list( prm::create_computation_parameter_list() );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {

        gPrecSolver == "Slepc" ? create_petsc_parameter_list( aParameterLists ) : create_trilinos_parameter_list( aParameterLists );

    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        aParameterLists( 0 ).set( "number_eigen_vectors", 5 );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists( 0 ).set( "Mesh_Type",  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists( 0 ).set( "Set_Names", tBulk );

        if ( tHaveStruct )
        {
            aParameterLists( 0 ).set( "Field_Names", "UX,UY" );
            aParameterLists( 0 ).set( "Field_Type", "NODAL,NODAL" );
            aParameterLists( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY" );
        }

        aParameterLists( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
    void
    create_trilinos_parameter_list( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::EIGEN_SOLVER ) );
        aParameterLists( 0 ).set( "Eigen_Algorithm", "EIGALG_BLOCK_DAVIDSON" );
        aParameterLists( 0 ).set( "Verbosity", false );
        aParameterLists( 0 ).set( "Which", "SM" );
        aParameterLists( 0 ).set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        aParameterLists( 0 ).set( "NumFreeDofs", 1000 );      // For 2D problem of rectangular elements number of free dofs = 2*node_x*node_y
        aParameterLists( 0 ).set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        aParameterLists( 0 ).set( "Num_Blocks", 2 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        aParameterLists( 0 ).set( "MaxSubSpaceDims", 75 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        aParameterLists( 0 ).set( "Initial_Guess", 0 );
        aParameterLists( 0 ).set( "MaxRestarts", 20 );
        aParameterLists( 0 ).set( "Convergence_Tolerance", 1e-05 );
        aParameterLists( 0 ).set( "Relative_Convergence_Tolerance", true );
        aParameterLists( 0 ).set( "preconditioners", "0" );

        // Ml Preconditioner parameters
        // aParameterLists( 0 ).set( "ml_prec_type", "NSSA"); // options: SA, NSSA, DD, DD-ML

        // Print eigenvector parameter
        //aParameterLists( 0 ).set( "Print_vector", "LINSOL_EXPORT_MATLAB" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "0" );
        aParameterLists( 1 ).set( "RHS_Matrix_Type", "MassMat" );    // MassMat or IdentityMat

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists( 5 ).set( "TSA_time_level_per_type", "UX,1;UY,1" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        aParameterLists( 6 ).set( "SOL_save_operator_to_matlab", "MassMat" );

        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
        // Ifpack Preconditioner parameters
        aParameterLists( 7 ).set( "ifpack_prec_type", "Amesos" );
        aParameterLists( 7 ).set( "amesos: solver type", gPrecSolver );    // Amesos_Umfpack or Amesos_Pardiso
        // Preconditioner parameters
        aParameterLists( 7 ).set( "overlap-level", 0 );
        aParameterLists( 7 ).set( "schwarz: combine mode", "add" );    // for Amesos_Umfpack and Amesos_Pardiso provide this parameter with "add" mode
    }

    //---------------------------------------------------------------------------

    void
    create_petsc_parameter_list( Module_Parameter_Lists& aParameterLists )
    {

        // 1 slpec solver and the associated linear solver object
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::PETSC ) );
        aParameterLists( 0 ).set( "KSPType", "gmres" );
        aParameterLists( 0 ).set( "preconditioners", "0" );    // 10 shift_invert

        // find max eigen value
        aParameterLists( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::SLEPC_SOLVER ) );
        aParameterLists( 0 ).set( "Eigen_Algorithm", "krylovschur" );
        aParameterLists( 0 ).set( "Which", std::string( "LM" ) );
        aParameterLists( 0 ).set( "Num_Eig_Vals", 5 );
        aParameterLists( 0 ).set( "STType", "shift_invert" );
        aParameterLists( 0 ).set( "sub_linear_solver", "0" );    // 10 shift_invert
        aParameterLists( 0 ).set( "is_symmetric", true );       // 10 shift_invert
        aParameterLists( 0 ).set( "Update_Flag", true );         // 10 shift_invert
        aParameterLists( 0 ).set( "Verbosity", false );

        // precondioerr
        aParameterLists( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC ) );
        aParameterLists( 7 ).set( "PCType", "mumps" );

        aParameterLists( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );
        aParameterLists( 1 ).set( "DLA_Linear_solver_algorithms", "1" );
        aParameterLists( 1 ).set( "RHS_Matrix_Type", "MassMat" );    // MassMat or IdentityMat

        aParameterLists( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        aParameterLists( 2 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists( 2 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists( 2 ).set( "NLA_max_iter", tNLA_max_iter );
        aParameterLists( 2 ).set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        aParameterLists( 3 ).set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        aParameterLists( 4 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists( 4 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

        aParameterLists( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        aParameterLists( 5 ).set( "TSA_DofTypes", "UX,UY" );
        aParameterLists( 5 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        aParameterLists( 5 ).set( "TSA_Output_Indices", "0" );
        aParameterLists( 5 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists( 5 ).set( "TSA_time_level_per_type", "UX,1;UY,1" );

        aParameterLists( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        // aParameterLists( 6 ).set( "SOL_save_operator_to_matlab", "MassMat" );
        aParameterLists( 6 ).set( "SOL_TPL_Type",  sol::MapType::Petsc ) ;

    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
