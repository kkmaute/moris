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
    std::string tOutputFileName = ( gPrecSolver == "Amesos_Pardiso" ) ? tName + "_Pardiso" + ".exo" : tName + "_Umfpack" + ".exo";

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
            const Vector< moris::real* >& aGeometryParameters )
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
            const Vector< moris::real* >& aGeometryParameters )
    {
        moris::Matrix< DDRMat > aReturnValue = { { 0.0 } };
        return aReturnValue;
    }


    /* ------------------------------------------------------------------------ */
    // PARAMETER LISTS
    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );
        tParameterlist( 1 ).resize( 0 );
        tParameterlist( 2 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = moris::prm::create_opt_problem_parameter_list();
        tParameterlist( 0 )( 0 ).set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_hmr_parameter_list();

        tParameterlist( 0 )( 0 ).set( "number_of_elements_per_dimension", tNumElemsPerDim );
        tParameterlist( 0 )( 0 ).set( "domain_dimensions", tDomainDims );
        tParameterlist( 0 )( 0 ).set( "domain_offset", tDomainOffset );
        tParameterlist( 0 )( 0 ).set( "domain_sidesets", tDomainSidesets );
        tParameterlist( 0 )( 0 ).set( "lagrange_output_meshes", std::string( "0" ) );

        tParameterlist( 0 )( 0 ).set( "lagrange_orders", tOrder );
        tParameterlist( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
        tParameterlist( 0 )( 0 ).set( "bspline_orders", tOrder );
        tParameterlist( 0 )( 0 ).set( "bspline_pattern", std::string( "0" ) );

        tParameterlist( 0 )( 0 ).set( "truncate_bsplines", 1 );

        tParameterlist( 0 )( 0 ).set( "use_number_aura", 1 );


        tParameterlist( 0 )( 0 ).set( "initial_refinement", "0" );

        tParameterlist( 0 )( 0 ).set( "initial_refinement_pattern", "0" );

        tParameterlist( 0 )( 0 ).set( "use_multigrid", 0 );
        tParameterlist( 0 )( 0 ).set( "severity_level", 0 );

        tParameterlist( 0 )( 0 ).set( "write_lagrange_output_mesh_to_exodus", "Cantilever_Eigen_HMR.exo" );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_xtk_parameter_list();
        tParameterlist( 0 )( 0 ).set( "decompose", true );
        tParameterlist( 0 )( 0 ).set( "decomposition_type", std::string( "conformal" ) );
        tParameterlist( 0 )( 0 ).set( "enrich", true );
        tParameterlist( 0 )( 0 ).set( "basis_rank", std::string( "bspline" ) );
        tParameterlist( 0 )( 0 ).set( "enrich_mesh_indices", std::string( "0" ) );
        tParameterlist( 0 )( 0 ).set( "ghost_stab", false );
        tParameterlist( 0 )( 0 ).set( "multigrid", false );
        tParameterlist( 0 )( 0 ).set( "verbose", true );
        tParameterlist( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).resize( 1 );

        // Main GEN parameter list
        tParameterlist( 0 )( 0 ) = prm::create_gen_parameter_list();
        tParameterlist( 0 )( 0 ).set( "output_mesh_file", tGENOutputFile );

        // init geometry counter
        uint tGeoCounter = 0;

        // Dummy Geometry
        tParameterlist( 1 ).push_back( prm::create_user_defined_geometry_parameter_list() );
        tParameterlist( 1 )( tGeoCounter ).set( "field_function_name", "Back_Wall" );
        tGeoCounter++;
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Vector< Vector< ParameterList > >& tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        // init property counter
        uint tPropCounter = 0;

        //------------------------------------------------------------------------------
        // MATERIAL PARAMETERS - STRUCTURE (ni-w-alloy?)
        //------------------------------------------------------------------------------

        // Density Shell
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropDensity" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tDensity );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // Youngs Modulus Shell
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", "PropYoungsModulus" );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tYoungsModulus );
        tParameterList( 0 )( tPropCounter ).set( "value_function", "Func_Const" );
        tPropCounter++;

        // Poisson Ratio Shell
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", std::string( "PropPoissonRatio" ) );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tPoissonRatio );
        tParameterList( 0 )( tPropCounter ).set( "value_function", std::string( "Func_Const" ) );
        tPropCounter++;

        //------------------------------------------------------------------------------
        // BOUNDARY CONDITIONS
        //------------------------------------------------------------------------------

        // pressure load
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", std::string( "PropTraction" ) );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tTraction );
        tParameterList( 0 )( tPropCounter ).set( "value_function", std::string( "Func_Neumann_U" ) );
        tPropCounter++;

        // pressure load
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", std::string( "PropPressure" ) );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", tTraction );
        tParameterList( 0 )( tPropCounter ).set( "value_function", std::string( "Func_Const" ) );
        tPropCounter++;

        // Dirichlet structure
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", std::string( "PropDirichletStruct" ) );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", "0.0;0.0" );
        tParameterList( 0 )( tPropCounter ).set( "value_function", std::string( "Func_Const" ) );
        tPropCounter++;

        // time continuity weights
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", std::string( "PropWeightCurrent" ) );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::string( tDensity ) );
        tParameterList( 0 )( tPropCounter ).set( "value_function", std::string( "Func_Const" ) );
        tPropCounter++;

        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", std::string( "PropWeightPrevious" ) );
        tParameterList( 0 )( tPropCounter ).set( "function_parameters", std::string( tDensity ) );
        tParameterList( 0 )( tPropCounter ).set( "value_function", std::string( "Func_Const" ) );
        tPropCounter++;

        // Initial Structure
        tParameterList( 0 ).push_back( prm::create_property_parameter_list() );
        tParameterList( 0 )( tPropCounter ).set( "property_name", std::string( "PropInitialConditionStruct" ) );
        tParameterList( 0 )( tPropCounter ).set( "value_function", std::string( "Func_Initial_Condition_Struct" ) );
        tPropCounter++;

        // init CM counter
        uint tCMCounter = 0;

        //------------------------------------------------------------------------------
        // LINEAR ELASTICITY
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // linear elasticity
            tParameterList( 1 ).push_back( prm::create_constitutive_model_parameter_list() );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_name", "CMStrucLinIso" );
            tParameterList( 1 )( tCMCounter ).set( "model_type", static_cast< uint >( fem::Model_Type::PLANE_STRESS ) );
            tParameterList( 1 )( tCMCounter ).set( "constitutive_type", static_cast< uint >( fem::Constitutive_Type::STRUC_LIN_ISO ) );
            tParameterList( 1 )( tCMCounter ).set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
            tParameterList( 1 )( tCMCounter ).set( "properties", std::string( "PropYoungsModulus,  YoungsModulus;" ) + std::string( "PropPoissonRatio,   PoissonRatio" ) );
            tCMCounter++;
        }

        // init SP counter
        uint tSPCounter = 0;

        //------------------------------------------------------------------------------
        // NITSCHE DIRICHLET
        //------------------------------------------------------------------------------

        if ( tHaveStruct )
        {
            // Displacements - Shell - back wall
            tParameterList( 2 ).push_back( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_name", "SPNitscheStruc" );
            tParameterList( 2 )( tSPCounter ).set( "stabilization_type", static_cast< uint >( fem::Stabilization_Type::DIRICHLET_NITSCHE ) );
            tParameterList( 2 )( tSPCounter ).set( "function_parameters", "100.0" );
            tParameterList( 2 )( tSPCounter ).set( "leader_properties", "PropYoungsModulus,Material" );
            tSPCounter++;
        }

        // init IWG counter
        uint tIWGCounter = 0;

        //------------------------------------------------------------------------------
        // BULK IWGs
        //------------------------------------------------------------------------------

        // linear elasticity
        if ( tHaveStruct )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGStructShell" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_BULK ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );

            {
                tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            }

            tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            // tParameterList( 3 )( tIWGCounter ).set( "leader_properties",          "PropBedding,Bedding" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulk );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // NEUMANN BCs - IWGs
        //------------------------------------------------------------------------------

        // pressure pulling on outside of Shell
        if ( tHaveStruct )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGNeumannPressure" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_NEUMANN ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );

            {
                tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            }

            tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropTraction,Traction" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tLeft );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // DIRICHLET BCS - IWGs
        //------------------------------------------------------------------------------

        // displacements - shell - back wall
        if ( tHaveStruct )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", "IWGDirichletStruct" );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", "UX,UY" );

            {
                tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            }

            tParameterList( 3 )( tIWGCounter ).set( "leader_properties", "PropDirichletStruct,Dirichlet" );
            tParameterList( 3 )( tIWGCounter ).set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
            tParameterList( 3 )( tIWGCounter ).set( "stabilization_parameters", "SPNitscheStruc,DirichletNitsche" );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBottom );
            tIWGCounter++;
        }

        //------------------------------------------------------------------------------
        // IWGs - TIME CONTINUITY
        //------------------------------------------------------------------------------

        // Time continuity Structure
        if ( tHaveStruct )
        {
            tParameterList( 3 ).push_back( prm::create_IWG_parameter_list() );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_name", std::string( "IWGTimeContinuityStruct" ) );
            tParameterList( 3 )( tIWGCounter ).set( "IWG_type", static_cast< uint >( fem::IWG_Type::TIME_CONTINUITY_DOF ) );
            tParameterList( 3 )( tIWGCounter ).set( "dof_residual", std::string( "UX,UY" ) );

            {
                tParameterList( 3 )( tIWGCounter ).set( "leader_dof_dependencies", "UX,UY" );
            }

            tParameterList( 3 )( tIWGCounter ).set( "leader_properties", std::string( "PropWeightCurrent,WeightCurrent;" ) + std::string( "PropWeightPrevious,WeightPrevious;" ) + std::string( "PropInitialConditionStruct,InitialCondition" ) );
            tParameterList( 3 )( tIWGCounter ).set( "mesh_set_names", tBulk );
            tParameterList( 3 )( tIWGCounter ).set( "time_continuity", tTimeContinuity );
            tIWGCounter++;
        }

        // init IQI counter
        uint tIQICounter = 0;

        //         Volume IQI - TotalDomain - use once to find total volume to compute max dof
        tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
        tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQITotalVolume" );
        tParameterList( 4 )( tIQICounter ).set( "IQI_type", static_cast< uint >( fem::IQI_Type::VOLUME ) );

        if ( tHaveStruct )
        {
            tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
        }
        tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulk );
        tIQICounter++;

        if ( tHaveStruct )
        {
            // X-displacement
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPX" );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::EIGEN_VECTOR ) );
            tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1" );
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
            tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 0 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulk );
            tIQICounter++;

            // Y-displacement
            tParameterList( 4 ).push_back( prm::create_IQI_parameter_list() );
            tParameterList( 4 )( tIQICounter ).set( "IQI_name", "IQIBulkDISPY" );
            tParameterList( 4 )( tIQICounter ).set( "IQI_type", (uint)( fem::IQI_Type::EIGEN_VECTOR ) );
            tParameterList( 4 )( tIQICounter ).set( "function_parameters", "1" );
            tParameterList( 4 )( tIQICounter ).set( "dof_quantity", "UX,UY" );
            tParameterList( 4 )( tIQICounter ).set( "leader_dof_dependencies", "UX,UY" );
            tParameterList( 4 )( tIQICounter ).set( "vectorial_field_index", 1 );
            tParameterList( 4 )( tIQICounter ).set( "mesh_set_names", tBulk );
            tIQICounter++;
        }

        // create computation parameter list
        tParameterList( 5 ).resize( 1 );
        tParameterList( 5 )( 0 ) = prm::create_computation_parameter_list();
    }

    void
    SOLParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 8 );
        for ( uint Ik = 0; Ik < 8; Ik++ )
        {
            tParameterlist( Ik ).resize( 1 );
        }

        tParameterlist( 0 )( 0 ) = moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::EIGEN_SOLVER );
        tParameterlist( 0 )( 0 ).set( "Eigen_Algorithm", "EIGALG_BLOCK_DAVIDSON" );
        tParameterlist( 0 )( 0 ).set( "Verbosity", false );
        tParameterlist( 0 )( 0 ).set( "Which", "SM" );
        tParameterlist( 0 )( 0 ).set( "Block_Size", 5 );          // Block Size should be same as Number of Eigen values
        tParameterlist( 0 )( 0 ).set( "NumFreeDofs", 1000 );      // For 2D problem of rectangular elements number of free dofs = 2*node_x*node_y
        tParameterlist( 0 )( 0 ).set( "Num_Eig_Vals", 5 );        // Number of Eigen values should be same as Block Size
        tParameterlist( 0 )( 0 ).set( "Num_Blocks", 2 );          // Number of Blocks should satisfy : Num_Blocks*Block_Size < InitVec Length
        tParameterlist( 0 )( 0 ).set( "MaxSubSpaceDims", 75 );    // Max Subspace Dimension = 3*Block_Size*Num_Eig_Vals
        tParameterlist( 0 )( 0 ).set( "Initial_Guess", 0 );
        tParameterlist( 0 )( 0 ).set( "MaxRestarts", 20 );
        tParameterlist( 0 )( 0 ).set( "Convergence_Tolerance", 1e-05 );
        tParameterlist( 0 )( 0 ).set( "Relative_Convergence_Tolerance", true );
        tParameterlist( 0 )( 0 ).set( "preconditioners", "0" );

        // Ml Preconditioner parameters
        // tParameterlist( 0 )( 0 ).set( "ml_prec_type", "NSSA"); // options: SA, NSSA, DD, DD-ML

        // Print eigenvector parameter
        //tParameterlist( 0 )( 0 ).set( "Print_vector", "LINSOL_EXPORT_MATLAB" );

        tParameterlist( 1 )( 0 ) = moris::prm::create_linear_solver_parameter_list();
        tParameterlist( 1 )( 0 ).set( "DLA_Linear_solver_algorithms", "0" );
        tParameterlist( 1 )( 0 ).set( "RHS_Matrix_Type", "MassMat" );    // MassMat or IdentityMat

        tParameterlist( 2 )( 0 ) = moris::prm::create_nonlinear_algorithm_parameter_list();
        tParameterlist( 2 )( 0 ).set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        tParameterlist( 2 )( 0 ).set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        tParameterlist( 2 )( 0 ).set( "NLA_max_iter", tNLA_max_iter );
        tParameterlist( 2 )( 0 ).set( "NLA_combined_res_jac_assembly", false );

        tParameterlist( 3 )( 0 ) = moris::prm::create_nonlinear_solver_parameter_list();
        tParameterlist( 3 )( 0 ).set( "NLA_DofTypes", "UX,UY" );


        tParameterlist( 4 )( 0 ) = moris::prm::create_time_solver_algorithm_parameter_list();
        tParameterlist( 4 )( 0 ).set( "TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        tParameterlist( 4 )( 0 ).set( "TSA_Time_Frame", tTSA_Time_Frame );

        tParameterlist( 5 )( 0 ) = moris::prm::create_time_solver_parameter_list();
        tParameterlist( 5 )( 0 ).set( "TSA_DofTypes", "UX,UY" );
        tParameterlist( 5 )( 0 ).set( "TSA_Initialize_Sol_Vec", "UX,0.0;UY,0.0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Indices", "0" );
        tParameterlist( 5 )( 0 ).set( "TSA_Output_Criteria", "Output_Criterion" );
        tParameterlist( 5 )( 0 ).set( "TSA_time_level_per_type", "UX,1;UY,1" );


        tParameterlist( 6 )( 0 ) = moris::prm::create_solver_warehouse_parameterlist();
        tParameterlist( 6 )( 0 ).set( "SOL_save_operator_to_matlab", "MassMat" );

        tParameterlist( 7 )( 0 ) = moris::prm::create_preconditioner_parameter_list(sol::PreconditionerType::IFPACK);
        // Ifpack Preconditioner parameters
        tParameterlist( 7 )( 0 ).set( "ifpack_prec_type", "Amesos" );
        tParameterlist( 7 )( 0 ).set( "amesos: solver type", gPrecSolver );    // Amesos_Umfpack or Amesos_Pardiso
        // Preconditioner parameters
        tParameterlist( 7 )( 0 ).set( "overlap-level", 0 );
        tParameterlist( 7 )( 0 ).set( "schwarz: combine mode", "add" );    // for Amesos_Umfpack and Amesos_Pardiso provide this parameter with "add" mode
    }

    void
    MSIParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_msi_parameter_list();
        tParameterlist( 0 )( 0 ).set( "number_eigen_vectors", 5 );
    }

    void
    VISParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).resize( 1 );

        tParameterlist( 0 )( 0 ) = prm::create_vis_parameter_list();
        tParameterlist( 0 )( 0 ).set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 )( 0 ).set( "Mesh_Type", static_cast< uint >( vis::VIS_Mesh_Type::STANDARD ) );
        tParameterlist( 0 )( 0 ).set( "Set_Names", tBulk );

        if ( tHaveStruct )
        {
            tParameterlist( 0 )( 0 ).set( "Field_Names", "UX,UY" );
            tParameterlist( 0 )( 0 ).set( "Field_Type", "NODAL,NODAL" );
            tParameterlist( 0 )( 0 ).set( "IQI_Names", "IQIBulkDISPX,IQIBulkDISPY" );
        }

        tParameterlist( 0 )( 0 ).set( "Save_Frequency", 1 );
    }

    void
    MORISGENERALParameterList( Vector< Vector< ParameterList > >& tParameterlist )
    {
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
