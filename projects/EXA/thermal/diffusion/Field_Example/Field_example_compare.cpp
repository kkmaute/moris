/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Field_example_compare.cpp
 *
 */

#include <string>
#include <iostream>
#include "paths.hpp"
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
#include "fn_norm.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

#ifdef  __cplusplus
extern "C"
{
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // Mesh Set Information

    // Phase 1: back  - Material 1
    // Phase 2: front - Material 2

    std::string tPhase1        = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tPhase2        = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    std::string tInterface     = "dbl_iside_p0_1_p1_0";

    std::string tLeftSurface   = "SideSet_4_n_p1,SideSet_4_c_p1";

    std::string tRightSurface   = "SideSet_2_n_p1";

    std::string tButtomSurface   = "SideSet_1_n_p1,SideSet_1_c_p1";

    std::string tPhase1Ghost   = "ghost_p1";
    std::string tPhase2Ghost   = "ghost_p0";

    std::string tTotalDomain   = tPhase1 + "," + tPhase2;

    /* ------------------------------------------------------------------------ */
    // geometry parameters

    /* ------------------------------------------------------------------------ */
    // boundary condition

    /* ------------------------------------------------------------------------ */
    // material parameters

    std::string  tCond = "1.0";

    /* ------------------------------------------------------------------------ */
    // Solver config

    moris::real tNLA_rel_res_norm_drop = 1.0e-08;
    moris::real tNLA_relaxation_parameter = 1.0;
    int tNLA_max_iter = 2;

    int tTSA_Num_Time_Steps = 1;
    moris::real tTSA_Time_Frame = 1.0e0;

    /* ------------------------------------------------------------------------ */
    // Minimum level set value
    bool tUseGhost = true;

    /* ------------------------------------------------------------------------ */
    // Output Config

    std::string tOutputFileName = "Field_example_compare.exo";

    std::string tPrefix = moris::get_base_moris_dir();
    std::string tFieldRefPath = tPrefix + "/projects/EXA/thermal/diffusion/Field_Example/Field_example_ref.hdf5";
    std::string tSolVecRefPath = tPrefix + "/projects/EXA/thermal/diffusion/Field_Example/Solution_Vector_Ref.hdf5";

    moris::real LevelSetFunction(
            const moris::Matrix< DDRMat >     & aCoordinates,
            const Vector< real > & aGeometryParameters )
    {
        //return norm(aCoordinates) - 0.01;
        return std::pow( aCoordinates( 0 ), 2 )/1 + std::pow( aCoordinates( 1 ), 2 )/2 - 0.331;
    }

    /* ------------------------------------------------------------------------ */
    // Constant function for properties

    void Func_Const(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    void Func_Field(
            moris::Matrix< moris::DDRMat >                 & aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
            moris::fem::Field_Interpolator_Manager         * aFIManager )
    {
        aPropMatrix = aFIManager->
                get_field_interpolators_for_type( mtk::Field_Type::FIELD_1 )->
                val();
    }

    /* ------------------------------------------------------------------------ */

    bool Output_Criterion( moris::tsa::Time_Solver * aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    void OPTParameterList( Vector< Submodule_Parameter_Lists > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_opt_problem_parameter_list() );

        tParameterlist( 0 ).set( "is_optimization_problem", false);
    }

    void HMRParameterList( Vector< Submodule_Parameter_Lists > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );

        tParameterlist( 0 ).set( "number_of_elements_per_dimension", "20,   20" );
        tParameterlist( 0 ).set( "domain_dimensions",                "2.0,   2.0" );
        tParameterlist( 0 ).set( "domain_offset",                    "0.0,  0.0"  );
        tParameterlist( 0 ).set( "domain_sidesets",                  "1,2,3,4");
        tParameterlist( 0 ).set( "lagrange_output_meshes",           "0");

        tParameterlist( 0 ).set( "lagrange_orders",  "2" );
        tParameterlist( 0 ).set( "lagrange_pattern",  "0" )  ;
        tParameterlist( 0 ).set( "bspline_orders",   "2" );
        tParameterlist( 0 ).set( "bspline_pattern",   "1" )  ;

        tParameterlist( 0 ).set( "lagrange_to_bspline", "0") ;

        tParameterlist( 0 ).set( "truncate_bsplines",  1 );
        tParameterlist( 0 ).set( "refinement_buffer",  1 );
        tParameterlist( 0 ).set( "staircase_buffer",   1 );
        tParameterlist( 0 ).set( "initial_refinement", "1,0" );
        tParameterlist( 0 ).set( "initial_refinement_pattern", "0,1" );

        tParameterlist( 0 ).set( "use_number_aura", 1);

        tParameterlist( 0 ).set( "use_multigrid",  0 );
        tParameterlist( 0 ).set( "severity_level", 0 );
    }

    void XTKParameterList( Vector< Submodule_Parameter_Lists > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
        tParameterlist( 0 ).set( "decompose",                 true );
        tParameterlist( 0 ).set( "decomposition_type",        "conformal") ;
        tParameterlist( 0 ).set( "enrich",                    true );
        tParameterlist( 0 ).set( "basis_rank",                "bspline") ;
        tParameterlist( 0 ).set( "enrich_mesh_indices",       "0") ;
        tParameterlist( 0 ).set( "ghost_stab",                tUseGhost );
        tParameterlist( 0 ).set( "multigrid",                 false );
        tParameterlist( 0 ).set( "verbose",                   true );
        tParameterlist( 0 ).set( "print_enriched_ig_mesh",    false );
        tParameterlist( 0 ).set( "exodus_output_XTK_ig_mesh", true );
        tParameterlist( 0 ).set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( Vector< Submodule_Parameter_Lists > & tParameterlist )
    {
        tParameterlist.resize( 3 );
        tParameterlist( 0 ).add_parameter_list( prm::create_gen_parameter_list() );


        // Geometry parameter lists
        tParameterlist( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        tParameterlist( 1 ).set( "field_function_name",       "LevelSetFunction" );
        tParameterlist( 1 ).set("number_of_refinements", 0);
        tParameterlist( 1 ).set("refinement_mesh_index", 0);
        tParameterlist( 1 ).set("discretization_mesh_index", 0);
        tParameterlist( 1 ).set("use_multilinear_interpolation", true);
    }

    void FEMParameterList( Vector< Submodule_Parameter_Lists > & tParameterList )
    {
        // create a cell of cell of parameter list for fem
        tParameterList.resize( 8 );

        //------------------------------------------------------------------------------

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name",            "PropConductivity");
        tParameterList( 0 ).set( "function_parameters",      tCond);
        tParameterList( 0 ).set( "value_function",           "Func_Const");

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name",            "PropTemperature");
        tParameterList( 0 ).set( "function_parameters",      "1.0" );
        tParameterList( 0 ).set( "value_function",           "Func_Const");

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name",            "PropSurfaceFlux") ;
        tParameterList( 0 ).set( "function_parameters",      "1.0" );
        tParameterList( 0 ).set( "value_function",           "Func_Const") ;

        tParameterList( 0 ).add_parameter_list( prm::create_property_parameter_list() );
        tParameterList( 0 ).set( "property_name",            "PropField") ;
        tParameterList( 0 ).set( "value_function",           "Func_Field") ;
        tParameterList( 0 ).set( "field_dependencies",       "FIELD_1" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - Inclusion
        tParameterList( 1 ).add_parameter_list( prm::create_constitutive_model_parameter_list() );
        tParameterList( 1 ).set( "constitutive_name", "CMDiffusion") ;
        tParameterList( 1 ).set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        tParameterList( 1 ).set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        tParameterList( 1 ).set( "properties",
                "PropConductivity , Conductivity");

        //------------------------------------------------------------------------------

        // Nitsche stabilization parameter for thermal problem
        tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
        tParameterList( 2 ).set( "stabilization_name",      "SPNitscheTemp") ;
        tParameterList( 2 ).set( "stabilization_type",       fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        tParameterList( 2 ).set( "function_parameters",     "100.0") ;
        tParameterList( 2 ).set( "leader_properties",       "PropConductivity,Material") ;

        if (tUseGhost)
        {
            // Ghost stabilization parameter for thermal problem
            tParameterList( 2 ).add_parameter_list( prm::create_stabilization_parameter_parameter_list() );
            tParameterList( 2 ).set( "stabilization_name",      "SPGhostTemp") ;
            tParameterList( 2 ).set( "stabilization_type",       fem::Stabilization_Type::GHOST_DISPL ) ;
            tParameterList( 2 ).set( "function_parameters",     "0.01") ;
            tParameterList( 2 ).set( "leader_properties",       "PropConductivity,Material") ;
            }

        //------------------------------------------------------------------------------
        // create IWG - bulk diffusion
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name",                   "IWGBulkTemp") ;
        tParameterList( 3 ).set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_BULK ) ;
        tParameterList( 3 ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion") ;
        tParameterList( 3 ).set( "mesh_set_names",             tPhase1 );

        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
        tParameterList( 3 ).set( "IWG_name",                   "IWGInletFlux") ;
        tParameterList( 3 ).set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        tParameterList( 3 ).set( "dof_residual",               "TEMP") ;
        tParameterList( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 3 ).set( "leader_properties",          "PropSurfaceFlux,Neumann") ;
        tParameterList( 3 ).set( "mesh_set_names",             "iside_b0_1_b1_0" );

         // create IWG - Dirichlet temp
        tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
         tParameterList( 3 ).set( "IWG_name",                   "IWGDirichletTemp") ;
         tParameterList( 3 ).set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
         tParameterList( 3 ).set( "dof_residual",               "TEMP") ;
         tParameterList( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
         tParameterList( 3 ).set( "leader_properties",          "PropTemperature,Dirichlet") ;
         tParameterList( 3 ).set( "leader_constitutive_models", "CMDiffusion,Diffusion") ;
         tParameterList( 3 ).set( "stabilization_parameters",   "SPNitscheTemp,DirichletNitsche") ;
         tParameterList( 3 ).set( "mesh_set_names",             tRightSurface );

        if (tUseGhost)
        {
            // create IWG - ghost temp
            tParameterList( 3 ).add_parameter_list( prm::create_IWG_parameter_list() );
            tParameterList( 3 ).set( "IWG_name",                   "IWGGhostTemp") ;
            tParameterList( 3 ).set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            tParameterList( 3 ).set( "dof_residual",               "TEMP") ;
            tParameterList( 3 ).set( "leader_dof_dependencies",    "TEMP") ;
            tParameterList( 3 ).set( "follower_dof_dependencies",     "TEMP") ;
            tParameterList( 3 ).set( "stabilization_parameters",   "SPGhostTemp,GhostSP") ;
            tParameterList( 3 ).set( "mesh_set_names",             tPhase1Ghost );
            }

        //------------------------------------------------------------------------------
        // Nodal Temperature IQI
        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name",                   "IQIBulkTEMP") ;
        tParameterList( 4 ).set( "IQI_type",                    fem::IQI_Type::DOF ) ;
        tParameterList( 4 ).set( "dof_quantity",               "TEMP");
        tParameterList( 4 ).set( "leader_dof_dependencies",    "TEMP") ;
        tParameterList( 4 ).set( "vectorial_field_index",      0 );
        tParameterList( 4 ).set( "mesh_set_names",             tPhase1 );

        tParameterList( 4 ).add_parameter_list( prm::create_IQI_parameter_list() );
        tParameterList( 4 ).set( "IQI_name",                   "IQIBulkL2Error") ;
        tParameterList( 4 ).set( "IQI_type",                    fem::IQI_Type::H1_ERROR ) ;
        tParameterList( 4 ).set( "dof_quantity",               "TEMP") ;
        tParameterList( 4 ).set( "function_parameters",        "1.0/0.0" );
        tParameterList( 4 ).set( "leader_properties",          "PropField,L2_Reference" );
        tParameterList( 4 ).set( "mesh_set_names",             tPhase1 );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        tParameterList( 5 ).add_parameter_list( prm::create_computation_parameter_list() );

        //------------------------------------------------------------------------------

        tParameterList( 6 ).add_parameter_list( prm::create_fem_field_parameter_list() );
        tParameterList( 6 ).set( "field_name",                "FieldNodalTEMP") ;
        tParameterList( 6 ).set( "field_entity_type",         "NODAL" );
        tParameterList( 6 ).set( "field_type",                "FIELD_1" );
        tParameterList( 6 ).set( "field_create_from_file",    tFieldRefPath) ;
        //tParameterList( 6 ).set( "IQI_Name",                  "IQIBulkTEMP") ;
        //tParameterList( 6 ).set( "field_output_to_file",      "Field_example_write.hdf5" );

    }

    void SOLParameterList( Vector< Submodule_Parameter_Lists > & tParameterlist )
    {
        tParameterlist.resize( 8 );

        tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );

        tParameterlist( 1 ).add_parameter_list( moris::prm::create_linear_solver_parameter_list() );

        tParameterlist( 2 ).add_parameter_list( moris::prm::create_nonlinear_algorithm_parameter_list() );
        tParameterlist( 2 ).set("NLA_rel_res_norm_drop",    tNLA_rel_res_norm_drop );
        tParameterlist( 2 ).set("NLA_relaxation_parameter", tNLA_relaxation_parameter  );
        tParameterlist( 2 ).set("NLA_max_iter",             tNLA_max_iter );

        tParameterlist( 3 ).add_parameter_list( moris::prm::create_nonlinear_solver_parameter_list() );
        tParameterlist( 3 ).set("NLA_DofTypes"      , "TEMP") ;

        tParameterlist( 4 ).add_parameter_list( moris::prm::create_time_solver_algorithm_parameter_list() );
        tParameterlist( 4 ).set("TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        tParameterlist( 4 ).set("TSA_Time_Frame",     tTSA_Time_Frame );

        tParameterlist( 5 ).add_parameter_list( moris::prm::create_time_solver_parameter_list() );
        tParameterlist( 5 ).set("TSA_DofTypes",           "TEMP") ;
        tParameterlist( 5 ).set("TSA_Initialize_Sol_Vec", "TEMP,0.0");
        tParameterlist( 5 ).set("TSA_Output_Indices",     "0") ;
        tParameterlist( 5 ).set("TSA_Output_Criteria",     "Output_Criterion") ;

        tParameterlist( 6 ).add_parameter_list( moris::prm::create_solver_warehouse_parameterlist() );
        //tParameterlist( 6 ).set( "SOL_save_final_sol_vec_to_file", tSolVecRefPath);
        tParameterlist( 6 ).set( "SOL_load_sol_vec_from_file", tSolVecRefPath);

        tParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE ) );
    }

    void MSIParameterList( Vector< Submodule_Parameter_Lists > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_msi_parameter_list() );
        tParameterlist( 0 ).set("order_adofs_by_host",false);
    }

    void VISParameterList( Vector< Submodule_Parameter_Lists > & tParameterlist )
    {
        tParameterlist.resize( 1 );
        tParameterlist( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        tParameterlist( 0 ).set( "File_Name"  , std::pair< std::string, std::string >( "./", tOutputFileName ) );
        tParameterlist( 0 ).set( "Mesh_Type"  ,  vis::VIS_Mesh_Type::STANDARD ) ;
        tParameterlist( 0 ).set( "Set_Names"  , tPhase1 );
        tParameterlist( 0 ).set( "Field_Names", "TEMP,DIFFERENCE") ;
        tParameterlist( 0 ).set( "Field_Type" , "NODAL,NODAL") ;
        tParameterlist( 0 ).set( "IQI_Names"  , "IQIBulkTEMP,IQIBulkL2Error") ;
        tParameterlist( 0 ).set( "Save_Frequency", 1 );
    }

    void MORISGENERALParameterList( Vector< Submodule_Parameter_Lists > & tParameterlist )
    {

    }

    /* ------------------------------------------------------------------------ */
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

