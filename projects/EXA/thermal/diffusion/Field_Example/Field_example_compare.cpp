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
#include "parameters.hpp"
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

    void OPTParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false);
    }

    void HMRParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", 20, 20 );
        aParameterLists.set( "domain_dimensions",                2.0, 2.0 );
        aParameterLists.set( "lagrange_output_meshes",           "0");

        aParameterLists.set( "lagrange_orders",  "2" );
        aParameterLists.set( "lagrange_pattern",  "0" )  ;
        aParameterLists.set( "bspline_orders",   "2" );
        aParameterLists.set( "bspline_pattern",   "1" )  ;

        aParameterLists.set( "lagrange_to_bspline", "0") ;

        aParameterLists.set( "truncate_bsplines",  1 );
        aParameterLists.set( "refinement_buffer",  1 );
        aParameterLists.set( "staircase_buffer",   1 );
        aParameterLists.set( "initial_refinement", "1,0" );
        aParameterLists.set( "initial_refinement_pattern", "0,1" );

        aParameterLists.set( "use_number_aura", 1);

        aParameterLists.set( "use_multigrid",  0 );
        aParameterLists.set( "severity_level", 0 );
    }

    void XTKParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists.set( "decompose",                 true );
        aParameterLists.set( "decomposition_type",        "conformal") ;
        aParameterLists.set( "enrich",                    true );
        aParameterLists.set( "basis_rank",                "bspline") ;
        aParameterLists.set( "enrich_mesh_indices",       "0") ;
        aParameterLists.set( "ghost_stab",                tUseGhost );
        aParameterLists.set( "multigrid",                 false );
        aParameterLists.set( "verbose",                   true );
        aParameterLists.set( "print_enriched_ig_mesh",    false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    void GENParameterList( Module_Parameter_Lists & aParameterLists )
    {
        // Geometry parameter lists
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name",       "LevelSetFunction" );
        aParameterLists.set("number_of_refinements", 0);
        aParameterLists.set("refinement_mesh_index", 0);
        aParameterLists.set("discretization_mesh_index", 0);
        aParameterLists.set("use_multilinear_interpolation", true);
    }

    void FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();
        // create a cell of cell of parameter list for fem

        //------------------------------------------------------------------------------

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name",            "PropConductivity");
        aParameterLists.set( "function_parameters",      tCond);
        aParameterLists.set( "value_function",           "Func_Const");

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name",            "PropTemperature");
        aParameterLists.set( "function_parameters",      "1.0" );
        aParameterLists.set( "value_function",           "Func_Const");

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name",            "PropSurfaceFlux") ;
        aParameterLists.set( "function_parameters",      "1.0" );
        aParameterLists.set( "value_function",           "Func_Const") ;

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name",            "PropField") ;
        aParameterLists.set( "value_function",           "Func_Field") ;
        aParameterLists.set( "field_dependencies",       "FIELD_1" );

        //------------------------------------------------------------------------------

        // create parameter list for constitutive model - Inclusion
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion") ;
        aParameterLists.set( "constitutive_type",  fem::Constitutive_Type::DIFF_LIN_ISO ) ;
        aParameterLists.set( "dof_dependencies",  std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity , Conductivity");

        //------------------------------------------------------------------------------

        // Nitsche stabilization parameter for thermal problem
        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name",      "SPNitscheTemp") ;
        aParameterLists.set( "stabilization_type",       fem::Stabilization_Type::DIRICHLET_NITSCHE ) ;
        aParameterLists.set( "function_parameters",     "100.0") ;
        aParameterLists.set( "leader_properties",       "PropConductivity,Material") ;

        if (tUseGhost)
        {
            // Ghost stabilization parameter for thermal problem
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name",      "SPGhostTemp") ;
            aParameterLists.set( "stabilization_type",       fem::Stabilization_Type::GHOST_DISPL ) ;
            aParameterLists.set( "function_parameters",     "0.01") ;
            aParameterLists.set( "leader_properties",       "PropConductivity,Material") ;
            }

        //------------------------------------------------------------------------------
        // create IWG - bulk diffusion
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name",                   "IWGBulkTemp") ;
        aParameterLists.set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_BULK ) ;
        aParameterLists.set( "dof_residual",               "TEMP") ;
        aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion") ;
        aParameterLists.set( "mesh_set_names",             tPhase1 );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name",                   "IWGInletFlux") ;
        aParameterLists.set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_NEUMANN ) ;
        aParameterLists.set( "dof_residual",               "TEMP") ;
        aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists.set( "leader_properties",          "PropSurfaceFlux,Neumann") ;
        aParameterLists.set( "mesh_set_names",             "iside_b0_1_b1_0" );

         // create IWG - Dirichlet temp
        aParameterLists( FEM::IWG ).add_parameter_list();
         aParameterLists.set( "IWG_name",                   "IWGDirichletTemp") ;
         aParameterLists.set( "IWG_type",                    fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE ) ;
         aParameterLists.set( "dof_residual",               "TEMP") ;
         aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
         aParameterLists.set( "leader_properties",          "PropTemperature,Dirichlet") ;
         aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion") ;
         aParameterLists.set( "stabilization_parameters",   "SPNitscheTemp,DirichletNitsche") ;
         aParameterLists.set( "mesh_set_names",             tRightSurface );

        if (tUseGhost)
        {
            // create IWG - ghost temp
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name",                   "IWGGhostTemp") ;
            aParameterLists.set( "IWG_type",                    fem::IWG_Type::GHOST_NORMAL_FIELD ) ;
            aParameterLists.set( "dof_residual",               "TEMP") ;
            aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
            aParameterLists.set( "follower_dof_dependencies",     "TEMP") ;
            aParameterLists.set( "stabilization_parameters",   "SPGhostTemp,GhostSP") ;
            aParameterLists.set( "mesh_set_names",             tPhase1Ghost );
            }

        //------------------------------------------------------------------------------
        // Nodal Temperature IQI
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name",                   "IQIBulkTEMP") ;
        aParameterLists.set( "IQI_type",                    fem::IQI_Type::DOF ) ;
        aParameterLists.set( "dof_quantity",               "TEMP");
        aParameterLists.set( "leader_dof_dependencies",    "TEMP") ;
        aParameterLists.set( "vectorial_field_index",      0 );
        aParameterLists.set( "mesh_set_names",             tPhase1 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name",                   "IQIBulkL2Error") ;
        aParameterLists.set( "IQI_type",                    fem::IQI_Type::H1_ERROR ) ;
        aParameterLists.set( "dof_quantity",               "TEMP") ;
        aParameterLists.set( "function_parameters",        "1.0/0.0" );
        aParameterLists.set( "leader_properties",          "PropField,L2_Reference" );
        aParameterLists.set( "mesh_set_names",             tPhase1 );

        //------------------------------------------------------------------------------
        // fill the computation part of the parameter list
        aParameterLists( FEM::COMPUTATION );

        //------------------------------------------------------------------------------

        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name",                "FieldNodalTEMP") ;
        aParameterLists.set( "field_entity_type",         "NODAL" );
        aParameterLists.set( "field_type",                "FIELD_1" );
        aParameterLists.set( "field_create_from_file",    tFieldRefPath) ;
        //aParameterLists.set( "IQI_Name",                  "IQIBulkTEMP") ;
        //aParameterLists.set( "field_output_to_file",      "Field_example_write.hdf5" );

    }

    void SOLParameterList( Module_Parameter_Lists & aParameterLists )
    {

        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set("NLA_rel_res_norm_drop",    tNLA_rel_res_norm_drop );
        aParameterLists.set("NLA_relaxation_parameter", tNLA_relaxation_parameter  );
        aParameterLists.set("NLA_max_iter",             tNLA_max_iter );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set("NLA_DofTypes"      , "TEMP") ;

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set("TSA_Num_Time_Steps", tTSA_Num_Time_Steps );
        aParameterLists.set("TSA_Time_Frame",     tTSA_Time_Frame );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set("TSA_DofTypes",           "TEMP") ;
        aParameterLists.set("TSA_Initialize_Sol_Vec", "TEMP,0.0");
        aParameterLists.set("TSA_Output_Indices",     "0") ;
        aParameterLists.set("TSA_Output_Criteria",     "Output_Criterion") ;

        aParameterLists( SOL::SOLVER_WAREHOUSE );
        //aParameterLists.set( "SOL_save_final_sol_vec_to_file", tSolVecRefPath);
        aParameterLists.set( "SOL_load_sol_vec_from_file", tSolVecRefPath);

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list(  sol::PreconditionerType::NONE );
    }

    void MSIParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists.set("order_adofs_by_host",false);
    }

    void VISParameterList( Module_Parameter_Lists & aParameterLists )
    {
        aParameterLists.set( "File_Name"  , std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type"  ,  vis::VIS_Mesh_Type::STANDARD ) ;
        aParameterLists.set( "Set_Names"  , tPhase1 );
        aParameterLists.set( "Field_Names", "TEMP,DIFFERENCE") ;
        aParameterLists.set( "Field_Type" , "NODAL,NODAL") ;
        aParameterLists.set( "IQI_Names"  , "IQIBulkTEMP,IQIBulkL2Error") ;
        aParameterLists.set( "Save_Frequency", 1 );
    }

    void MORISGENERALParameterList( Module_Parameter_Lists & aParameterLists )
    {

    }

    /* ------------------------------------------------------------------------ */
}

//------------------------------------------------------------------------------
#ifdef  __cplusplus
}
#endif

